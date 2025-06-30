import sys
from collections import defaultdict
from itertools import product
from Bio import SeqIO
import argparse

# =========================
# Utility functions
# =========================
def revcomp(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]

def canonical_kmer(seq):
    rc = revcomp(seq)
    return min(seq, rc)

def hamming_neighbors(kmer):
    """Generate all k-mers with Hamming distance 1."""
    neighbors = set()
    bases = 'ACGT'
    for i in range(len(kmer)):
        for b in bases:
            if b != kmer[i]:
                neighbors.add(canonical_kmer(kmer[:i] + b + kmer[i+1:]))
    return neighbors

def extract_kmers_from_fasta(fasta_path, k, methylation=False):
    kmers = set()
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        if methylation:
            seq = seq.replace('C', 'T')  # Convert to bisulfite-treated form
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if 'N' not in kmer:
                kmers.add(canonical_kmer(kmer))
    return kmers

# =========================
# Main classification
# =========================
def classify_kmers(host_kmers, graft_kmers):
    all_kmers = host_kmers | graft_kmers
    classification = defaultdict(int)

    for kmer in all_kmers:
        in_host = kmer in host_kmers
        in_graft = kmer in graft_kmers

        if in_host and in_graft:
            classification['both'] += 1
        elif in_host:
            classification['host'] += 1
        elif in_graft:
            classification['graft'] += 1

    # Classify weak kmers
    weak_host = 0
    weak_graft = 0
    print("Scanning for weak k-mers... (this may take a while)")
    for kmer in host_kmers - graft_kmers:
        for neighbor in hamming_neighbors(kmer):
            if neighbor in graft_kmers:
                weak_host += 1
                break

    for kmer in graft_kmers - host_kmers:
        for neighbor in hamming_neighbors(kmer):
            if neighbor in host_kmers:
                weak_graft += 1
                break

    classification['weak_host'] = weak_host
    classification['weak_graft'] = weak_graft

    return classification, len(all_kmers)

# =========================
# Argument parsing and driver
# =========================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate k-mer classification summary like Xengsort Table 2.")
    parser.add_argument("--host", required=True, help="Path to host (e.g. mouse) genome FASTA.")
    parser.add_argument("--graft", required=True, help="Path to graft (e.g. human) genome FASTA.")
    parser.add_argument("-k", type=int, default=25, help="k-mer size (default: 25)")
    parser.add_argument("--methylation", action="store_true", help="Convert all Cs to Ts in reference (simulate bisulfite-treated 3-letter genome).")
    args = parser.parse_args()

    print(f"Extracting {args.k}-mers from host (methylation mode: {args.methylation})...")
    host_kmers = extract_kmers_from_fasta(args.host, args.k, methylation=args.methylation)
    print(f"Found {len(host_kmers):,} host kmers.")

    print(f"Extracting {args.k}-mers from graft (methylation mode: {args.methylation})...")
    graft_kmers = extract_kmers_from_fasta(args.graft, args.k, methylation=args.methylation)
    print(f"Found {len(graft_kmers):,} graft kmers.")

    result, total = classify_kmers(host_kmers, graft_kmers)

    print("\n--- K-mer Classification Summary ---")
    print(f"Total k-mers      : {total:,}")
    for label in ['host', 'graft', 'both', 'weak_host', 'weak_graft']:
        count = result[label]
        pct = count / total * 100
        print(f"{label:<15}: {count:,} ({pct:.2f}%)")
