#!/usr/bin/env python3

import sys
import gzip

def open_fastq(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def get_read_ids(fastq_file, genome_tag):
    """
    Extract read IDs from a FASTQ file that contain the specified genome tag.
    """
    ids = set()
    ids_wrong = set()
    with open_fastq(fastq_file) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            # Remove '@' and any trailing whitespace
            read_id = header.strip().split()[0][1:]
            if genome_tag in read_id:
                ids.add(read_id)
            else: 
                ids_wrong.add(read_id)
    return ids, ids_wrong

def main():
    if len(sys.argv) != 4:
        print("Usage: python stat_accuracy.py <true_total> <host_R1.fastq.gz>,<host_tag> <graft_R1.fastq.gz>,<graft_tag>")
        sys.exit(1)

    true_total = int(sys.argv[1])
    host_arg = sys.argv[2]
    graft_arg = sys.argv[3]

    host_fastq, host_tag = host_arg.split(',')
    graft_fastq, graft_tag = graft_arg.split(',')

    # Get true host and graft read IDs
    true_host_ids, true_host_ids_wrong = get_read_ids(host_fastq, host_tag)
    true_graft_ids, true_graft_ids_wrong = get_read_ids(graft_fastq, graft_tag)

    # Calculate stats
    n_host = len(true_host_ids)
    n_graft = len(true_graft_ids)
    n_total = n_host + n_graft

    print(f"Total true reads: {true_total}")
    print(f"Host ({host_tag}) reads found: {n_host}")
    print(f"Graft ({graft_tag}) reads found: {n_graft}")
    print(f"Sum found: {n_total}")
    print(f"Fraction host: {n_host/true_total:.4f}")
    print(f"Fraction graft: {n_graft/true_total:.4f}")
    print(f"Fraction host (wrong IDs): {len(true_host_ids_wrong)/n_host:.4f}")
    print(f"Fraction graft (wrong IDs): {len(true_graft_ids_wrong)/n_graft:.4f}")

    # Sensitivity: correctly assigned / all true
    sensitivity_host = n_host / true_total if true_total else 0
    sensitivity_graft = n_graft / true_total if true_total else 0

    # Specificity: correctly rejected / all negatives
    specificity_host = 1 - (len(true_graft_ids_wrong) / n_graft) if n_graft else 0
    specificity_graft = 1 - (len(true_host_ids_wrong) / n_host) if n_host else 0

    print(f"Sensitivity (host): {sensitivity_host:.4f}")
    print(f"Specificity (host): {specificity_host:.4f}")
    print(f"Sensitivity (graft): {sensitivity_graft:.4f}")
    print(f"Specificity (graft): {specificity_graft:.4f}")

if __name__ == "__main__":
    main()