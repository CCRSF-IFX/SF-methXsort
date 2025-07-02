#!/usr/bin/env python3

import toolshed
from toolshed import nopen, reader, is_newer_b
import argparse
import sys
import os
from itertools import groupby
import gzip
import pysam
from multiprocessing import Pool, cpu_count
import subprocess
import shlex
import pathlib

def wrap(seq, width=60):
    return [seq[i:i+width] for i in range(0, len(seq), width)]

def fasta_iter(fasta_name):
    fh = nopen(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        yield header, "".join(s.strip() for s in next(faiter)).upper()
        

def convert_fasta(ref_fasta, out_fa=None):
    """
    Code originally copied from: bwa-meth and adapted for use in methxsort.py
    Convert a fasta file to a bwameth compatible format.
    """
    if out_fa is None:
        out_fa = ref_fasta + ".converted"
    msg = "c2t/g2a in %s to %s" % (ref_fasta, out_fa)
    if is_newer_b(ref_fasta, out_fa):
        sys.stderr.write("already converted: %s\n" % msg)
        return out_fa
    sys.stderr.write("converting %s\n" % msg)
    try:
        fh = open(out_fa, "w")
        for header, seq in fasta_iter(ref_fasta):
            ########### Reverse ######################
            fh.write(">r%s\n" % header)
            for line in wrap(seq.replace("G", "A")):
                fh.write(line + '\n')

            ########### Forward ######################
            fh.write(">f%s\n" % header)
            for line in wrap(seq.replace("C", "T")):
                fh.write(line + '\n')
        fh.close()
    except:
        try:
            fh.close()
        except UnboundLocalError:
            pass
        os.unlink(out_fa)
        raise
    return out_fa


def convert_reads_c2t_r1_g2a_r2(read, read2=None, out=None, out2=None):
    """
    Fast conversion using awk+sed+gzip for large FASTQ files.
    C->T in read (single or read1), G->A in read2 (if paired-end).
    Store the original sequence in the header line as an extra field.
    """
    if out is None:
        out = read + ".meth"
    if read2 and out2 is None:
        out2 = read2 + ".meth"

    # For read1: C->T, add original sequence to header (corrected awk command)
    cmd1 = (
        f"zcat {read} | "
        "awk 'NR%4==1{{getline seq; print $0 \" ORIGINAL_SEQ:\" seq; print seq; getline; print $0; getline; print $0}}' | "
        "sed '2~4s/C/T/g;2~4s/c/t/g' | gzip > {out}"
    ).format(read=read, out=out)
    print(f"[convert-reads] CMD1: {cmd1}", file=sys.stdout)

    # For read2: G->A, add original sequence to header (corrected awk command)
    if read2 and out2:
        cmd2 = (
            f"zcat {read2} | "
            "awk 'NR%4==1{{getline seq; print $0 \" ORIGINAL_SEQ:\" seq; print seq; getline; print $0; getline; print $0}}' | "
            "sed '2~4s/G/A/g;2~4s/g/a/g' | gzip > {out2}"
        ).format(read2=read2, out2=out2)
        print(f"[convert-reads] CMD2: {cmd2}", file=sys.stdout)
        # Run both commands in parallel
        p1 = subprocess.Popen(cmd1, shell=True)
        p2 = subprocess.Popen(cmd2, shell=True)
        p1.wait()
        p2.wait()
        if p1.returncode != 0 or p2.returncode != 0:
            raise subprocess.CalledProcessError(p1.returncode if p1.returncode != 0 else p2.returncode, cmd1 if p1.returncode != 0 else cmd2)
        return out, out2
    else:
        subprocess.check_call(cmd1, shell=True)
        return out

# Old Python implementation (commented out for reference)
# def convert_reads_c2t_r1_g2a_r2(read, read2=None, out=None, out2=None):
#     """
#     Convert C->T in read (single or read1) and G->A in read2 (if paired-end).
#     Store the original sequence in the header line as an extra field.
#     """
#     def convert_r1(seq):
#         return seq.replace("C", "T").replace("c", "t")
#     def convert_r2(seq):
#         return seq.replace("G", "A").replace("g", "a")
#     if out is None:
#         out = read + ".meth"
#     if read2 and out2 is None:
#         out2 = read2 + ".meth"

#     def open_out(filename):
#         if filename and filename.endswith('.gz'):
#             return gzip.open(filename, "wt")
#         else:
#             return open(filename, "w")

#     with nopen(read) as r1, open_out(out) as o1:
#         while True:
#             header = r1.readline()
#             if not header:
#                 break
#             seq = r1.readline()
#             plus = r1.readline()
#             qual = r1.readline()
#             # Add original sequence to header
#             header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
#             o1.write(header)
#             o1.write(convert_r1(seq))
#             o1.write(plus)
#             o1.write(qual)
#     if read2:
#         with nopen(read2) as r2, open_out(out2) as o2:
#             while True:
#                 header = r2.readline()
#                 if not header:
#                     break
#                 seq = r2.readline()
#                 plus = r2.readline()
#                 qual = r2.readline()
#                 header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
#                 o2.write(header)
#                 o2.write(convert_r2(seq))
#                 o2.write(plus)
#                 o2.write(qual)
#     return out, out2 if read2 else out

def bam_to_fastq_with_original_seq(bamfile, out1, out2=None):
    """
    Convert BAM to FASTQ, using the original sequence if present in the header (e.g., as ORIGINAL_SEQ:...), 
    otherwise use 'N' * read length. Output read 1 to out1, read 2 to out2 (if paired-end).
    """
    def open_out(filename):
        if filename and filename.endswith('.gz'):
            return gzip.open(filename, "wt")
        else:
            return open(filename, "w")

    fq1 = open_out(out1)
    fq2 = open_out(out2) if out2 else None

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            # Determine if read is read1 or read2
            if read.is_read1:
                fq = fq1
            elif read.is_read2 and fq2:
                fq = fq2
            else:
                fq = fq1
            parts = [read.query_name]
            # Try to extract original sequence from header
            orig_seq = None
            if "ORIGINAL_SEQ:" in read.query_name:
                # e.g. @name ORIGINAL_SEQ:ACGT...
                parts = read.query_name.split(" ORIGINAL_SEQ:")
                if len(parts) > 1:
                    orig_seq = parts[1].strip()
            # If not found, use N's of the same length as the read
            if not orig_seq or orig_seq == "":
                orig_seq = "N" * read.query_length
            fq.write(f"@{parts[0]}\n")
            fq.write(f"{orig_seq}\n")
            fq.write("+\n")
            qual = read.qual if read.qual else "I" * read.query_length
            if read.is_reverse:
                qual = qual[::-1]
            fq.write(f"{qual}\n")

    fq1.close()
    if fq2:
        fq2.close()

def process_chunk(records):
    out = []
    for header, seq, plus, qual in records:
        header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
        seq = seq.replace("C", "T").replace("c", "t")
        out.extend([header, seq, plus, qual])
    return out

def read_fastq_chunk(fh, chunk_size):
    records = []
    for _ in range(chunk_size):
        header = fh.readline()
        if not header:
            break
        seq = fh.readline()
        plus = fh.readline()
        qual = fh.readline()
        records.append((header, seq, plus, qual))
    return records

def process_chunk_r1(records):
    out = []
    for header, seq, plus, qual in records:
        header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
        seq = seq.replace("C", "T").replace("c", "t")
        out.extend([header, seq, plus, qual])
    return out

def process_chunk_r2(records):
    out = []
    for header, seq, plus, qual in records:
        header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
        seq = seq.replace("G", "A").replace("g", "a")
        out.extend([header, seq, plus, qual])
    return out

def parallel_convert_fastq(infile, outfile, chunk_size=100000, mode="r1"):
    opener = gzip.open if infile.endswith('.gz') else open
    process_chunk = process_chunk_r1 if mode == "r1" else process_chunk_r2
    with opener(infile, 'rt') as fh, gzip.open(outfile, 'wt') if outfile.endswith('.gz') else open(outfile, 'w') as out_fh:
        pool = Pool(cpu_count())
        while True:
            chunk = read_fastq_chunk(fh, chunk_size)
            if not chunk:
                break
            results = pool.map(process_chunk, [chunk])
            for recs in results:
                out_fh.writelines(recs)
        pool.close()
        pool.join()

def process_bam_chunk(records):
    out = []
    for read in records:
        if read.is_unmapped:
            continue
        parts = [read.query_name]
        # Try to extract original sequence from header
        orig_seq = None
        if "ORIGINAL_SEQ:" in read.query_name:
            parts = read.query_name.split(" ORIGINAL_SEQ:")
            if len(parts) > 1:
                orig_seq = parts[1].strip()
        if not orig_seq or orig_seq == "":
            sys.stderr.write(f"Warning: No original sequence found for read {read.query_name}.\n")
            sys.exit(1)
            #orig_seq = "N" * read.query_length
        fq_entry = f"@{parts[0]}\n{orig_seq}\n+\n"
        qual = read.qual if read.qual else "I" * read.query_length
        if read.is_reverse:
            qual = qual[::-1]
        fq_entry += f"{qual}\n"
        out.append((read.is_read1, fq_entry))
    return out

def read_bam_chunk(bam, chunk_size):
    records = []
    try:
        for _ in range(chunk_size):
            read = next(bam)
            records.append(read)
    except StopIteration:
        pass
    return records

# def parallel_bam_to_fastq(bamfile, out1, out2=None, chunk_size=100000):
#     def open_out(filename):
#         if filename and filename.endswith('.gz'):
#             return gzip.open(filename, "wt")
#         else:
#             return open(filename, "w")

#     fq1 = open_out(out1)
#     fq2 = open_out(out2) if out2 else None

#     with pysam.AlignmentFile(bamfile, "rb") as bam:
#         for read in bam.fetch(until_eof=True):
#             if read.is_unmapped:
#                 continue
#             # Determine if read is read1 or read2
#             if read.is_read1:
#                 fq = fq1
#             elif read.is_read2 and fq2:
#                 fq = fq2
#             else:
#                 fq = fq1
#             parts = [read.query_name]
#             # Try to extract original sequence from header
#             orig_seq = None
#             if "ORIGINAL_SEQ:" in read.query_name:
#                 parts = read.query_name.split(" ORIGINAL_SEQ:")
#                 if len(parts) > 1:
#                     orig_seq = parts[1].strip()
#             if not orig_seq or orig_seq == "":
#                 orig_seq = "N" * read.query_length
#             fq.write(f"@{parts[0]}\n")
#             fq.write(f"{orig_seq}\n")
#             fq.write("+\n")
#             qual = read.qual if read.qual else "I" * read.query_length
#             if read.is_reverse:
#                 qual = qual[::-1]
#             fq.write(f"{qual}\n")

#     fq1.close()
#     if fq2:
#         fq2.close()

def run_bbsplit(read1, read2, host, graft, out_host, out_graft, 
                bbsplit_index_build=1, bbsplit_index_path="bbsplit_index",
                bbsplit_path="bbsplit.sh", bbsplit_extra=""):
    """
    Run bbsplit.sh to split reads and require BAM output files.
    # ambiguous=best ambiguous2=split
    # AMBIGUOUS_bbsplit_graft.bam  AMBIGUOUS_bbsplit_host.bam
    """
    if not out_host.lower().endswith(".bam") or not out_graft.lower().endswith(".bam"):
        raise ValueError("Output files must have .bam suffix.")

    cmd = (
        f"{bbsplit_path} build={bbsplit_index_build} "
        f"path={bbsplit_index_path} "
        f"in={read1} "
        f"{f'in2={read2} ' if read2 else ''}"
        f"out_{host}={out_host} "
        f"out_{graft}={out_graft} "
        f"{bbsplit_extra}"
    )
    print(f"Running: {cmd}", file=sys.stderr)
    subprocess.check_call(shlex.split(cmd))

def build_bbsplit_index(host_fa, graft_fa, host_name, graft_name, bbsplit_idx_dir="bbsplit_idx_convert", bbsplit_path="bbsplit.sh", build="1"):
    """
    Build bbsplit index using the provided converted reference FASTA files and custom names.
    """
    cmd = (
        f"{bbsplit_path} path={bbsplit_idx_dir} build={build} "
        f"ref_{host_name}={host_fa} ref_{graft_name}={graft_fa}"
    )
    print(f"[bbsplit-build] CMD: {cmd}", file=sys.stdout)
    subprocess.check_call(cmd, shell=True)

def ensure_bbsplit_index_structure(bbsplit_index_path, build="1", host=None, graft=None):
    """
    1. Ensure the bbsplit index folder structure exists:
    <bbsplit_index_path>/ref/genome/<build>
    <bbsplit_index_path>/ref/index/<build>
    2. Also ensure that host and graft names exist in namelist.txt.
    """
    base = pathlib.Path(bbsplit_index_path)
    genome_dir = base / "ref" / "genome" / str(build)
    index_dir = base / "ref" / "index" / str(build)
    namelist_file = genome_dir / "namelist.txt"

    # Check directory structure
    if not genome_dir.is_dir() or not index_dir.is_dir():
        raise FileNotFoundError(
            f"Required bbsplit index directories do not exist:\n"
            f"  {genome_dir}\n"
            f"  {index_dir}\n"
            "Please build the bbsplit index first."
        )

    # Check host and graft names in namelist.txt if provided
    if host or graft:
        if not namelist_file.is_file():
            raise FileNotFoundError(f"namelist.txt not found: {namelist_file}")
        with open(namelist_file) as f:
            names = set(line.strip() for line in f if line.strip())
        missing = [n for n in (host, graft) if n and n not in names]
        if missing:
            raise ValueError(
                f"The following names are missing in {namelist_file}: {', '.join(missing)}"
            )

def fastq_count(filename):
    """
    Fast count of reads in a FASTQ file (gzipped or plain).
    """
    import subprocess
    import os

    if filename.endswith('.gz'):
        cmd = f"zcat {filename} | wc -l"
    else:
        cmd = f"wc -l < {filename}"
    n_lines = int(subprocess.check_output(cmd, shell=True).decode().strip())
    return n_lines // 4

def stat_split(raw_fastq, host_fastq, graft_fastq):
    """
    Output CSV with sample name, raw read number, host read number, graft read number, percent graft.
    Uses fastq_count for efficient counting.
    """
    import os

    sample_name = os.path.basename(raw_fastq).split('_')[0]
    n_raw = fastq_count(raw_fastq)
    n_host = fastq_count(host_fastq)
    n_graft = fastq_count(graft_fastq)
    percent_graft = (n_graft / n_raw * 100) if n_raw else 0

    print("Sample_name,raw_read_number,host_read_number,graft_read_number,percent_graft")
    print(f"{sample_name},{n_raw},{n_host},{n_graft},{percent_graft:.2f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sort reads into host and graft categories")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand 1: convert-ref
    parser_ref = subparsers.add_parser("convert-ref", help="Convert reference genome C->T (fwd) and G->A (rev)")
    parser_ref.add_argument("ref_fasta", help="Reference FASTA file")
    parser_ref.add_argument("-o", "--output", help="Output FASTA file (default: <ref_fasta>.c2t)")

    # Subcommand 2: convert-reads
    parser_reads = subparsers.add_parser("convert-reads", help="Convert C->T in read and G->A in read2 (if paired)")
    parser_reads.add_argument("--read", required=True, help="FASTQ file for read (single or read1)")
    parser_reads.add_argument("--read2", help="FASTQ file for read 2 (optional)")
    parser_reads.add_argument("--out", help="Output file for converted read")
    parser_reads.add_argument("--out2", help="Output file for converted read 2")

    # Subcommand 3: bam-to-fastq
    parser_bam2fq = subparsers.add_parser("bam-to-fastq", help="Convert BAM to FASTQ with sequence replaced by 'ORIGINAL_SEQ' or the original sequence if present in the header")
    parser_bam2fq.add_argument("bamfile", help="Input BAM file")
    parser_bam2fq.add_argument("--out", required=True, help="Output FASTQ file for read 1 (can be .gz)")
    parser_bam2fq.add_argument("--out2", help="Output FASTQ file for read 2 (can be .gz, optional)")

    # Subcommand 4: bbsplit
    parser_bbsplit = subparsers.add_parser("bbsplit", help="Invoke bbsplit.sh to split reads and convert to BAM")
    parser_bbsplit.add_argument("--read", required=True, help="FASTQ file for read 1")
    parser_bbsplit.add_argument("--read2", help="FASTQ file for read 2 (optional)")
    parser_bbsplit.add_argument("--host", required=True, help="Reference genome name for host")
    parser_bbsplit.add_argument("--graft", required=True, help="Reference genome name for graft")
    parser_bbsplit.add_argument("--out_host", required=True, help="Output BAM file for host")
    parser_bbsplit.add_argument("--out_graft", required=True, help="Output BAM file for graft")
    parser_bbsplit.add_argument("--bbsplit_path", default="bbsplit.sh", help="Path to bbsplit.sh")
    parser_bbsplit.add_argument("--bbsplit_extra", 
                                default="ambiguous=best ambiguous2=toss", 
                                help="Extra parameters for bbsplit.sh")
    parser_bbsplit.add_argument("--bbsplit_index_build", default=1, help="bbsplit index build (default: 1)")
    parser_bbsplit.add_argument("--bbsplit_index_path", default="bbsplit_index", help="bbsplit index path (default: bbsplit_index)")

    # Subcommand 5: bbsplit-index
    parser_bbsplit_build = subparsers.add_parser("bbsplit-index", help="Build bbsplit index from converted reference FASTA files")
    parser_bbsplit_build.add_argument("--host", required=True, help="Converted host reference FASTA file")
    parser_bbsplit_build.add_argument("--graft", required=True, help="Converted graft reference FASTA file")
    parser_bbsplit_build.add_argument("--host_name", required=True, help="Name for host reference (e.g. mm)")
    parser_bbsplit_build.add_argument("--graft_name", required=True, help="Name for graft reference (e.g. hs)")
    parser_bbsplit_build.add_argument("--bbsplit_path", default="bbsplit.sh", help="Path to bbsplit.sh")
    parser_bbsplit_build.add_argument("--bbsplit_index_path", default="bbsplit_idx_convert", help="bbsplit index directory")
    parser_bbsplit_build.add_argument("--bbsplit_index_build", default="1", help="Build number for bbsplit index (default: 1)")

    # Subcommand: stat-split
    parser_stat_split = subparsers.add_parser("stat-split", help="Output split statistics as CSV")
    parser_stat_split.add_argument("--raw", required=True, help="Raw FASTQ file (R1)")
    parser_stat_split.add_argument("--host", required=True, help="Host FASTQ file (R1)")
    parser_stat_split.add_argument("--graft", required=True, help="Graft FASTQ file (R1)")

    args = parser.parse_args()

    if args.command == "convert-ref":
        out_fa = args.output if args.output else args.ref_fasta + ".c2t"
        out_fa = convert_fasta(args.ref_fasta, out_fa=out_fa)
        print(out_fa)
    elif args.command == "convert-reads":
        # Old parallel_convert_fastq code (commented out)
        # out = args.out if args.out else args.read + ".meth"
        # parallel_convert_fastq(args.read, out, mode="r1")
        # if args.read2:
        #     out2 = args.out2 if args.out2 else args.read2 + ".meth"
        #     parallel_convert_fastq(args.read2, out2, mode="r2")

        # Use the fast awk/sed/gzip implementation
        out = args.out if args.out else args.read + ".meth"
        out2 = args.out2 if args.read2 and args.out2 else (args.read2 + ".meth" if args.read2 else None)
        convert_reads_c2t_r1_g2a_r2(args.read, args.read2, out, out2)
    elif args.command == "bam-to-fastq":
        bam_to_fastq_with_original_seq(args.bamfile, args.out, args.out2)
    elif args.command == "bbsplit":
        ensure_bbsplit_index_structure(
            args.bbsplit_index_path, args.bbsplit_index_build, host=args.host, graft=args.graft
        )
        run_bbsplit(
            args.read, args.read2, args.host, args.graft,
            args.out_host, args.out_graft,
            bbsplit_path=args.bbsplit_path,
            bbsplit_extra=args.bbsplit_extra,
            bbsplit_index_build=args.bbsplit_index_build,
            bbsplit_index_path=args.bbsplit_index_path
        )
    elif args.command == "bbsplit-index":
        build_bbsplit_index(
            host_fa=args.host,
            graft_fa=args.graft,
            host_name=args.host_name,
            graft_name=args.graft_name,
            bbsplit_idx_dir=args.bbsplit_index_path,
            bbsplit_path=args.bbsplit_path,
            build=args.bbsplit_index_build
        )
    elif args.command == "stat-split":
        stat_split(args.raw, args.host, args.graft)

