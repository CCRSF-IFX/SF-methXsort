#!/usr/bin/env python3

import toolshed
from toolshed import nopen, reader, is_newer_b
import argparse
import sys
import os
from itertools import groupby
import gzip
import pysam

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
    Convert C->T in read (single or read1) and G->A in read2 (if paired-end).
    Store the original sequence in the header line as an extra field.
    """
    def convert_r1(seq):
        return seq.replace("C", "T").replace("c", "t")
    def convert_r2(seq):
        return seq.replace("G", "A").replace("g", "a")
    if out is None:
        out = read + ".meth"
    if read2 and out2 is None:
        out2 = read2 + ".meth"

    def open_out(filename):
        if filename and filename.endswith('.gz'):
            return gzip.open(filename, "wt")
        else:
            return open(filename, "w")

    with nopen(read) as r1, open_out(out) as o1:
        while True:
            header = r1.readline()
            if not header:
                break
            seq = r1.readline()
            plus = r1.readline()
            qual = r1.readline()
            # Add original sequence to header
            header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
            o1.write(header)
            o1.write(convert_r1(seq))
            o1.write(plus)
            o1.write(qual)
    if read2:
        with nopen(read2) as r2, open_out(out2) as o2:
            while True:
                header = r2.readline()
                if not header:
                    break
                seq = r2.readline()
                plus = r2.readline()
                qual = r2.readline()
                header = header.rstrip() + f" ORIGINAL_SEQ:{seq.strip()}\n"
                o2.write(header)
                o2.write(convert_r2(seq))
                o2.write(plus)
                o2.write(qual)
    return out, out2 if read2 else out

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

    args = parser.parse_args()

    if args.command == "convert-ref":
        out_fa = args.output if args.output else args.ref_fasta + ".c2t"
        out_fa = convert_fasta(args.ref_fasta, out_fa=out_fa)
        print(out_fa)
    elif args.command == "convert-reads":
        convert_reads_c2t_r1_g2a_r2(args.read, args.read2, args.out, args.out2)
    elif args.command == "bam-to-fastq":
        bam_to_fastq_with_original_seq(args.bamfile, args.out, args.out2)
