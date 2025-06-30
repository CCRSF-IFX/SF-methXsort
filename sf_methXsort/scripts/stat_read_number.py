#!/usr/bin/env python3

import sys
import gzip
import os

def open_fastq(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def count_reads(fastq_file):
    """
    Count the number of reads in a FASTQ file.
    """
    count = 0
    with open_fastq(fastq_file) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            count += 1
    return count

def main():
    if len(sys.argv) != 4:
        print("Usage: python stat_read_number.py <raw_R1.fastq.gz> <host_R1.fastq.gz> <graft_R1.fastq.gz>")
        sys.exit(1)

    raw_fastq = sys.argv[1]
    host_fastq = sys.argv[2]
    graft_fastq = sys.argv[3]

    # Use the base name of the raw fastq file as the sample name
    sample_name = os.path.basename(raw_fastq).split('_')[0]

    n_raw = count_reads(raw_fastq)
    n_host = count_reads(host_fastq)
    n_graft = count_reads(graft_fastq)
    percent_graft = (n_graft / n_raw * 100) if n_raw else 0

    # Output CSV header and values
    print("Sample_name,raw_read_number,host_read_number,graft_read_number,percent_graft")
    print(f"{sample_name},{n_raw},{n_host},{n_graft},{percent_graft:.2f}")

if __name__ == "__main__":
    main()