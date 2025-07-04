# SF-methXsort

**methXsort** is a command-line toolkit for sorting, converting, and analyzing bisulfite sequencing reads in xenograft experiments. It provides fast, parallelized methods for reference and read conversion, species-specific read assignment, BAM-to-FASTQ extraction, and summary/statistics utilities for workflow integration.

---

## Installation

Clone this repository and ensure all dependencies (Python 3, `pysam`, `toolshed`, `gzip`, and required external tools like `bbsplit.sh`, `filterbyname.sh`, and `xengsort`) are available in your environment.

```bash
git clone <repo_url>
cd methXsort
```

---

## Usage

All commands are run via the main script:

```bash
python methXsort.py <subcommand> [options]
```

### Subcommands

#### 1. Convert Reference Genome

Convert a reference FASTA for bisulfite mapping (C→T and G→A):

```bash
python methXsort.py convert-ref <ref_fasta> [-o OUTPUT]
```

#### 2. Convert Reads

Convert reads for bisulfite mapping (C→T for R1, G→A for R2):

```bash
python methXsort.py convert-reads --read <R1.fastq.gz> [--read2 <R2.fastq.gz>] [--out <R1_out>] [--out2 <R2_out>] [--with_orig_seq]
```
- `--with_orig_seq`: Store the original sequence in the header (slower, but traceable).

#### 3. Filter FASTQ by BAM

Extract reads from FASTQ that are present in a BAM file (e.g., after bbsplit):

```bash
python methXsort.py filter-fastq-by-bam --read <R1.fastq.gz> [--read2 <R2.fastq.gz>] --bam <file.bam> --out <R1_out> [--out2 <R2_out>] [--filterbyname_path <path>]
```

#### 4. Run bbsplit

Split reads into host and graft using bbsplit:

```bash
python methXsort.py bbsplit --read <R1.fastq.gz> [--read2 <R2.fastq.gz>] --host <host_name> --graft <graft_name> --out_host <host.bam> --out_graft <graft.bam> [--bbsplit_path <path>] [--bbsplit_extra <extra>]
```

#### 5. Build bbsplit Index

```bash
python methXsort.py bbsplit-index --host <host.fa> --graft <graft.fa> --host_name <host> --graft_name <graft> [--bbsplit_path <path>] [--bbsplit_index_path <dir>]
```

#### 6. Split Statistics

Output CSV statistics for split reads:

```bash
python methXsort.py stat-split --raw <raw_R1.fastq.gz> --host <host_R1.fastq.gz> --graft <graft_R1.fastq.gz>
```

#### 7. Build xengsort Index

```bash
python methXsort.py xengsort-index --host <host.fa> --graft <graft.fa> --index <index_dir> [-n N] [--fill FILL] [--statistics STAT] [-k K] [--xengsort_path <path>] [--xengsort_extra <extra>]
```

#### 8. Classify Reads with xengsort

```bash
python methXsort.py xengsort-classify --read <R1.fastq.gz> [--read2 <R2.fastq.gz>] --index <index_dir> --out_prefix <prefix> --threads <N> [--xengsort_path <path>] [--xengsort_extra <extra>]
```

#### 9. Restore FASTQ from xengsort Output

Restore original sequences in FASTQ files classified by xengsort:

```bash
python methXsort.py restore-fastq --read <classified_R1.fq.gz> --out <restored_R1.fq.gz> [--read2 <classified_R2.fq.gz> --out2 <restored_R2.fq.gz>]
```

---

## Example Workflow

1. **Convert reference genomes:**
    ```bash
    python methXsort.py convert-ref mm10.fa -o mm10_converted.fa
    python methXsort.py convert-ref hg38.fa -o hg38_converted.fa
    ```

2. **Build bbsplit and xengsort indices:**
    ```bash
    python methXsort.py bbsplit-index --host mm10_converted.fa --graft hg38_converted.fa --host_name mm --graft_name hs
    python methXsort.py xengsort-index --host mm10_converted.fa --graft hg38_converted.fa --index xengsort_index_7B
    ```

3. **Convert reads:**
    ```bash
    python methXsort.py convert-reads --read sample_R1.fastq.gz --read2 sample_R2.fastq.gz --with_orig_seq
    ```

4. **Run bbsplit or xengsort:**
    ```bash
    python methXsort.py bbsplit --read sample_R1.meth.gz --read2 sample_R2.meth.gz --host mm --graft hs --out_host host.bam --out_graft graft.bam
    python methXsort.py xengsort-classify --read sample_R1.meth.gz --read2 sample_R2.meth.gz --index xengsort_index_7B --out_prefix sample_xengsort --threads 8
    ```

5. **Restore original FASTQ:**
    ```bash
    python methXsort.py restore-fastq --read sample_xengsort-graft.1.fq.gz --out sample_graft_R1_restored.fq.gz --read2 sample_xengsort-graft.2.fq.gz --out2 sample_graft_R2_restored.fq.gz
    ```

---

## Notes

- For all subcommands, use `-h` or `--help` to see detailed options.
- Make sure all required external tools (`bbsplit.sh`, `filterbyname.sh`, `xengsort`) are in your `PATH` or specify their locations with the appropriate options.
- For paired-end data, always provide both `--read2` and `--out2` where required.

---

## Contact

Email: ccrsfifx@nih.gov

