hg38 = config["hg38"]
mm10 = config["mm10"]
outdir = config["outdir"]
read_number = config["read_pair_number"]
bbsplit_path = config["bbsplit_idx_path"]
bbsplit_idx = config["bbsplit_idx"]

sherman_path = "/mnt/ccrsf-ifx/Software/tools/Sherman/v0.1.8/Sherman"

rule all:
    input:
        expand(os.path.join(outdir, "fastq/bbsplit_{genome}_{pair}.fastq.gz"), 
               genome=["hg38", "mm10"], 
               pair =["R1", "R2"]),

rule sherman: 
    input: 
        hg38 = hg38,
        mm10 = mm10,
    params:
        od_hg38 = os.path.join(outdir, "hg38/"),
        od_mm10 = os.path.join(outdir, "mm10/"),
    output:
        hg38 = os.path.join(outdir, "hg38/sherman.log"), 
        mm10 = os.path.join(outdir, "mm10/sherman.log"),
        sim_hg38_R1 = os.path.join(outdir, "sim_reads/simulated_bs_hg38_R1.fastq.gz"),
        sim_hg38_R2 = os.path.join(outdir, "sim_reads/simulated_bs_hg38_R2.fastq.gz"),
        sim_mm10_R1 = os.path.join(outdir, "sim_reads/simulated_bs_mm10_R1.fastq.gz"),
        sim_mm10_R2 = os.path.join(outdir, "sim_reads/simulated_bs_mm10_R2.fastq.gz"),
    shell:
        """
ln -s {input.hg38} {params.od_hg38}/hg38.fa
ln -s {input.mm10} {params.od_mm10}/mm10.fa
mkdir -p {params.od_hg38} && cd  {params.od_hg38} && \
 {sherman_path} -l 150 -n {read_number} --genome_folder {params.od_hg38} -pe > {output.hg38} 2>&1 & 
mkdir -p {params.od_mm10} && cd  {params.od_mm10} && \
 {sherman_path} -l 150 -n {read_number} --genome_folder {params.od_mm10} -pe > {output.mm10} 2>&1 &
wait
cd ../ && 
cat {params.od_hg38}/simulated_1.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\\1hg38_/' |gzip -c > {output.sim_hg38_R1}
cat {params.od_hg38}/simulated_2.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\\1hg38_/' |gzip -c > {output.sim_hg38_R2}
cat {params.od_mm10}/simulated_1.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\\1mm10_/' |gzip -c > {output.sim_mm10_R1}
cat {params.od_mm10}/simulated_2.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\\1mm10_/' |gzip -c > {output.sim_mm10_R2}
"""

rule merge_hg38_mm10:
    input:
        hg38_R1 = rules.sherman.output.sim_hg38_R1,
        hg38_R2 = rules.sherman.output.sim_hg38_R2,
        mm10_R1 = rules.sherman.output.sim_mm10_R1,
        mm10_R2 = rules.sherman.output.sim_mm10_R2,
    output:
        merged_R1 = os.path.join(outdir, "sim_reads/merged_simulated_bs_R1.fastq.gz"),
        merged_R2 = os.path.join(outdir, "sim_reads/merged_simulated_bs_R2.fastq.gz"),
        merged_cvt_R1 = os.path.join(outdir, "sim_reads/merged_simulated_bscvt_R1.fastq.gz"),
        merged_cvt_R2 = os.path.join(outdir, "sim_reads/merged_simulated_bscvt_R2.fastq.gz"),
    shell:
        """
cat {input.hg38_R1} {input.mm10_R1} > {output.merged_R1}
cat {input.hg38_R2} {input.mm10_R2} > {output.merged_R2}
python /mnt/ccrsf-ifx/Software/github/methXsort/methxsort.py convert-reads \
    --read {output.merged_R1} --read2 {output.merged_R2} \
    --out {output.merged_cvt_R1} --out2 {output.merged_cvt_R2}
"""

rule bbsplit:
    input:
        reads_R1 = rules.merge_hg38_mm10.output.merged_cvt_R1,
        reads_R2 = rules.merge_hg38_mm10.output.merged_cvt_R2,
    output:
        bbsplit_log = os.path.join(outdir, "bbsplit.log"),
        bam_hg38 = os.path.join(outdir, "bbsplit/bbsplit_hg38.bam"),
        bam_mm10 = os.path.join(outdir, "bbsplit/bbsplit_mm10.bam"),
        scafstats = os.path.join(outdir, "bbsplit/scafstats.log"),
        refstats = os.path.join(outdir, "bbsplit/refstats.log"),
    shell:
        """
bbsplit.sh build={bbsplit_idx} path={bbsplit_path} \
    in={input.reads_R1} in2={input.reads_R2} \
    out_mm10={output.bam_mm10} out_hg38={output.bam_hg38}  \
    minhits=1 ambiguous2=all \
    scafstats={output.scafstats} refstats={output.refstats} \
    > {output.bbsplit_log} 2>&1
"""

rule convert_bam_to_fastq:
    input:
        bam_hg38 = rules.bbsplit.output.bam_hg38,
        bam_mm10 = rules.bbsplit.output.bam_mm10,
    output:
        fastq_hg38_R1 = os.path.join(outdir, "fastq/bbsplit_hg38_R1.fastq.gz"),
        fastq_hg38_R2 = os.path.join(outdir, "fastq/bbsplit_hg38_R2.fastq.gz"),
        fastq_mm10_R1 = os.path.join(outdir, "fastq/bbsplit_mm10_R1.fastq.gz"),
        fastq_mm10_R2 = os.path.join(outdir, "fastq/bbsplit_mm10_R2.fastq.gz"),
    shell:
        """
python /mnt/ccrsf-ifx/Software/github/methXsort/methxsort.py bam-to-fastq --out {output.fastq_hg38_R1} --out2 {output.fastq_hg38_R2} {input.bam_hg38}
python /mnt/ccrsf-ifx/Software/github/methXsort/methxsort.py bam-to-fastq --out {output.fastq_mm10_R1} --out2 {output.fastq_mm10_R2} {input.bam_mm10}
"""

