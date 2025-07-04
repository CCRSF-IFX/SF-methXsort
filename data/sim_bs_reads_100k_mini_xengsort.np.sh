Building DAG of jobs...
Job counts:
	count	jobs
	1	all
	1	convert_ref
	1	convert_restore_fastq
	1	merge_graft_host
	1	sherman
	1	stat_accuracy
	1	stat_read_number
	1	xengsort_classify
	1	xengsort_idx
	9

[Fri Jul  4 00:03:17 2025]
rule convert_ref:
    input: /mnt/ccrsf-ifx/RefGenomes/SingleCell_REF/10X/hg38/refdata-gex-GRCh38-2024-A/fasta/genome.fa, /mnt/ccrsf-ifx/RefGenomes/mm10_all/mm10.fdb.fa
    output: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/graft_cvt.fa, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/host_cvt.fa
    jobid: 4


python /mnt/ccrsf-ifx/Software/github/methXsort/methXsort.py convert-ref     /mnt/ccrsf-ifx/RefGenomes/SingleCell_REF/10X/hg38/refdata-gex-GRCh38-2024-A/fasta/genome.fa --out /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/graft_cvt.fa 
python /mnt/ccrsf-ifx/Software/github/methXsort/methXsort.py convert-ref     /mnt/ccrsf-ifx/RefGenomes/mm10_all/mm10.fdb.fa --out /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/host_cvt.fa


[Fri Jul  4 00:03:17 2025]
rule sherman:
    input: /mnt/ccrsf-ifx/RefGenomes/SingleCell_REF/10X/hg38/refdata-gex-GRCh38-2024-A/fasta/genome.fa, /mnt/ccrsf-ifx/RefGenomes/mm10_all/mm10.fdb.fa
    output: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/graft/sherman.log, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/host/sherman.log, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_graft_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_graft_R2.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_host_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_host_R2.fastq.gz
    jobid: 7


mkdir -p /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/graft && ln -s /mnt/ccrsf-ifx/RefGenomes/SingleCell_REF/10X/hg38/refdata-gex-GRCh38-2024-A/fasta/genome.fa /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/graft/graft.fa
mkdir -p /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/host && ln -s /mnt/ccrsf-ifx/RefGenomes/mm10_all/mm10.fdb.fa /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/host/host.fa
cd  /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/graft &&  /mnt/ccrsf-ifx/Software/tools/Sherman/v0.1.8/Sherman --conversion_rate 60 -l 150 -n 100000 --genome_folder /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/graft -pe > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/graft/sherman.log 2>&1 & 
cd  /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/host &&  /mnt/ccrsf-ifx/Software/tools/Sherman/v0.1.8/Sherman --conversion_rate 60 -l 150 -n 100000 --genome_folder /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/host -pe > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/host/sherman.log 2>&1 &
wait
cd ../ && 
cat /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/graft/simulated_1.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\1graft_/' |gzip -c > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_graft_R1.fastq.gz
cat /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/graft/simulated_2.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\1graft_/' |gzip -c > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_graft_R2.fastq.gz
cat /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/host/simulated_1.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\1host_/' |gzip -c > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_host_R1.fastq.gz
cat /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/host/simulated_2.fastq |sed '1~4s/^@\([0-9][0-9]*_\)/@\1host_/' |gzip -c > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_host_R2.fastq.gz


[Fri Jul  4 00:03:17 2025]
rule merge_graft_host:
    input: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_graft_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_graft_R2.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_host_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_host_R2.fastq.gz
    output: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bs_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bs_R2.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R2.fastq.gz
    jobid: 5


cat /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_graft_R1.fastq.gz /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_host_R1.fastq.gz > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bs_R1.fastq.gz
cat /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_graft_R2.fastq.gz /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/simulated_bs_host_R2.fastq.gz > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bs_R2.fastq.gz
python /mnt/ccrsf-ifx/Software/github/methXsort/methXsort.py convert-reads --with_orig_seq     --read /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bs_R1.fastq.gz --read2 /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bs_R2.fastq.gz     --out /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R1.fastq.gz --out2 /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R2.fastq.gz


[Fri Jul  4 00:03:17 2025]
rule xengsort_idx:
    input: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/graft_cvt.fa, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/host_cvt.fa
    output: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/xengsort_index.log
    jobid: 1


python /mnt/ccrsf-ifx/Software/github/methXsort/methXsort.py xengsort-index     --host /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/host_cvt.fa --graft /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/graft_cvt.fa     --index /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/xengsort_index     -n 7_000_000_000 --fill 0.88 --statistics summary -k 25     > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/xengsort_index.log 2>&1


[Fri Jul  4 00:03:17 2025]
rule xengsort_classify:
    input: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R2.fastq.gz
    output: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/classify.log, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-graft.1.fq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-graft.2.fq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-host.1.fq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-host.2.fq.gz
    jobid: 8


mkdir -p /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify && cd /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify && python /mnt/ccrsf-ifx/Software/github/methXsort/methXsort.py xengsort-classify     --read /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R1.fastq.gz --read2 /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R2.fastq.gz     --index /mnt/ccrsf-ifx/RefGenomes/dragen_ref/xenoxsort/methxsort_xengsort_index_converted/methxsort_index_xs_convert     --out_prefix /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify     --threads 33     > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/classify.log 2>&1


[Fri Jul  4 00:03:17 2025]
rule convert_restore_fastq:
    input: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-graft.1.fq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-graft.2.fq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-host.1.fq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-host.2.fq.gz
    output: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_graft_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_graft_R2.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_host_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_host_R2.fastq.gz
    jobid: 6


python /mnt/ccrsf-ifx/Software/github/methXsort/methXsort.py restore-fastq     --read /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-graft.1.fq.gz --read2 /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-graft.2.fq.gz     --out /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_graft_R1.fastq.gz --out2 /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_graft_R2.fastq.gz 
python /mnt/ccrsf-ifx/Software/github/methXsort/methXsort.py restore-fastq     --read /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-host.1.fq.gz --read2 /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/xengsort_classify/xengsort_classify-host.2.fq.gz     --out /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_host_R1.fastq.gz --out2 /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_host_R2.fastq.gz


[Fri Jul  4 00:03:17 2025]
rule stat_read_number:
    input: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_graft_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_host_R1.fastq.gz
    output: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/read_number_stat.txt
    jobid: 2


python /mnt/ccrsf-ifx/Software/github/methXsort/methXsort.py stat-split --raw /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R1.fastq.gz              --graft /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_graft_R1.fastq.gz --host /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_host_R1.fastq.gz               > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/read_number_stat.txt


[Fri Jul  4 00:03:17 2025]
rule stat_accuracy:
    input: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/sim_reads/merged_simulated_bscvt_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_graft_R1.fastq.gz, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_host_R1.fastq.gz
    output: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/accuracy_stat.txt
    jobid: 3


python /mnt/ccrsf-ifx/Software/github/methXsort/scripts/simulate_reads/stat_accuracy.py 100000              /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_host_R1.fastq.gz,host  /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/fastq/bbsplit_graft_R1.fastq.gz,graft               > /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/accuracy_stat.txt


[Fri Jul  4 00:03:17 2025]
localrule all:
    input: /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/ref_idx/xengsort_index.log, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/read_number_stat.txt, /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k_xengsort/accuracy_stat.txt
    jobid: 0

Job counts:
	count	jobs
	1	all
	1	convert_ref
	1	convert_restore_fastq
	1	merge_graft_host
	1	sherman
	1	stat_accuracy
	1	stat_read_number
	1	xengsort_classify
	1	xengsort_idx
	9
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
