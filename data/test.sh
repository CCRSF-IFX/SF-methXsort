



 paste <(zcat sim_reads/merged_simulated_bs_R1.fastq.gz) <(zcat sim_reads/merged_simulated_bs_R2.fastq.gz)  |awk 'NR%4==1{h1=$1; h2=$2} NR%4==2{s1=$1; s2=$2} NR%4==3{p1=$1} NR%4==0{q1=$1; print h1 " ORIGINAL_SEQ_R1:" s1 " ORIGINAL_SEQ_R2:" s2; print s1; print p1; print q1}' |head -1000 | sed '2~4s/G/A/g;2~4s/g/a/g' | gzip > test_R1.fq.gz

zcat sim_reads/merged_simulated_bs_R2.fastq.gz |sed '2~4s/G/A/g;2~4s/g/a/g' |head -1000 |gzip > test_R2.fq.gz

python /mnt/ccrsf-ifx/Software/github/methXsort/methXsort.py bbsplit    \
	 --read test_R1.fq.gz \ 
         --read2  test_R2.fq.gz \
	  --host ecoli --graft ehec     --bbsplit_index_build 1  \ 
         --bbsplit_index_path /mnt/ccrsf-ifx/Software/github/methXsort/data/sherman_reads_100k/ref_idx/bbsplit_index    \
	 --out_host test_host.bam     \
	 --out_graft test_graft.bam    \
         --bbsplit_extra "scafstats=scafstats.log refstats=refstats.log"  > bbsplit.log 2>&1
