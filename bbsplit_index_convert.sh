#!/bin/bash
#SBATCH --partition=norm
#SBATCH --job-name=wgsim 
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-user=xies4@nih.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=96:00:00

module load miniconda
conda activate /mnt/ccrsf-ifx/Software/tools/methXsort/v0.1.0/
python /mnt/ccrsf-ifx/Software/github/methXsort/methxsort.py convert-ref -o mm10_converted.fa /mnt/ccrsf-ifx/RefGenomes/mm10_all/mm10.fdb.fa > mm10_converted.fa.log 2>&1 & 
python /mnt/ccrsf-ifx/Software/github/methXsort/methxsort.py convert-ref -o hg38_converted.fa /mnt/ccrsf-ifx/RefGenomes/SingleCell_REF/10X/hg38/refdata-gex-GRCh38-2024-A/fasta/genome.fa > hg38_converted.fa.log 2>&1 &
wait

bbsplit.sh path=bbsplit_idx_convert/ build=1  \
        ref_mm=mm10_converted.fa \
        ref_hs=hg38_converted.fa > bbsplit.log 2>&1 
