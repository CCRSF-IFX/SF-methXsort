#!/bin/bash
#SBATCH --partition=norm
#SBATCH --job-name=test_methxsort 
#SBATCH --nodes=1
#SBATCH --ntasks=34
#SBATCH --mail-user=xies4@nih.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=96:00:00

module load miniconda
conda activate /mnt/ccrsf-ifx/Software/tools/methXsort/v0.1.0
 snakemake -s ../scripts/simulate_reads/sim_bs_reads_sherman_xengsort.smk --configfile ../scripts/simulate_reads/sim_bs_reads_100k_mini_xengsort.yaml -np > sim_bs_reads_100k_mini_xengsort.np.sh

  snakemake -s ../scripts/simulate_reads/sim_bs_reads_sherman_xengsort.smk --configfile ../scripts/simulate_reads/sim_bs_reads_100k_mini_xengsort.yaml -j 20 > sim_bs_reads_100k_mini_xengsort.log 2>&1 
