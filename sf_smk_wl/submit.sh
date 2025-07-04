#!/bin/bash
#SBATCH --partition=norm
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=4g
#SBATCH --mail-user=xies4@nih.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=96:00:00
#SBATCH --no-requeue
#SBATCH --job-name=CS039503

#2019/12/27 shent2: --no-requeue is added
export PATH=/mnt/nasapps/development/python/3.7.1/bin:$PATH
snakemake --jobname '39503.{jobid}.{rulename}' --latency-wait 600 -j 100 --rerun-incomplete --keep-going --stats snakemake.stats --printshellcmds --cluster-config cluster.json --cluster "sbatch --partition={cluster.partition} --nodes=1 --ntasks={cluster.threads} --mem={cluster.mem}  --time={cluster.time} --no-requeue"   >& snakemake.log
