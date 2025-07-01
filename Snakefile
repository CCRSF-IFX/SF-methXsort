#Snakefile for sorting methylation reads into
# human reads (graft) and mouse reads (host)
shell.executable("/bin/bash")
shell.prefix("source /etc/profile.d/modules.sh; ")

import config
import program
import reference
from snakemake.utils import R
import glob
import os, re, sys
from os import listdir
from os.path import isfile, isdir, join
import os.path as path
from xml.dom import minidom

copydir = "/mnt/ccrsf-ifx/Report_archive/report_archive_illumina"
makefilebase = config.analysis + "/fastq/"
project_name = os.path.basename(config.analysis)
searchterms = project_name.split("_")
date = searchterms[-1]
pi_name = searchterms[0]
csas = '00000'
for searchterm in searchterms:
    if searchterm.startswith("CS0"):
        csas = searchterm[-5:]

analysis = config.analysis
unaligned = config.unaligned
unalignedfolder = "/" + "/".join(unaligned.split("/")[1:len(unaligned.split("/")) - 1])
one_up = path.abspath(path.join(os.getcwd(),"../"))
include: program.runParametersImport
run_name = analysis.split('/')[-2]
flowcell = run_name[-9:]
try:
    (RTAVersion, flowcell, workFlowType, flowcellMode, chemistry, chemistryVersion, xmlRunParametersPath, xmlRunInfoPath) = runParametersXmlPath(run_name)
except IOError as error:
    xmlRunParametersPath = "Unknown"
    xmlRunInfoPath = "Unknown"
    sys.stdout.write("\n\n\nExecution failed: " + str(error) +"\n")
    sys.stdout.write("No RunParameters.xml and RunInfo.xml will be archied and flowcell ID is geneated from the analysis folder\n\n\n")
except:
    sys.stdout.write("Unexpected error:" + str(sys.exc_info()[0]) +"\n")
report_result = one_up + "/" + project_name + "_" + flowcell + ".xlsx"
wreport_result = one_up + "/" + project_name + "_" + flowcell + ".docx"
copy_result = one_up + "/" + project_name + "_" + flowcell + "_copy.txt"

if hasattr(config, 'forcestrand'):
	forcestrand = config.forcestrand
else:
	forcestrand = ""

#archtar = "/is2/projects/CCR-SF/archive/illumina/" + run_name + "_" + project_name + ".tar.gz"
metalog = one_up + "/csafe_" + project_name + ".log"
scperr = one_up + "/" + project_name + "_scp.err"
log2 = one_up + "/csafe_" + project_name + "_2.log"

sample = [os.path.basename(file).split('.')[0] for file in glob.glob(makefilebase+'/*')]
samps = []
i=1
for item in sample:
        newvar = item.split("_R1")
        othervar = item.split("_R2")
        samps.append(newvar[0])
new = []
for item in samps:
        if '_R2_' not in item:
                new.append(item)
samples = [s.replace('Sample_', '') for s in new]
print(samples)
with open(config.analysis + '/cluster.json') as file:
    clusterConfig = json.load(file)

rule all:
    input: 
        expand("Sample_{sample}/{sample}_read_number.txt", sample=samples), 
        one_up + "/meta2json_complete.txt",
        
        
rule convert_reads:
    input: 
        R1 = makefilebase + "{sample}_R1_001.fastq.gz", 
        R2 = makefilebase + "{sample}_R2_001.fastq.gz"
    output: 
        R1 = "Sample_{sample}/{sample}_converted_R1.fastq.gz", 
        R2 = "Sample_{sample}/{sample}_converted_R2.fastq.gz", 
    params: batch = "-l nodes=1:ppn=16,mem=64g", prefix = "Sample_{sample}/logs/" 
    shell: 
        """
module load miniconda
conda activate {program.methXsort_env_path}
python methXsort.py convert-reads --read {input.R1} --read2 {input.R2} \
              --out {output.R1} --out2 {output.R2} 
"""

rule bbsplit:
    input: 
        reads_R1 = "Sample_{sample}/{sample}_converted_R1.fastq.gz", 
        reads_R2 = "Sample_{sample}/{sample}_converted_R2.fastq.gz"
    output: 
        bbsplit_log = "Sample_{sample}/logs/{sample}_bbsplit.log",
        bam_hs = "Sample_{sample}/{sample}_hs.bam",
        bam_mm = "Sample_{sample}/{sample}_mm.bam",
        scafstats = "Sample_{sample}/{sample}_scafstats.txt",
        refstats = "Sample_{sample}/{sample}_refstats.txt"
    params: 
        batch = "-l nodes=1:ppn=16,mem=64g", 
        bbsplit_idx = reference.bbsplit_idx,
        bbsplit_path = reference.bbsplit_path,
    shell: 
        """
        module load miniconda
        conda activate {program.methXsort_env_path}
        python methXsort.py bbsplit \
            --read {input.reads_R1} \
            --read2 {input.reads_R2} \
            --host mm \
            --graft hs \
            --out_host {output.bam_mm} \
            --out_graft {output.bam_hs} \
            --bbsplit_path {params.bbsplit_path} \
            --bbsplit_index_build 1 \
            --bbsplit_index_path {params.bbsplit_idx} \
            --bbsplit_extra "scafstats={output.scafstats} refstats={output.refstats}" \
            > {output.bbsplit_log} 2>&1
        """

rule convert_bam_to_fastq:
    input:
        bam_hs = "Sample_{sample}/{sample}_hs.bam",
        bam_mm = "Sample_{sample}/{sample}_mm.bam"
    output:
        fastq_hs_R1 = "Sample_{sample}/{sample}_hs_R1_001.fastq.gz",
        fastq_hs_R2 = "Sample_{sample}/{sample}_hs_R2_001.fastq.gz",
        fastq_mm_R1 = "Sample_{sample}/{sample}_mm_R1_001.fastq.gz",
        fastq_mm_R2 = "Sample_{sample}/{sample}_mm_R2_001.fastq.gz"
    params: batch = "-l nodes=1:ppn=16,mem=64g",
    shell:
        """
        module load miniconda
        conda activate {program.methXsort_env_path}
        python methXsort.py bam-to-fastq --out {output.fastq_hs_R1} --out2 {output.fastq_hs_R2} {input.bam_hs}
        python methXsort.py bam-to-fastq --out {output.fastq_mm_R1} --out2 {output.fastq_mm_R2} {input.bam_mm}
        """

rule stat_read_number: 
    input: 
        raw_R1 = makefilebase + "{sample}_R1_001.fastq.gz", 
        fastq_hs_R1 = rules.convert_bam_to_fastq.output.fastq_hs_R1,
        fastq_mm_R1 = rules.convert_bam_to_fastq.output.fastq_mm_R1
    output: 
        "Sample_{sample}/{sample}_read_number.txt"
    params: batch = "-l nodes=1:ppn=16,mem=64g"
    shell: 
        """
        module load miniconda
        conda activate {program.methXsort_env_path}
        python methXsort.py stat-split {input.raw_R1} {input.fastq_hs_R1} {input.fastq_mm_R1} \
           > {output}
        """
	
rule meta2json:
    input: 
        expand(makefilebase + "{sample}_R1_001.fastq.gz", sample=samples), 
        expand(rules.stat_read_number.output, sample=samples)
    output: one_up + "/meta2json_complete.txt"
    params: batch = "-l nodes=1:ppn=4,mem=8g", prefix = project_name + "_" + flowcell + "_Metadata.txt"
    run:
        import os, re, sys
        command = f'cd {one_up}; python /mnt/ccrsf-ifx/Software/scripts/bin/run_meta2json.py -m {params.prefix} -uf {unalignedfolder} -r {run_name} -af {project_name}'
        shell(f'{command}')
        command = f'cd {one_up}; echo meta2json is complete > meta2json_complete.txt'
        with open(f'{one_up}/csafe_{project_name}_fastqFileList.log', 'r') as LOG:
            for line in LOG:
                if 'Status: FAILED' in line:
                    command = f'cd {one_up}; false'
        shell(f'{command}')
