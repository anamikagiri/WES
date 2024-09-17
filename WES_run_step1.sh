#!/bin/bash
# To run this script type sh run_step1.sh in the command line. DO NOT USE SBATCH OR SRUN HERE!

# Make sure you use the right path to the folder containing the fastq.gz files

for f in Fastq_files/*_R1.fastq.gz
do
	sbatch WES_step1_generate_gvcf_files.sh $f
done