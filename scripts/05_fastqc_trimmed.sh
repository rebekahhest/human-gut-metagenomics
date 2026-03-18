#!/bin/bash
#SBATCH --job-name=quality_control_post-trimming
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

# Assess Quality using FASTQC
module load fastqc

mkdir -p data/fastqc_trimmed

for i in data/trimmed/*_paired.fastq
do
  echo "Running FASTQC on ${i}" 

  fastqc ${i} -o data/fastqc_trimmed -t $SLURM_CPUS_PER_TASK

done