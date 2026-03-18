#!/bin/bash
#SBATCH --job-name=quality_control
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

# Assess Quality using FASTQC
module load fastqc

mkdir -p data/fastqc

for i in data/raw/*.fastq
do
  echo "Running FASTQC on ${i}" 

  fastqc ${i} -o data/fastqc -t $SLURM_CPUS_PER_TASK

done