#!/bin/bash
#SBATCH --job-name=trim_adapters
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

# Trim adapters for QC

cd data
mkdir -p trimmed

module load trimmomatic

for i in raw/*_1.fastq
do
  samp=$(basename ${i} _1.fastq)

  echo "Trimming ${samp}"

  java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
  raw/${samp}_1.fastq raw/${samp}_2.fastq \
  trimmed/${samp}_1_paired.fastq trimmed/${samp}_1_unpaired.fastq \
  trimmed/${samp}_2_paired.fastq trimmed/${samp}_2_unpaired.fastq \
  ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/NexteraPE-PE.fa:2:40:15 \
  -threads $SLURM_CPUS_PER_TASK
  LEADING:3 \
  TRAILING:3 \
  MINLEN:36 \
  SLIDINGWINDOW:4:20
  

done
