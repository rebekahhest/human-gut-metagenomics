#!/bin/bash
#SBATCH --job-name=taxonomic_classification
#SBATCH --account=def-cottenie
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

# Taxonomic Classification
module load kraken2

cd data
mkdir -p kraken2

for i in trimmed/*_1_paired.fastq
do
  samp=$(basename ${i} _1_paired.fastq)

  echo "Running Kraken2 on ${samp}" 

  kraken2 \
  --db kraken2_db \
  --confidence 0.15 \
  --use-names \
  --threads $SLURM_CPUS_PER_TASK \
  --paired trimmed/${samp}_1_paired.fastq trimmed/${samp}_2_paired.fastq \
  --output kraken2/${samp}.kraken \
  --report kraken2/${samp}.report

done
