#!/bin/bash
#SBATCH --job-name=taxonomic_classification
#SBATCH --account=def-cottenie
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

# Abundance Re-Estimation
module load bracken

cd data
mkdir -p bracken

for i in kraken2/*.report
do
  samp=$(basename ${i} .report)

  echo "Running Bracken on ${samp}" 

  bracken \
  -d kraken2_db \
  -i ${i} \
  -r 150 \
  -l S \
  -o bracken/${samp}.bracken \
  -w bracken/${samp}.report 

done