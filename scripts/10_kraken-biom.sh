#!/bin/bash
#SBATCH --job-name=create_biom_table
#SBATCH --account=def-cottenie
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

# Abundance Re-Estimation
#source kraken_env/bin/activate 

cd data/bracken

kraken-biom *.report \
-o bracken.biom \
--fmt json

