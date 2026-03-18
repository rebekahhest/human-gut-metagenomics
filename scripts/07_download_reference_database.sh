#!/bin/bash
#SBATCH --job-name=download_reference_database
#SBATCH --account=def-cottenie
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

# Download Reference Database
mkdir -p kraken2_db
cd kraken2_db

wget https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz

tar -xvf k2_core_nt_20251015.tar.gz
  