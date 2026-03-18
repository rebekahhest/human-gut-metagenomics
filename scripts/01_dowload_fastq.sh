#!/bin/bash
#SBATCH --job-name=download_SRR_files
#SBATCH --account=def-cottenie
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

# Download FASTQ files
module load sra-toolkit

mkdir -p data/raw
mkdir -p tmp
cd data/raw

for i in $(cat sra.txt); 
do
  echo "Downloading ${i}"

  prefetch ${i}

  fasterq-dump ${i} --split-files --threads $SLURM_CPUS_PER_TASK --mem 8G --temp tmp

  rm -r ${i}
  
done
cd ../../