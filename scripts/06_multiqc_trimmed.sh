#!/bin/bash
#SBATCH --job-name=QC_report_post-trimming
#SBATCH --account=def-cottenie
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

# Generate Quality Report using MultiQC
module load multiqc

cd data

multiqc fastqc_trimmed -o fastqc_trimmed
