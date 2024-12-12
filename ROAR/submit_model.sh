#!/bin/bash

#SBATCH --account=TXW19
#SBATCH --partition=sla-prio
#SBATCH --time=21-00:00:00
#SBATCH --job-name=jsPGA_sim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100
#SBATCH --mem-per-cpu=4gb


cd $SLURM_SUBMIT_DIR

module load r/4.3.1 

Rscript jspga_model.R

