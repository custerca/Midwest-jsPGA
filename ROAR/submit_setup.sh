#!/bin/bash

#SBATCH --account=TXW19
#SBATCH --partition=sla-prio
#SBATCH --time=1:00:00
#SBATCH --job-name=jsPGA_sim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1


cd $SLURM_SUBMIT_DIR

module load r/4.3.1 

Rscript jspga_setup.R