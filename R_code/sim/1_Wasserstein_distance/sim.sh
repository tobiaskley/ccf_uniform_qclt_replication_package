#!/bin/bash 
#SBATCH --job-name=W1
#SBATCH --array 1-600
#SBATCH --output out/W1-%A_%a.out
#SBATCH --error out/W1-%A_%a.err
#SBATCH --time=0-6:00:00
#SBATCH --mem=4G
#SBATCH --partition=medium

module purge
module load gcc/14.2.0
module load openmpi/4.1.7
module load r/4.5.2

R CMD BATCH --vanilla "sim.R" "out/$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID.Routput"

