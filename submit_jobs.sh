#!/bin/bash
#SBATCH -N 1
#SBATCH -t 00:59:00

module purge
module load 2020
module load R/4.0.2-intel-2020a

cp -r "$HOME"/GMM_ord "$TMPDIR"
cd "$TMPDIR"/GMM_ord

echo $SLURM_ARRAY_TASK_ID

Rscript --vanilla Simulation_LISA.R $SLURM_ARRAY_TASK_ID

cp -r ./*.RDS "$HOME"/GMM_ord

