#!/bin/bash
#SBATCH --job-name=monte_carlo
#SBATCH --error=error.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=0-9%50


start_index=$(( 1000 * $SLURM_ARRAY_TASK_ID ))

python /home/harris.se/rmg/rmg_tools/monte_carlo_sim/run_all_sims.py $start_index
