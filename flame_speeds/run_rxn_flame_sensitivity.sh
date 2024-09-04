#!/bin/bash
#SBATCH --job-name=rxn_flames
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=26
#SBATCH --array=0-70

# the upper limit of the array * 50 should be sized to the number of reactions.
# for example, if there are 3481 reactions, do ~70 tasks

# cpus-per-task is also sized such that the array of 51 temperatures can be done in two cycles

# this 50 must match the 50 in the for loop (line 229) of make_rxn_delay_npys.py
start_index=$(($SLURM_ARRAY_TASK_ID * 50))

python /work/westgroup/harris.se/autoscience/reaction_calculator/fs_sensitivity/make_rxn_fs_npys.py $1 $start_index
