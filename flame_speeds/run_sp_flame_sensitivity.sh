#!/bin/bash
#SBATCH --job-name=sp_flames
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=26

python /work/westgroup/harris.se/autoscience/reaction_calculator/fs_sensitivity/make_sp_fs_npys.py $1
