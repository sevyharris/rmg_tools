#!/bin/bash
#SBATCH --job-name=kinetics_gen
#SBATCH --error=error.log
#SBATCH --nodes=1
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --mincpus=16
#SBATCH --partition=west


python /home/harris.se/rmg/rmg_tools/autotst/kinetics/reaction_kinetics.py


