#!/bin/bash
#SBATCH --job-name=rxn236
#SBATCH --error=error.log
#SBATCH --nodes=1
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --mincpus=32


python /home/harris.se/rmg/rmg_tools/autotst/kinetics/run_kinetics.py 236


