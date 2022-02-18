#!/bin/bash
#SBATCH --job-name=arkane
#SBATCH --nodes=1
#SBATCH --array=0-10

indices=(7 17 58 68 69 71 73 77 79 80 84)

index=${indices[$SLURM_ARRAY_TASK_ID]}


RUN_i=$(printf "%04.0f" $index)

cd "/work/westgroup/harris.se/autoscience/dft/thermo/species_${RUN_i}/conformers/arkane/"
python /home/harris.se/rmg/RMG-Py/Arkane.py input.py

