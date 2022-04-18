#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=CO_oxidation_Pt111
#SBATCH --partition=west,short
#SBATCH --mem=20Gb
#SBATCH --ntasks=1


python /home/harris.se/rmg/RMG-Py/rmg.py input.py



