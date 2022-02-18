#!/bin/bash
#SBATCH --job-name=arkane
#SBATCH --nodes=1

#indices=(7 17 58 68 71 73 77 79 80 84)


#for i in ${indices[@]}
#do
    
#    echo $i
#done


python /home/harris.se/rmg/rmg_tools/autotst/thermo/write_arkane.py

#python /home/harris.se/rmg/rmg_tools/autotst/thermo/run_thermo.py 71


