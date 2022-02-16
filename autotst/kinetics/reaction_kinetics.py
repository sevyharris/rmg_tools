import os

import rmgpy.species
import rmgpy.chemkin
import autotst.species
import autotst.reaction
import autotst.calculator.gaussian
import logging
from hotbit import Hotbit
import ase.io
import glob
import job_manager

def get_label(rmg_rxn):
    label = ''
    for reactant in rmg_rxn.reactants:
        label += f'{reactant.smiles}+'
    label = label[:-1]
    label += '_'
    for product in rmg_rxn.products:
        label += f'{product.smiles}+'
    label = label[:-1]
    return label


print("START!!!!")

# load the model
chemkin_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/chem_annotated.inp"
dictionary_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/species_dictionary.txt"
transport_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/tran.dat"
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(
    chemkin_path,
    dictionary_path=dictionary_path,
    transport_path=transport_path
)
print(f"Loaded model with {len(species_list)} species and {len(reaction_list)} reactions")


rxn_idx = 186  # R recombo
rxn_idx = 194
rxn_idx = 274 # first H abstraction

# rxn_idx = 236 another H abstraction




rmg_rxn = reaction_list[rxn_idx]
print(rmg_rxn)
print(rmg_rxn.family)

kinetics_base_dir = '/work/westgroup/harris.se/autoscience/dft/kinetics'
reaction_base_dir = os.path.join(kinetics_base_dir, f'reaction_{rxn_idx:04}')
os.makedirs(kinetics_base_dir, exist_ok=True)


label = get_label(rmg_rxn)
print("label", label)

reaction = autotst.reaction.Reaction(label=label, rmg_reaction=rmg_rxn)
print("reaction", reaction)

reaction.get_labeled_reaction()

print("Got labeled reaction")
transition_states = reaction.ts["reverse"]


print("ts0", transition_states[0])

print("rxn.ts", reaction.ts)
print("About to start hotbit")
logging.warning("Danger zone 0")
# reaction.generate_conformers_all(ase_calculator=Hotbit())
reaction.generate_conformers(ase_calculator=Hotbit())

print("ran hotbit")

for ts in reaction.ts['forward']:
    print(ts)


# overall calc
ts_dir = os.path.join(reaction_base_dir, 'ts')
# write Gaussian input files
for i, ts in enumerate(reaction.ts['forward']):
    gaussian = autotst.calculator.gaussian.Gaussian(conformer=ts)
    calc = gaussian.get_overall_calc()
    calc.label = f'fwd_ts_{i:04}'

    calc.directory = ts_dir
    calc.parameters.pop('scratch')
    calc.parameters.pop('multiplicity')
    calc.parameters['mult'] = ts.rmg_molecule.multiplicity
    
    calc.write_input(ts.ase_molecule)

print("done")


n_ts = len(reaction.ts['forward'])
print(f'{n_ts} found with Hotbit')

# Make a file to run Gaussian
slurm_run_file = os.path.join(ts_dir, 'run.sh')
slurm_settings = {
    '--job-name': f'g16_ts_{rxn_idx}',
    '--error': 'error.log',
    '--output': 'output.log',
    '--nodes': 1,
    '--partition': 'west,short',
    '--mem': '20Gb',
    '--time': '24:00:00',
    '--cpus-per-task': 16,
    '--array': f'0-{n_conformers}%40',
}

slurm_file_writer = job_manager.SlurmJobFile(
    full_path=slurm_run_file,
)
slurm_file_writer.settings = slurm_settings

slurm_file_writer.content = [
    'export GAUSS_SCRDIR=/scratch/harris.se/guassian_scratch\n',
    'mkdir -p $GAUSS_SCRDIR\n',
    'module load gaussian/g16\n',
    'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

    'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
    'fname="fwd_ts_${RUN_i}.com"\n\n',

    'g16 $fname\n',
]
slurm_file_writer.write_file()

# submit the job
gaussian_conformers_job = job_manager.SlurmJob()
slurm_cmd = f"sbatch {slurm_run_file}"
gaussian_conformers_job.submit(my_cmd)






