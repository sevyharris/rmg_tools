import os
import job_manager
import time



working_dir = '/scratch/harris.se/guassian_scratch/bac'




N = 10  # num runners
start_dir = os.getcwd()
currently_running = []

for i in range(20, 421):
    script_file = os.path.join(working_dir, f'species_{i:04}', 'run_orca.sh')
    while len(currently_running) >= N:
        time.sleep(10)
        jobs_to_remove = []
        for j in range(len(currently_running)):
            if currently_running[j].completed() or currently_running[j].failed():
                jobs_to_remove.append(currently_running[j])
        for job in jobs_to_remove:
            currently_running.remove(job)

    # add one job
    orca_job = job_manager.SlurmJob()
    slurm_cmd = f"sbatch {script_file}"
    os.chdir(os.path.dirname(script_file))
    orca_job.submit(slurm_cmd)
    print(f'submitted job')
    time.sleep(4)
    print(f'job ID is {orca_job.job_id}')
    currently_running.append(orca_job)
    print(f'{len(currently_running)} jobs currently running')
    os.chdir(start_dir)
    print()

