import job_manager
import time
import glob
import os


# working_dir = os.path.join(os.getcwd(),"uncertainty_output_folder/rmg_run_scripts/")
working_dir = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/rmg_run_scripts/"
print("Collecting SLURM scripts")
slurm_scripts = glob.glob(os.path.join(working_dir, "rmg_runs_*.sh"))

slurm_scripts.sort()

for i, script in enumerate(slurm_scripts):
    print(f"{i}/{len(slurm_scripts)}\tRunning job {script}")

    rmg_job = job_manager.SlurmJob()
    my_cmd = f'sbatch {script}'
    print(my_cmd)
    rmg_job.submit(my_cmd)

    # wait for job
    # rmg_job.wait()
    rmg_job.wait_all()
    print(i)




