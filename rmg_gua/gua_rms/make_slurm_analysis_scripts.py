# script to do the 5000 RMG runs on discovery
# Sevy Harris 9/8/2021

import numpy as np
import job_manager
import os
import glob

def make_slurm_analysis_scripts(
    unc_folder, 
    working_dir, 
    conda_path,
    output_name="rms_analysis.csv", 
    M=10, 
    N=50,
    ): 
    # WARNING - this will fail if M%N != 0
    skip_completed_runs = False  # set to false to overwrite RMG runs that completed
    remove_old_slurm_scripts = True

    # check if "rms_run_scripts" folder exists
    slurm_script_dir = os.path.join(working_dir, "rms_run_scripts")
    if not os.path.exists(slurm_script_dir):
        os.mkdir(slurm_script_dir)



    # if remove old slurm scripts is true, will clear directory before 
    # proceeding
    if remove_old_slurm_scripts:
        rm_rms_files = glob.glob(os.path.join(working_dir, "rms_run_scripts/rms_runs_*.sh"))
        for file in rm_rms_files:
            try:
                os.remove(file)
            except OSError:
                print(f"file doesn't exist at {file}")
    
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    # M = 5000  # total number of times to run RMG
    # N = 50  # number of jobs to run at a time
    for i in range(0, M, N):
        sbatch_index = int(i / N)
        range_max = np.amin([i + N, M])
        last_index = range_max - 1
        job_indices = [a for a in range(i, range_max)]

        # max slurm array index is 1000, so after that, subtract multiples of 1000
        task_id_offset = int(i/1000) * 1000
        
        # Write the job file
        fname = f'rms_runs_{i}-{last_index}.sh'
        script_folder = "rms_run_scripts"
        jobfile = job_manager.SlurmJobFile(full_path=os.path.join(working_dir, script_folder, fname))
        jobfile.settings['--array'] = f'{i - task_id_offset}-{last_index - task_id_offset}'
        jobfile.settings['--job-name'] = fname
        jobfile.settings['--error'] = os.path.join(working_dir, script_folder, f'rms_job_error{sbatch_index}.log')
        jobfile.settings['--output'] = os.path.join(working_dir, script_folder, f'rms_job_output{sbatch_index}.log')
        jobfile.settings['--mem'] = f'20Gb'
        jobfile.settings['--cpus-per-task'] = '4'
        
        
        content = ['# activate conda environment\n']
        content.append(f'source activate {conda_path}\n')
        content.append('# Define useful bash variables\n')
        
        content.append(f'SLURM_TASK_ID_OFFSET={task_id_offset}\n')
        content.append('RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID + $SLURM_TASK_ID_OFFSET)))\n')
        rmg_run_dir = os.path.join(working_dir, "run_${RUN_i}")

        content.append(f'CSV_FILE="{rmg_run_dir}/rms/{output_name}"\n')

        # how do we get it to autodetect the latest RMS run? could we just have rmg do it? 
        # use ls | tail -1 to get the latest file
        content.append(f'MECH_FILE=$(ls "{rmg_run_dir}/rms/" | tail -1)\n')
        
        # skip if csv file already exists
        # this will not work properly if there is a file with a new name, revise
        # to remove any results file? 
        if skip_completed_runs:
            content.append('if test -f "$CSV_FILE"; then\n')
            content.append('echo "skipping completed run ${RUN_i}"; exit 0\n')
            content.append('fi\n\n')
        else: 
            # run the analysis script
            content.append('# remove the old CSV file\n')
            content.append('rm -f $CSV_FILE\n') 
            
        

        content.append('# Run the analysis\n')    
        content.append(f'python-jl {unc_folder + "gua_rms/run_reactor.py"} {rmg_run_dir}/rms/$MECH_FILE {output_name} \n')
        jobfile.content = content
        jobfile.write_file()
    
 
