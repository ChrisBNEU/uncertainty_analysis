#!/bin/bash
#SBATCH --job-name=rmg_runs_0-99.sh
#SBATCH --error=/scratch/blais.ch/methanol_results_2022_04_06/error0.log
#SBATCH --output=/scratch/blais.ch/methanol_results_2022_04_06/output0.log
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=0-99


# Define useful bash variables
SLURM_TASK_ID_OFFSET=0
RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID + $SLURM_TASK_ID_OFFSET)))
DATABASE_n=$(printf "%04.0f" $(($(($SLURM_ARRAY_TASK_ID + $SLURM_TASK_ID_OFFSET)) % 100)))

python copy_rmg_database.py /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/ /work/westgroup/ChrisB/Scratch/  5
# Copy the files from the full database to the mostly symbolic one
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction_Beta/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction_Beta/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction_Beta_double_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction_Beta_double_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction_Beta_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction_Beta_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction_Single_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction_Single_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Addition_Single_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Addition_Single_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Abstraction_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Abstraction_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Bidentate/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Bidentate/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Dissociative/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Dissociative/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Dissociative_Double/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Dissociative_Double/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Double/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Double/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Single/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Single/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Bidentate_Dissociation/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Bidentate_Dissociation/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation_Beta/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation_Beta/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation_Beta_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation_Beta_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation_Double/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation_Double/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation_Double_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation_Double_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_DoubleBond_to_Bidentate/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_DoubleBond_to_Bidentate/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dual_Adsorption_vdW/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Dual_Adsorption_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_EleyRideal_Addition_Multiple_Bond/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_EleyRideal_Addition_Multiple_Bond/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Migration/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_Migration/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_vdW_to_Bidentate/rules_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/families/Surface_vdW_to_Bidentate/rules.py"

cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/thermo/groups/adsorptionPt111_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/thermo/groups/adsorptionPt111.py"

rm "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/thermo/libraries/surfaceThermoPt111.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/thermo/libraries/surfaceThermoPt111_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/thermo/libraries/surfaceThermoPt111.py"
rm "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/libraries/Surface/CPOX_Pt/Deutschmann2006_adjusted/reactions.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/libraries/Surface/CPOX_Pt/Deutschmann2006_adjusted/reactions_${RUN_i}.py" "/scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/kinetics/libraries/Surface/CPOX_Pt/Deutschmann2006_adjusted/reactions.py"

# Prepare the directory for the RMG run
mkdir "/scratch/blais.ch/methanol_results_2022_04_06/run_${RUN_i}"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_rmg/input.py" "/scratch/blais.ch/methanol_results_2022_04_06/run_${RUN_i}/input.py"
echo "database.directory /scratch/blais.ch/methanol_results_2022_04_06/db_${DATABASE_n}/input/" > "/scratch/blais.ch/methanol_results_2022_04_06/run_${RUN_i}/rmgrc"
# Run RMG
cd /scratch/blais.ch/methanol_results_2022_04_06/run_${RUN_i} 
export PYTHONPATH="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-Py/:$PYTHONPATH"
export PATH="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-Py/:$PATH"
export RMG="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-Py//rmg.py"
source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda/
python-jl $RMG input.py
