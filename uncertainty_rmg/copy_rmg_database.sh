#!/bin/bash
#SBATCH --job-name=copy_rmg_database
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c5003
#SBATCH --mem=8Gb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --error=dberror.log


cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_rmg
source activate ../conda


RMG_db_folder="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/"
output_path="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/"
unc_folder="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
N=3 # db number

python-jl ./copy_rmg_database.py  ${RMG_db_folder} ${output_path} ${unc_folder} ${N}

# Define useful bash variables
RUN_i=$(printf "%04.0f" $(($N)))
DATABASE_n=$(printf "%04.0f" $(($(($N)) % 10)))
# Copy the files from the full database to the mostly symbolic one
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction_Beta/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction_Beta/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction_Beta_double_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction_Beta_double_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction_Beta_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction_Beta_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction_Single_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction_Single_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Abstraction_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Abstraction_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Addition_Single_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Addition_Single_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Abstraction_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Abstraction_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Bidentate/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Bidentate/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Dissociative/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Dissociative/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Dissociative_Double/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Dissociative_Double/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Double/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Double/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_Single/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_Single/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Adsorption_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Adsorption_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Bidentate_Dissociation/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Bidentate_Dissociation/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation_Beta/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation_Beta/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation_Beta_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation_Beta_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation_Double/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation_Double/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation_Double_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation_Double_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dissociation_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Dissociation_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_DoubleBond_to_Bidentate/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_DoubleBond_to_Bidentate/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Dual_Adsorption_vdW/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Dual_Adsorption_vdW/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_EleyRideal_Addition_Multiple_Bond/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_EleyRideal_Addition_Multiple_Bond/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_Migration/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_Migration/rules.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/families/Surface_vdW_to_Bidentate/rules_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/families/Surface_vdW_to_Bidentate/rules.py"

cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/thermo/groups/adsorptionPt111_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/thermo/groups/adsorptionPt111.py"

rm "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/thermo/libraries/surfaceThermoPt111.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/thermo/libraries/surfaceThermoPt111_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/thermo/libraries/surfaceThermoPt111.py"
rm "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/libraries/Surface/CPOX_Pt/Deutschmann2006_adjusted/reactions.py"
cp "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database/input/kinetics/libraries/Surface/CPOX_Pt/Deutschmann2006_adjusted/reactions_${RUN_i}.py" "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/db_${DATABASE_n}/input/kinetics/libraries/Surface/CPOX_Pt/Deutschmann2006_adjusted/reactions.py"




