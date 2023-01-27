#!/bin/bash
source ~/.zshrc
conda activate rmg_env

GRAAFDIR="/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/experimental_data/graaf"
GRABOWDIR="/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/experimental_data/grabow"
python generate_yaml.py $GRAAFDIR $GRABOWDIR

YAMLFILE="/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/gua_cantera/all_experiments.yaml"
python preprocess.py $YAMLFILE "sbr"
