#!/bin/bash

source /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda/bin/activate

GRAAFDIR="../../experimental_data/graaf"
GRABOWDIR="../../experimental_data/grabow"
python generate_yaml.py $GRAAFDIR $GRABOWDIR

YAMLFILE="../../experimental_data/all_experiments.yaml"
python preprocess.py $YAMLFILE "sbr"
