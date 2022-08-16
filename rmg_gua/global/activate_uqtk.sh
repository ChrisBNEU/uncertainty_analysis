#!/bin/bash
# load in the initialization script
source /work/westgroup/ChrisB/_01_MeOH_repos/UQTK_repos/path_config.sh

# setup UQTk install path as environment variable
export UQTK_INS=/work/westgroup/ChrisB/_01_MeOH_repos/UQTK_repos/UQTk-install
export PATH=${UQTK_INS}/bin:$PATH
