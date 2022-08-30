#!/bin/bash -e


#=====================================================================================
# Need to have UQTK_INS defined and bin added to the path.
# export UQTK_INS=${HOME}/research/UQTk-install
# export PATH=${UQTK_INS}/bin:$PATH
#=====================================================================================

# Script location
export UQPC=${UQTK_INS}/examples/uqpc

source activate uqtk

# Plot model-vs-surrogate
${UQPC}/plot.py dm training validation

# Plot total sensitivities
${UQPC}/plot.py sens total

${UQPC}/plot.py sens main
# Enjoy the .eps files!

