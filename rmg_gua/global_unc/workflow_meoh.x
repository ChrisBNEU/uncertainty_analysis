#!/bin/bash -e


#=====================================================================================
# Need to have UQTK_INS defined and bin added to the path.
# export UQTK_INS=${HOME}/research/UQTk-install
# export PATH=${UQTK_INS}/bin:$PATH
#=====================================================================================

# Script location
export UQPC=${UQTK_INS}/examples/uqpc

source activate uqtk

# Use all inputs as training, and no validation
# cp Input.txt ptrain.dat
# #cp output_CH4.txt ytrain.dat # Or better work with the logarithms
# awk '{print log($1), log($2), log($3)}' output_CH4.txt > ytrain.dat

cp parameter_names_meoh.txt pnames.txt

TRAIN=3410
VAL=1137

# this needs to be extensible, so it just takes 10% of the runs for training
head -n$TRAIN Input_meoh.txt > ptrain.dat
tail -n$VAL Input_meoh.txt > pval.dat

head -n$TRAIN outputs_meoh.txt > ytrain.dat
tail -n$VAL outputs_meoh.txt > yval.dat

# Scale the inputs
${UQPC}/scale.x ptrain.dat from parameter_ranges_meoh.txt qtrain.dat
${UQPC}/scale.x pval.dat from parameter_ranges_meoh.txt qval.dat

# Get the number of training samples and validation samples
NSAM=`echo | awk 'END{print NR}' ptrain.dat`
NVAL=`echo | awk 'END{print NR}' pval.dat`

echo $NSAM
echo $NVAL

# Build surrogates
${UQPC}/uq_pc.py -r offline_post -p parameter_ranges_meoh.txt -m bcs -s rand -n $NSAM -v $NVAL -t 3 -e 1.e-7

# Plot model-vs-surrogate
${UQPC}/plot.py dm training validation

# Plot total sensitivities
${UQPC}/plot.py sens total

# Enjoy the .eps files!

