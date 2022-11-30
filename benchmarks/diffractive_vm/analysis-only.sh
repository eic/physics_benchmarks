#!/bin/bash

## make sure we launch this script from the project root directory
PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/../..
pushd ${PROJECT_ROOT}

echo "Running the diffractive_vm benchmarks"
source parse_cmd.sh $@
source benchmarks/diffractive_vm/env.sh

GEN_FILE=${INPUT_PATH}/gen-${CONFIG}_${JUGGLER_N_EVENTS}.hepmc

SIM_FILE=${TMP_PATH}/sim-${CONFIG}.edm4hep.root
SIM_LOG=${TMP_PATH}/sim-${CONFIG}.log


REC_FILE=${TMP_PATH}/rec-${CONFIG}.root
REC_LOG=${TMP_PATH}/sim-${CONFIG}.log

PLOT_TAG=${CONFIG}

## =============================================================================
## Step 4: Analysis
## write a temporary configuration file for the analysis script
echo "Running analysis"
CONFIG="${TMP_PATH}/${PLOT_TAG}.json"
cat << EOF > ${CONFIG}
{
  "rec_file": "${REC_FILE}",
  "detector": "${DETECTOR}",
  "output_prefix": "${RESULTS_PATH}/${PLOT_TAG}",
  "test_tag": "${BEAM_TAG}"
}
EOF
#cat ${CONFIG}
export VM_TYPE_TAG=1
export MC_TYPE_TAG=1
root -b -q "benchmarks/diffractive_vm/analysis/diffractive_vm_analysis.cxx+(\"${CONFIG}\",${VM_TYPE_TAG},${MC_TYPE_TAG})"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running rec_diffractive_vm_analysis script"
  exit 1
fi

echo "Diffractive VM benchmarks analysis complete"