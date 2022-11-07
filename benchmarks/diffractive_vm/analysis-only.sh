#!/bin/bash

## make sure we launch this script from the project root directory
PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/../..
pushd ${PROJECT_ROOT}

echo "Running the diffractive_vm benchmarks"

## =============================================================================
## Step 4: Analysis
## write a temporary configuration file for the analysis script
echo "Running analysis"
CONFIG="${TMP_PATH}/${PLOT_TAG}.json"
cat << EOF > ${CONFIG}
{
  "rec_file": "${REC_FILE}",
  "detector": "${JUGGLER_DETECTOR}",
  "output_prefix": "${RESULTS_PATH}/${PLOT_TAG}",
  "ebeam": ${EBEAM},
  "pbeam": ${PBEAM},
  "minq2": ${MINQ2},
  "test_tag": "${BEAM_TAG}"
}
EOF
#cat ${CONFIG}
export VM_TYPE_TAG=2
export MC_TYPE_TAG=1
root -b -q "benchmarks/diffractive_vm/analysis/diffractive_vm_analysis.cxx+(\"${CONFIG}\",${VM_TYPE_TAG},${MC_TYPE_TAG})"
#root -b -q "benchmarks/diffractive_vm/analysis/test_phi_analysis.cxx+(\"${CONFIG}\",${VM_TYPE_TAG},${MC_TYPE_TAG})"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running rec_diffractive_vm_analysis script"
  exit 1
fi

echo "Diffractive Phi benchmarks analysis complete"
