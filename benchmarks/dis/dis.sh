#!/bin/bash
source strict-mode.sh

## =============================================================================
## Run the DIS benchmarks in 5 steps:
## 1. Parse the command line and setup environment
## 2. Detector simulation through ddsim
## 3. Digitization and reconstruction through Juggler
## 4. Root-based Physics analyses
## 5. Finalize
## =============================================================================

## make sure we launch this script from the project root directory
PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/../..
pushd ${PROJECT_ROOT}

echo "Running the DIS benchmarks"

## =============================================================================
## Step 1: Setup the environment variables
##
## First parse the command line flags.
## This sets the following environment variables:
## - CONFIG:   The specific generator configuration
## - EBEAM:    The electron beam energy
## - PBEAM:    The ion beam energy
export REQUIRE_MINQ2=true
source ${LOCAL_PREFIX}/bin/parse_cmd.sh $@

## To run the reconstruction, we need the following global variables:
## - JUGGLER_INSTALL_PREFIX: Install prefix for Juggler (simu/recon)
## - DETECTOR:       the detector package we want to use for this benchmark
## - DETECTOR_PATH:          full path to the detector definitions
##
## defined in common_bench repo
## You can ready bin/env.sh for more in-depth explanations of the variables
## and how they can be controlled.

## We also need the following benchmark-specific variables:
##
## - BENCHMARK_TAG: Unique identified for this benchmark process.
## - BEAM_TAG:      Identifier for the chosen beam configuration
## - INPUT_PATH:    Path for generator-level input to the benchmarks
## - TMP_PATH:      Path for temporary data (not exported as artifacts)
## - RESULTS_PATH:  Path for benchmark output figures and files
##
## You can read dvmp/env.sh for more in-depth explanations of the variables.
source benchmarks/dis/env.sh

## Get a unique file names based on the configuration options
GEN_FILE=${INPUT_PATH}/gen-${CONFIG}_${JUGGLER_N_EVENTS}.hepmc

SIM_FILE=${TMP_PATH}/sim-${CONFIG}.edm4hep.root
SIM_LOG=${TMP_PATH}/sim-${CONFIG}.log

# JUGGLER_REC_FILE_BASE= ${TMP_PATH}/rec-${CONFIG}
REC_FILE=${TMP_PATH}/rec-${CONFIG}.root
REC_LOG=${TMP_PATH}/sim-${CONFIG}.log

PLOT_TAG=${CONFIG}

## =============================================================================
## Step 2: Run the simulation
echo "Running Geant4 simulation"
if [ ! -f ${SIM_FILE} ] ; then
ddsim --runType batch \
      --part.minimalKineticEnergy 1000*GeV  \
      --filter.tracker edep0 \
      -v WARNING \
      --numberOfEvents ${JUGGLER_N_EVENTS} \
      --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml \
      --inputFiles ${GEN_FILE} \
      --outputFile ${SIM_FILE}
if [ "$?" -ne "0" ] ; then
  echo "ERROR running ddsim"
  exit 1
fi
fi

## =============================================================================
## Step 3: Run digitization & reconstruction
echo "Running the digitization and reconstruction"
if [ ${RECO} == "eicrecon" ] ; then
  /usr/bin/time -v eicrecon ${SIM_FILE} -Ppodio:output_file=${REC_FILE}
  if [ "$?" -ne "0" ] ; then
    echo "ERROR running eicrecon"
    exit 1
  fi
fi

if [[ ${RECO} == "juggler" ]] ; then
  export JUGGLER_SIM_FILE=${SIM_FILE}
  export JUGGLER_REC_FILE=${REC_FILE}
  gaudirun.py options/reconstruction.py || [ $? -eq 4 ]
  if [ "$?" -ne "0" ] ; then
    echo "ERROR running juggler"
    exit 1
  fi
fi


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
  "ebeam": ${EBEAM},
  "pbeam": ${PBEAM},
  "minq2": ${MINQ2},
  "test_tag": "${BEAM_TAG}"
}
EOF

root -b -q "benchmarks/dis/analysis/dis_electrons.cxx+g(\"${CONFIG}\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running dis_electron script"
  exit 1
fi

python benchmarks/dis/analysis/kinematics_correlations.py --rec_file ${REC_FILE} --config ${PLOT_TAG}_${DETECTOR_CONFIG} --results_path ${RESULTS_PATH} --nevents ${JUGGLER_N_EVENTS}
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running kinematics_correlations script"
  exit 1
fi

python benchmarks/dis/analysis/truth_reconstruction.py --rec_file ${REC_FILE} --config ${PLOT_TAG}_${DETECTOR_CONFIG} --results_path ${RESULTS_PATH} --nevents ${JUGGLER_N_EVENTS}
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running truth_reconstruction script"
  exit 1
fi


## =============================================================================
## Step 5: finalize
echo "Finalizing DIS benchmark"

## Move over reconsturction artifacts as long as we don't have
## too many events
if [ "${JUGGLER_N_EVENTS}" -lt "500" ] ; then 
  cp ${REC_FILE} ${RESULTS_PATH}
fi

## =============================================================================
## All done!
echo "DIS benchmarks complete"
