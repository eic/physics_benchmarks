#!/bin/bash
source strict-mode.sh

## =============================================================================
## Run the Diffractive VMP benchmarks in 5 steps:
## 1. Parse the command line and setup environment
## 2. Detector simulation through ddsim
## 3. Digitization and reconstruction through Juggler
## 4. Root-based Physics analyses
## 5. Finalize
## =============================================================================

## make sure we launch this script from the project root directory
PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/../..
pushd ${PROJECT_ROOT}

echo "Running the Diffractive VMP benchmarks"


## =============================================================================
## Step 1: Setup the environment variables
##
## First parse the command line flags.
## This sets the following environment variables:
## - CONFIG:   The specific generator configuration
## - EBEAM:    The electron beam energy
## - PBEAM:    The ion beam energy
## - LEADING:  Leading particle of interest (J/psi)
export REQUIRE_LEADING=1
source parse_cmd.sh $@

## We also need the following benchmark-specific variables:
##
## - BENCHMARK_TAG: Unique identified for this benchmark process.
## - BEAM_TAG:      Identifier for the chosen beam configuration
## - INPUT_PATH:    Path for generator-level input to the benchmarks
## - TMP_PATH:      Path for temporary data (not exported as artifacts)
## - RESULTS_PATH:  Path for benchmark output figures and files
##
## You can read dvmp/env.sh for more in-depth explanations of the variables.
source benchmarks/diffractive_vm/env.sh

## Get a unique file names based on the configuration options
GEN_FILE=${INPUT_PATH}/gen-${CONFIG}_${LEADING}_${JUGGLER_N_EVENTS}.hepmc

SIM_FILE=${TMP_PATH}/sim-${CONFIG}.edm4hep.root
SIM_LOG=${TMP_PATH}/sim-${CONFIG}.log


REC_FILE=${TMP_PATH}/rec-${CONFIG}.root
REC_LOG=${TMP_PATH}/sim-${CONFIG}.log

PLOT_TAG=${CONFIG}

## =============================================================================
## Step 2: Run the simulation
echo "Running Geant4 simulation"
ls -lrth 
ls -lrth input
echo ${TMP_PATH}
ls -lrth ${TMP_PATH}
ddsim --runType batch \
      --part.minimalKineticEnergy 100*GeV  \
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

## =============================================================================
## Step 3: Run digitization & reconstruction
echo "Running the digitization and reconstruction"
## FIXME Need to figure out how to pass file name to juggler from the commandline
## the tracker_reconstruction.py options file uses the following environment
## variables:
## - JUGGLER_SIM_FILE:    input detector simulation
## - JUGGLER_REC_FILE:    output reconstructed data
## - JUGGLER_N_EVENTS:    number of events to process (part of global environment)
## - DETECTOR:    detector package (part of global environment)
export JUGGLER_SIM_FILE=${SIM_FILE}
export JUGGLER_REC_FILE=${REC_FILE}
if [ ${RECO} == "eicrecon" ] ; then
  eicrecon ${JUGGLER_SIM_FILE} -Ppodio:output_file=${JUGGLER_REC_FILE}
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running eicrecon"
    exit 1
  fi
fi

if [[ ${RECO} == "juggler" ]] ; then
  gaudirun.py options/reconstruction.py || [ $? -eq 4 ]
  if [ "$?" -ne "0" ] ; then
    echo "ERROR running juggler"
    exit 1
  fi
fi
## =============================================================================
## Step 4: Analysis

## write a temporary configuration file for the analysis script
CONFIG="${TMP_PATH}/${PLOT_TAG}.json"
cat << EOF > ${CONFIG}
{
  "rec_file": "${REC_FILE}",
  "vm_name": "${LEADING}",
  "detector": "${DETECTOR_CONFIG}",
  "ebeam": ${EBEAM},
  "pbeam": ${PBEAM},
  "output_prefix": "${RESULTS_PATH}/${PLOT_TAG}",
  "test_tag": "${LEADING}_${BEAM_TAG}"
}
EOF
#cat ${CONFIG}

## run the analysis script with this configuration
root -b -q "benchmarks/diffractive_vm/analysis/diffractive_vm.cxx+(\"${CONFIG}\")"
if [ "$?" -ne "0" ] ; then
  echo "ERROR running vm_mass script"
  exit 1
fi

## =============================================================================
## Step 5: finalize
echo "Finalizing Diffractive VMP benchmark"

## Copy over reconstruction artifacts as long as we don't have
## too many events
if [ "${JUGGLER_N_EVENTS}" -lt "500" ] ; then 
  cp ${REC_FILE} ${RESULTS_PATH}
fi

## cleanup output files
#rm -f ${REC_FILE} ${SIM_FILE} ## --> not needed for CI

## =============================================================================
## All done!
echo "Diffractive VMP benchmarks complete"
