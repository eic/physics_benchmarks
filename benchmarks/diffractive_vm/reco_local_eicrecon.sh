#!/bin/bash

## =============================================================================
## Run the DIS benchmarks in 5 steps:
## 1. Parse the command line and setup environment
## 2. Detector simulation through npsim
## 3. Digitization and reconstruction through Juggler
## 4. Root-based Physics analyses
## 5. Finalize
## =============================================================================

## make sure we launch this script from the project root directory
PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/../..
pushd ${PROJECT_ROOT}

echo "Running the Diffractive VM benchmarks"

## =============================================================================
## Step 1: Setup the environment variables
##
## First parse the command line flags.
## This sets the following environment variables:
## - CONFIG:   The specific generator configuration
## - EBEAM:    The electron beam energy
## - PBEAM:    The ion beam energy
source ${LOCAL_PREFIX}/bin/parse_cmd.sh $@

## To run the reconstruction, we need the following global variables:
## - JUGGLER_INSTALL_PREFIX: Install prefix for Juggler (simu/recon)
## - JUGGLER_DETECTOR:       the detector package we want to use for this benchmark
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
source benchmarks/diffractive_vm/env.sh

## Get a unique file names based on the configuration options
GEN_FILE=${INPUT_PATH}/gen-${CONFIG}_${JUGGLER_N_EVENTS}.hepmc

SIM_FILE=${TMP_PATH}/sim-${CONFIG}.edm4hep.root
SIM_LOG=${TMP_PATH}/sim-${CONFIG}.log

REC_FILE=${TMP_PATH}/rec-${CONFIG}.root
REC_LOG=${TMP_PATH}/sim-${CONFIG}.log

PLOT_TAG=${CONFIG}

# Disable now for reco only for EICrecon
## =============================================================================
## Step 2: Run the simulation
# echo "Running Geant4 simulation"
# ddsim --runType batch \
#       --part.minimalKineticEnergy 1000*GeV  \
#       -v INFO \
#       --numberOfEvents ${JUGGLER_N_EVENTS} \
#       --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml \
#       --inputFiles ${GEN_FILE} \
#       --outputFile ${SIM_FILE}
# if [ "$?" -ne "0" ] ; then
#   echo "ERROR running ddsim"
#   exit 1
# fi

## =============================================================================
## Step 3: Run digitization & reconstruction
echo "Running the digitization and reconstruction"
run_eicrecon_reco_flags.py ${SIM_FILE} $(basename ${REC_FILE} .root)
if [ "$?" -ne "0" ] ; then
  echo "ERROR running eicrecon"
  exit 1
fi
mv $(basename ${REC_FILE} .root).tree.edm4eic.root ${TMP_PATH}/$(basename ${REC_FILE} .root).tree.edm4eic.root

echo "Diffractive VM benchmarks reco complete."

