#!/bin/bash
source strict-mode.sh

## =============================================================================
## Standin for a proper pythia generation process, similar to how we
## generate events for DVMP
## Runs in 5 steps:
##   1. Parse the command line and setup the environment
##   2. Check if we can download the file
##   3. Finalize
## =============================================================================
## =============================================================================

## make sure we launch this script from the project root directory
PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/../..
pushd ${PROJECT_ROOT}

## =============================================================================
## Step 1: Setup the environment variables
## First parse the command line flags.
## This sets the following environment variables:
## - CONFIG:   The specific generator configuration --> not currenlty used FIXME
## - EBEAM:    The electron beam energy --> not currently used FIXME
## - PBEAM:    The ion beam energy --> not currently used FIXME
export REQUIRE_LEADING=1
source parse_cmd.sh $@

## To run the generator, we need the following global variables:
##
## - LOCAL_PREFIX:      Place to cache local packages and data
## - JUGGLER_N_EVENTS:  Number of events to process
## - JUGGLER_RNG_SEED:  Random seed for event generation.
##
## defined in common_bench repo
## You can ready bin/env.sh for more in-depth explanations of the variables
## and how they can be controlled.


## We also need the following benchmark-specific variables:
##
## - BENCHMARK_TAG: Unique identified for this benchmark process.
## - INPUT_PATH:    Path for generator-level input to the benchmarks
## - TMP_PATH:      Path for temporary data (not exported as artifacts)
##
## You can read dvmp/env.sh for more in-depth explanations of the variables.
source benchmarks/diffractive_vm/env.sh

## Get a unique file name prefix based on the configuration options
GEN_TAG=gen-${CONFIG}_${LEADING}_${JUGGLER_N_EVENTS} ## Generic file prefix

## =============================================================================
## Step 2: Check if we can find the file
if [ -f "${INPUT_PATH}/${GEN_TAG}.hepmc" ]; then
  echo "Found cached generator output for $GEN_TAG, no need to rerun"
  exit 0
fi

## =============================================================================
## Step 3: Copy the file
events_per_file=100
nfiles=$(( (${JUGGLER_N_EVENTS} + ${events_per_file} - 1) / ${events_per_file} ))
for ix in $(seq -f "%03g" 0 10); do
  DATA_URL=S3/eictest/EPIC/EVGEN/EXCLUSIVE/DIFFRACTIVE_${LEADING^^}_ABCONV/Sartre/Coherent/sartre_bnonsat_Au_${LEADING}_ab_eAu_1_${ix}.hepmc.gz
  mc config host add S3 https://eics3.sdcc.bnl.gov:9000 ${S3_ACCESS_KEY} ${S3_SECRET_KEY}
  mc cp ${DATA_URL} ${TMP_PATH}/${GEN_TAG}.hepmc.gz
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR downloading file"
    exit 1
  fi
  gunzip -c ${TMP_PATH}/${GEN_TAG}.hepmc.gz >> ${TMP_PATH}/${GEN_TAG}.hepmc
done

## =============================================================================
## Step 4: Finally, move relevant output into the artifacts directory and clean up
## =============================================================================
echo "Moving generator output to ${INPUT_PATH}/${GEN_TAG}.hepmc"
mv ${TMP_PATH}/${GEN_TAG}.hepmc ${INPUT_PATH}/${GEN_TAG}.hepmc
## this step only matters for local execution
echo "Cleaning up"
## does nothing

## =============================================================================
## All done!
echo "$BENCHMARK_TAG event generation complete"
