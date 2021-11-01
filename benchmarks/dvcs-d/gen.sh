#!/bin/bash

## Kong's dummy generator step. The generator files are all produced 
## Here we run this gen.sh step is to match the tags,paths,and formats for later steps.

## =============================================================================
## Standin for a proper pythia generation process, similar to how we
## generate events for DVMP
## Runs in 5 steps:
##   1. Parse the command line and setup the environment
##   2. Check if we can load the requested file from the cache
##   3. Build generator exe 
##   4. Run the actual generator
##   5. Finalize
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
source benchmarks/dvcs-d/env.sh

## Get a unique file name prefix based on the configuration options
GEN_TAG=gen-${CONFIG}_${JUGGLER_N_EVENTS} ## Generic file prefix

## =============================================================================
## Step 2: Check if we really need to run, or can use the cache.
if [ -f "${INPUT_PATH}/${GEN_TAG}.hepmc" ]; then
  echo "Found cached generator output for $GEN_TAG, no need to look for where the input file"
  exit 0
fi

echo "Generator output for $GEN_TAG not found in cache, need to copy generator files over"

## =============================================================================
## Step 3: skip
##         

## =============================================================================
## Step 4: Copy the event generator file over

# echo "copying from /gpfs02/eic/DATA/EpIC/DVCS_10x100_1M/ ..."
# echo "Copying the generator file. Warning: this is a local BNL directory  "
# cp /gpfs02/eic/DATA/EpIC/DVCS_10x100_1M/DVCS.1.hepmc ${TMP_PATH}/${GEN_TAG}.hepmc
# if [[ "$?" -ne "0" ]] ; then
#   echo "ERROR copying EpIC"
#   exit 1
# fi

DATA_URL="S3/eictest/ATHENA/EVGEN/EXCLUSIVE/DVCSD_ABCONV/10x100/DVCS-eD-ab-hiAcc-10x100.hepmc"
mc -C . config host add S3 https://dtn01.sdcc.bnl.gov:9000 $S3_ACCESS_KEY $S3_SECRET_KEY
mc -C . cat  --insecure ${DATA_URL} > "${TMP_PATH}/${GEN_TAG}.hepmc"
if [[ "$?" -ne "0" ]] ; then
  echo "Failed to download hepmc file"
  exit 1
fi

## =============================================================================
## Step 5: Finally, move relevant output into the artifacts directory and clean up
## =============================================================================
echo "Moving generator output into ${INPUT_PATH}"
mv ${TMP_PATH}/${GEN_TAG}.hepmc ${INPUT_PATH}/${GEN_TAG}.hepmc
## this step only matters for local execution
echo "Cleaning up"
## does nothing

## =============================================================================
## All done!
echo "$BENCHMARK_TAG event generation complete"
