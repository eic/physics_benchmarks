#!/bin/bash

## =============================================================================
## Global configuration variables for the benchmark scripts
##
## This script is responsible for CI/local runner infrastructure defaults such as
## local data storage and PATH setup. Benchmark tuning values live in snakemake.yml
## and are applied by Snakemake's shell prefix.
##
## It defines the following additional variables for internal usage:
##  - LOCAL_PREFIX:           prefix for packages installed during the benchmark
##  - LOCAL_DATA_PATH:        local storage for pipeline jobs
##  - DETECTOR_PATH:          root path to locally installed detector definition xml files
##
## Finally, it makes sure LOCAL_PREFIX and JUGGLER_PREFIX are added to PATH
## and LD_LIBRARY_PATH
## =============================================================================

echo "Setting up the Physics Benchmarks environment"

## Location of local data for pass data from job to job within pipeline.
## Not saved as artifacts.
## Local /scratch directory is presumed to be writable. 
if [ ! -n  "${LOCAL_DATA_PATH}" ] ; then
  if [ -w /scratch ] ; then
    export LOCAL_DATA_PATH="/scratch/${CI_PROJECT_NAME}_${CI_PIPELINE_ID}"
  else
    echo "/scratch not writable; using $PWD/scratch"
    export LOCAL_DATA_PATH="$PWD/scratch/${CI_PROJECT_NAME}_${CI_PIPELINE_ID}"
  fi
fi
mkdir -p "${LOCAL_DATA_PATH}"
if [ ! -d "${LOCAL_DATA_PATH}" ]; then 
  echo "LOCAL_DATA_PATH (${LOCAL_DATA_PATH}) does not exist!!"
  echo "Creating LOCAL_DATA_PATH=$(pwd)/local_data "
  export LOCAL_DATA_PATH="$(pwd)/local_data"
  mkdir -p "${LOCAL_DATA_PATH}"
fi

## =============================================================================
## Other utility variables that govern how some of the dependent packages
## are built and installed. You should not have to change these.

## local prefix to be used for local storage of packages
## downloaded/installed during the benchmark process
LOCAL_PREFIX=".local"
mkdir -p "${LOCAL_PREFIX}"
export LOCAL_PREFIX=`realpath ${LOCAL_PREFIX}`

## detector path: root path to locally installed detector definition xml files
## The detector installation itself is configured by the Snakemake shell prefix
## using DETECTOR_PREFIX/thisepic.sh.
export DETECTOR_PATH="${LOCAL_PREFIX}/share/${DETECTOR:-epic}"

## build dir for ROOT to put its binaries etc.
export ROOT_BUILD_DIR=$LOCAL_PREFIX/root_build

## =============================================================================
## Setup PATH and LD_LIBRARY_PATH to include our prefixes
echo "Adding LOCAL_PREFIX and repo bin/ to PATH and LD_LIBRARY_PATH"
## Determine the directory of this script so we can add the repo's bin/ and include/ to PATH
_ENV_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_REPO_ROOT="$(dirname "${_ENV_SCRIPT_DIR}")"
export PATH=${LOCAL_PREFIX}/bin:${_ENV_SCRIPT_DIR}:${PATH}
export LD_LIBRARY_PATH=${LOCAL_PREFIX}/lib:${LD_LIBRARY_PATH}

## Include the repo's own include/ dir so '#include "common_bench/..."' resolves
## regardless of the working directory when ROOT is invoked
export ROOT_INCLUDE_PATH=${_REPO_ROOT}/include:${LOCAL_PREFIX}/include:${ROOT_INCLUDE_PATH}

# Local field maps
mkdir -p ${LOCAL_DATA_PATH}/fieldmaps
ln -sf ${LOCAL_DATA_PATH}/fieldmaps

## =============================================================================
## That's all!
echo "Environment setup complete."
