#!/bin/bash

## =============================================================================
## Environment setup for the Physics Benchmarks CI jobs
##
## This script defines the following environment variables:
##
##  - LOCAL_DATA_PATH:        local storage for passing data between pipeline
##                            jobs (not saved as artifacts)
##
## Everything describing *what* the benchmarks compute -- the detector, its
## configuration, event counts, thread counts and RNG seeds -- lives in
## snakemake.yml instead, and is read by the workflow via Snakemake's config.
## Override those with "snakemake --config KEY=VALUE", not with environment
## variables.
## =============================================================================

echo "Setting up the Physics Benchmarks environment"

## Location of local data for passing data from job to job within a pipeline.
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

# Local field maps
mkdir -p ${LOCAL_DATA_PATH}/fieldmaps
ln -sf ${LOCAL_DATA_PATH}/fieldmaps

## =============================================================================
## That's all!
echo "Environment setup complete."
