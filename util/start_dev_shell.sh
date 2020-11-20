#!/bin/bash

## =============================================================================
## Setup (if needed) and start a development shell environment on Linux or MacOS
## =============================================================================

## make sure we launch this script from the project root directory
PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/..
pushd ${PROJECT_ROOT}

## We do not load the global development environment here, as this script is 
## to be executed on a "naked" system outside of any container 

## =============================================================================
## Step 1: Parse command line options

## do we want to force-update the container (only affects Linux)
## default: we do not want to do this.
FORCE_UPDATE=

function print_the_help {
  echo "USAGE:    ./util/start_dev_shell [-f]"
  echo "OPTIONS:"
  echo "          -f,--force    Force-update container (Only affects Linux)"
  echo "          -h,--help     Print this message"
  echo ""
  echo "  This script will setup and launch a containerized development
  environment"
  exit
}
while [ $# -gt 0 ]
do
  key="$1"
  case $key in
    -f|--force)
      FORCE_UPDATE="true"
      shift # past value
      ;;
    -h|--help)
      print_the_help
      shift
      ;;
    *)    # unknown option
      echo "unknown option $1"
      exit 1
      ;;
  esac
done

## get OS type
OS=`uname -s`

## =============================================================================
## Step 2: Update container and launch shell
echo "Launching a containerized development shell"

case ${OS} in
  Linux)
    echo "  - Detected OS: Linux"
    ## Use the same prefix as we use for other local packages
    export PREFIX=.local/lib
    if [ ! -f $PREFIX/juggler_latest.sif ] || [ ! -z ${FORCE_UPDATE} ]; then
      echo "  - Fetching singularity image"
      mkdir -p $PREFIX
      wget https://eicweb.phy.anl.gov/eic/juggler/-/jobs/artifacts/master/raw/build/juggler.sif?job=docker:singularity
      -O $PREFIX/juggler_latest.sif
    fi
    echo "  - Using singularity to launch shell..."
    singularity exec $PREFIX/juggler_latest.sif eic-shell
    ;;
  Darwin)
    echo "  - Detector OS: MacOS"
    echo "  - Syncing docker container"
    docker pull sly2j/juggler:latest
    echo "  - Using docker to launch shell..."
    docker run -v /Users:/Users -w=$PWD -i -t --rm sly2j/juggler:latest eic-shell
    ;;
  *)
    echo "ERROR: dev shell not available for this OS (${OS})"
    exit 1
esac

## =============================================================================
## Step 3: All done
echo "Exiting development environment..."
