#!/bin/bash
set -Eu
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR
IFS=$'\n\t'

function print_the_help {
  echo "USAGE: ${0} [--rec] [--sim] [--analysis] [--all] "
  echo "    The default options are to run all steps (sim,rec,analysis) "
  echo "OPTIONS: "
  echo "  --data-init     download the input event data"
  echo "  --sim,-s        Runs the Geant4 simulation"
  echo "  --rec,-r        Run the juggler reconstruction"
  echo "  --analysis,-a   Run the analysis scripts"
  echo "  --all           (default) Do all steps. Argument is included so usage can convey intent."
  exit 
}

DO_ALL=1
DATA_INIT=
DO_SIM=
DO_REC=
DO_ANALYSIS=
EBEAM=
PBEAM=
TAG=

POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    -h|--help)
      shift # past argument
      print_the_help
      ;;
    --all)
      DO_ALL=2
      if [[ ! "${DO_REC}${DO_SIM}${DO_ANALYSIS}" -eq "" ]] ; then
        echo "Error: cannot use --all with other arguments." 1>&2
        print_the_help
        exit 1
      fi
      shift # past value
      ;;
    --tag)
      shift # past argument
      TAG=$1
      shift # past value
      ;;
    --pbeam)
      shift # past argument
      PBEAM=$1
      shift # past value
      ;;
    --ebeam)
      shift # past argument
      EBEAM=$1
      shift # past value
      ;;
    -s|--sim)
      DO_SIM=1
      DO_ALL=
      shift # past value
      ;;
    --data-init)
      DATA_INIT=1
      DO_ALL=
      shift # past value
      ;;
    -r|--rec)
      DO_REC=1
      DO_ALL=
      shift # past value
      ;;
    -a|--analysis)
      DO_ANALYSIS=1
      DO_ALL=
      shift # past value
      ;;
    *)    # unknown option
      #POSITIONAL+=("$1") # save it in an array for later
      echo "unknown option $1"
      print_the_help
      shift # past argument
      ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# assuming something like .local/bin/env.sh has already been sourced.
print_env.sh

export BENCHMARK_TAG="tcs"
export BEAM_TAG="${EBEAM}x${PBEAM}"

FILE_NAME_TAG="${BENCHMARK_TAG}_${BEAM_TAG}_${JUGGLER_N_EVENTS}"
DATA_URL="S3/eictest/ATHENA/EVGEN/EXCLUSIVE/TCS_ABCONV/${EBEAM}x${PBEAM}/hel_minus/TCS_gen_ab_hiAcc_${EBEAM}x${PBEAM}m_${TAG}_novtx.hepmc.gz"

export TMP_PATH="${LOCAL_DATA_PATH}/tmp/${BENCHMARK_TAG}/${BEAM_TAG}"
mkdir -p "${TMP_PATH}"

export INPUT_PATH="data/${BENCHMARK_TAG}/${BEAM_TAG}"
mkdir -p "${INPUT_PATH}"

export SIM_OUTPUT_PATH="sim_output/${BENCHMARK_TAG}/${BEAM_TAG}"
mkdir -p "${SIM_OUTPUT_PATH}"

export REC_OUTPUT_PATH="results/${BENCHMARK_TAG}/${BEAM_TAG}"
mkdir -p "${REC_OUTPUT_PATH}"

export RESULTS_PATH="results/${BENCHMARK_TAG}/${BEAM_TAG}"
mkdir -p ${RESULTS_PATH}

export JUGGLER_MC_FILE="${INPUT_PATH}/mc_${FILE_NAME_TAG}.hepmc"
export JUGGLER_SIM_FILE="${SIM_OUTPUT_PATH}/sim_${FILE_NAME_TAG}.edm4hep.root"
export JUGGLER_REC_FILE="${REC_OUTPUT_PATH}/rec_${FILE_NAME_TAG}.root"

echo "FILE_NAME_TAG       = ${FILE_NAME_TAG}"
echo "JUGGLER_N_EVENTS    = ${JUGGLER_N_EVENTS}"
echo "JUGGLER_DETECTOR    = ${JUGGLER_DETECTOR}"


## To run the reconstruction, we need the following global variables:
## - JUGGLER_INSTALL_PREFIX: Install prefix for Juggler (simu/recon)
## - JUGGLER_DETECTOR:       the detector package we want to use for this benchmark
## - DETECTOR_PATH:          full path to the detector definitions

## Step 1. Get the data
if [[ -n "${DATA_INIT}" || -n "${DO_ALL}" ]] ; then
  mc -C . config host add S3 https://dtn01.sdcc.bnl.gov:9000 $S3_ACCESS_KEY $S3_SECRET_KEY
  mc -C . cat --insecure ${DATA_URL} | gunzip -c | head -n $((20+10*JUGGLER_N_EVENTS)) |  sanitize_hepmc3 > "${JUGGLER_MC_FILE}"
  if [[ "$?" -ne "0" ]] ; then
    echo "Failed to download hepmc file"
    exit 1
  else
    echo "Downloaded ${JUGGLER_MC_FILE}"
    ls -lrtha
  fi
fi

### Step 2. Run the simulation (geant4)
if [[ -n "${DO_SIM}" || -n "${DO_ALL}" ]] ; then
  ## run geant4 simulations
  ls -lrtha
  ddsim --runType batch \
    --part.minimalKineticEnergy 1000*GeV  \
    --filter.tracker edep0 \
    -v ERROR \
    --numberOfEvents ${JUGGLER_N_EVENTS} \
    --compactFile ${DETECTOR_PATH}/${JUGGLER_DETECTOR_CONFIG}.xml \
    --inputFiles "${JUGGLER_MC_FILE}" \
    --outputFile "${JUGGLER_SIM_FILE}"
  if [ "$?" -ne "0" ] ; then
    echo "ERROR running ddsim"
    exit 1
  else
    echo "Simulated ${JUGGLER_SIM_FILE}"
    ls -lrtha
  fi
fi

### Step 3. Run the reconstruction (juggler)
export PBEAM
if [[ -n "${DO_REC}" || -n "${DO_ALL}" ]] ; then
  for rec in options/*.py ; do
    unset tag
    [[ $(basename ${rec} .py) =~ (.*)\.(.*) ]] && tag=".${BASH_REMATCH[2]}"
    JUGGLER_REC_FILE=${JUGGLER_REC_FILE/.root/${tag:-}.root} \
      gaudirun.py ${rec} || [ $? -eq 4 ]
    rootls -t "${JUGGLER_REC_FILE}"
  done

  root_filesize=$(stat --format=%s "${JUGGLER_REC_FILE}")
  if [[ "${JUGGLER_N_EVENTS}" -lt "500" ]] ; then 
    # file must be less than 10 MB to upload
    if [[ "${root_filesize}" -lt "10000000" ]] ; then 
      cp "${JUGGLER_REC_FILE}" results/.
    fi
  fi
fi

### Step 4. Run the analysis code
if [[ -n "${DO_ANALYSIS}" || -n "${DO_ALL}" ]] ; then
  echo "Running analysis scripts"
  rootls -t "${JUGGLER_REC_FILE}"

  # Store all plots here (preferribly png and pdf files)
  mkdir -p results/tcs

  # here you can add as many scripts as you want.
  root -b -q "benchmarks/tcs/analysis/tcs_tests.cxx+(\"${JUGGLER_REC_FILE}\")"
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running root script"
    exit 1
  fi
fi



