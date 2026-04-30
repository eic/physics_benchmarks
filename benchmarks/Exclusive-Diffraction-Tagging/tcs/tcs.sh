#!/bin/bash
source strict-mode.sh
source benchmarks/Exclusive-Diffraction-Tagging/options.sh

function print_the_help {
  echo "USAGE: ${0} --ebeam E --pbeam P --tag T [--sim] [--rec] [--analysis] [--all]"
  echo "SCRIPT OPTIONS:"
  echo "  --ebeam E       Electron beam energy"
  echo "  --pbeam P       Ion beam energy"
  echo "  --tag T         Event file tag"
  print_step_options_help
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

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -h|--help)
      print_the_help
      ;;
    --all)
      DO_ALL=1
      if [[ -n "${DO_REC}${DO_SIM}${DO_ANALYSIS}" ]]; then
        echo "Error: cannot use --all with other step arguments." 1>&2
        print_the_help
        exit 1
      fi
      shift
      ;;
    --tag)
      TAG="$2"
      shift; shift
      ;;
    --pbeam)
      PBEAM="$2"
      shift; shift
      ;;
    --ebeam)
      EBEAM="$2"
      shift; shift
      ;;
    -s|--sim)
      DO_SIM=1
      DO_ALL=
      shift
      ;;
    --data-init)
      DATA_INIT=1
      DO_ALL=
      shift
      ;;
    -r|--rec)
      DO_REC=1
      DO_ALL=
      shift
      ;;
    -a|--analysis)
      DO_ANALYSIS=1
      DO_ALL=
      shift
      ;;
    *)
      echo "unknown option '$1'" 1>&2
      print_the_help
      exit 1
      ;;
  esac
done

# assuming something like .local/bin/env.sh has already been sourced.
print_env.sh

FILE_NAME_TAG="tcs_${EBEAM}x${PBEAM}m_${TAG}"
XROOTD_BASEURL="root://dtn-eic.jlab.org//volatile/eic/EPIC"
INPUT_FILE="EVGEN/EXCLUSIVE/TCS_ABCONV/${EBEAM}x${PBEAM}/hel_minus/TCS_gen_ab_hiAcc_${EBEAM}x${PBEAM}m_${TAG}.hepmc3.tree.root"

export JUGGLER_MC_FILE="${XROOTD_BASEURL}/${INPUT_FILE}"
export JUGGLER_SIM_FILE="${LOCAL_DATA_PATH}/sim_${FILE_NAME_TAG}.edm4hep.root"
export JUGGLER_REC_FILE="${LOCAL_DATA_PATH}/rec_${FILE_NAME_TAG}.root"

echo "FILE_NAME_TAG       = ${FILE_NAME_TAG}"
echo "JUGGLER_N_EVENTS    = ${JUGGLER_N_EVENTS}"
echo "DETECTOR    = ${DETECTOR}"


## To run the reconstruction, we need the following global variables:
## - DETECTOR:       the detector package we want to use for this benchmark
## - DETECTOR_PATH:          full path to the detector definitions

### Step 1. Run the simulation (geant4)
if [[ -n "${DO_SIM}" || -n "${DO_ALL}" ]] ; then
  ## run geant4 simulations
  npsim --runType batch \
    --part.minimalKineticEnergy 1000*GeV  \
    --filter.tracker edep0 \
    -v ERROR \
    --numberOfEvents ${JUGGLER_N_EVENTS} \
    --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml \
    --inputFiles "${JUGGLER_MC_FILE}" \
    --outputFile  ${JUGGLER_SIM_FILE}
  if [ "$?" -ne "0" ] ; then
    echo "ERROR running npsim"
    exit 1
  fi
fi

### Step 2. Run the reconstruction (eicrecon)
export PBEAM
if [[ -n "${DO_REC}" || -n "${DO_ALL}" ]] ; then
  eicrecon ${JUGGLER_SIM_FILE} -Ppodio:output_file=${JUGGLER_REC_FILE}
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running eicrecon"
    exit 1
  fi

  root_filesize=$(stat --format=%s "${JUGGLER_REC_FILE}")
  if [[ "${JUGGLER_N_EVENTS}" -lt "500" ]] ; then 
    # file must be less than 10 MB to upload
    if [[ "${root_filesize}" -lt "10000000" ]] ; then 
      cp ${JUGGLER_REC_FILE} results/.
    fi
  fi
fi

### Step 3. Run the analysis code
if [[ -n "${DO_ANALYSIS}" || -n "${DO_ALL}" ]] ; then
  echo "Running analysis scripts"
  rootls -t  ${JUGGLER_REC_FILE}

  # Store all plots here (preferribly png and pdf files)
  mkdir -p results/tcs

  # here you can add as many scripts as you want.
  root -b -q "benchmarks/Exclusive-Diffraction-Tagging/tcs/analysis/tcs_tests.cxx+(\"${JUGGLER_REC_FILE}\")"
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running root script"
    exit 1
  fi
fi



