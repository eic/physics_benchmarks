#!/bin/bash

## =============================================================================
## Only read in the full simulations from s3 and do analysis
## =============================================================================

## make sure we launch this script from the project root directory
PROJECT_ROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/../..
pushd ${PROJECT_ROOT}

source .local/bin/env.sh

FILE_NAME_TAG="dvcs-d"
REC_FILE="${FILE_NAME_TAG}_input"
JUGGLER_DETECTOR="ATHENA"
RESULTS_PATH="results/s3_full/${FILE_NAME_TAG}"
BEAM_TAG="10X100"

mkdir -p "s3_full/${FILE_NAME_TAG}"
mkdir -p "results/s3_full/${FILE_NAME_TAG}"

INPUT_PATH_FROM_S3_TAG="s3_full/${FILE_NAME_TAG}"

echo "Local print out: "
echo "-----------------"
echo "This is running analysis from root files on S3. Many local variables are set"
echo "Input will be copy to ${INPUT_PATH_FROM_S3_TAG}"
echo "Output will be stored at ${RESULTS_PATH}" 

echo "Running the dvcs-d analysis"

if [ -f "${INPUT_PATH_FROM_S3_TAG}/${REC_FILE}" ]; then
  echo "Found file and already downloaded."
else
  mc -C . config host add S3 https://dtn01.sdcc.bnl.gov:9000 $S3_ACCESS_KEY $S3_SECRET_KEY
  for i in {10..20}
  do

    DATA_URL="S3/eictest/ATHENA/RECO/acadia-v2.1/EXCLUSIVE/DVCS_ABCONV/10x100/DVCS.1.ab.hiAcc.10x100_novtx.00${i}.root"
    mc -C . cp --insecure ${DATA_URL} "${INPUT_PATH_FROM_S3_TAG}/${REC_FILE}_${i}.root"
  done
  hadd -f "${REC_FILE}.root" "${INPUT_PATH_FROM_S3_TAG}/${REC_FILE}_"*".root"
fi


## =============================================================================
## Step 4: Analysis
## write a temporary configuration file for the analysis script
echo "Running analysis"
CONFIG="${INPUT_PATH_FROM_S3_TAG}/${FILE_NAME_TAG}.json"
cat << EOF > ${CONFIG}
{
  "rec_file": "${INPUT_PATH_FROM_S3_TAG}/${REC_FILE}.root",
  "detector": "${JUGGLER_DETECTOR}",
  "output_prefix": "${RESULTS_PATH}/${FILE_NAME_TAG}",
  "test_tag": "${BEAM_TAG}"
}
EOF
#cat ${CONFIG}
root -b -q "benchmarks/dvcs-d/analysis/dvcs_d_analysis.cxx+(\"${CONFIG}\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running rec_dvcs-d_analysis script"
  exit 1
fi

echo "DVCS-d benchmarks complete"
