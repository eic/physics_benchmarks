#!/bin/bash

source $(dirname $0)/common.sh $*

# Analyze
/usr/bin/time -v \
root -l -b -q "benchmarks/single/analysis/analyze.cxx+(\"${JUGGLER_REC_FILE}\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR analysis failed"
  exit 1
fi

python benchmarks/dis/analysis/truth_reconstruction.py --rec_file ${JUGGLER_REC_FILE} --config ${JUGGLER_FILE_NAME_TAG}_${DETECTOR_CONFIG} --results_path ${RESULTS_PATH} --nevents ${JUGGLER_N_EVENTS}
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running truth_reconstruction script"
  exit 1
fi
