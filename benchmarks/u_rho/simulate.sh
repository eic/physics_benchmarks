#!/bin/bash
source strict-mode.sh

source benchmarks/u_rho/setup.config $*

JUGGLER_IN_FILE=$1
JUGGLER_OUT_FILE=$2
JUGGLER_N_EVENTS=$3

if [ -f ${JUGGLER_IN_FILE} ]; then
  echo "ERROR: Input simulation file does ${JUGGLER_IN_FILE} not exist."
else
  echo "GOOD: Input simulation file ${JUGGLER_IN_FILE} exists!"
fi

# Simulate
ddsim --runType batch \
      -v WARNING \
      --numberOfEvents ${JUGGLER_N_EVENTS} \
      --part.minimalKineticEnergy 100*GeV  \
      --filter.tracker edep0 \
      --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml \
      --inputFiles ${JUGGLER_IN_FILE} \
      --outputFile  ${JUGGLER_OUT_FILE}
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running ddsim"
  exit 1
fi
