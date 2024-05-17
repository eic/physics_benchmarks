#!/bin/bash
source strict-mode.sh

source benchmarks/u_rho/setup.config $*

if [ -f ${JUGGLER_IN_FILE} ]; then
  echo "Input simulation file does ${JUGGLER_IN_FILE} not exist."
fi

# Simulate
/usr/bin/time -v \
ddsim --runType run \
      --printLevel WARNING \
      --numberOfEvents ${JUGGLER_N_EVENTS} \
      --part.minimalKineticEnergy 1*TeV  \
      --filter.tracker edep0 \
      --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml \
      --inputFile ${JUGGLER_IN_FILE} \
      --outputFile  ${JUGGLER_OUT_FILE}
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running ddsim"
  exit 1
fi
