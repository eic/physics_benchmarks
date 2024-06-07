#!/bin/bash
source strict-mode.sh

source benchmarks/u_rho/setup.config $*

if [ -f ${INPUT_FILE} ]; then
  echo "ERROR: Input simulation file does ${INPUT_FILE} not exist."
else
  echo "GOOD: Input simulation file ${INPUT_FILE} exists!"
fi

# Simulate
ddsim --runType batch \
      -v WARNING \
      --numberOfEvents ${N_EVENTS} \
      --part.minimalKineticEnergy 100*GeV  \
      --filter.tracker edep0 \
      --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml \
      --inputFiles ${INPUT_FILE} \
      --outputFile  ${OUTPUT_FILE}
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR running ddsim"
  exit 1
fi
