#!/bin/bash

source $(dirname $0)/common.sh $*

# Reconstruct
JANA_HOME=/usr/local/lib/EICrecon eicrecon -Ppodio:output_file=${JUGGLER_REC_FILE} ${JUGGLER_SIM_FILE}
if [ "$?" -ne "0" ] ; then
  echo "ERROR running eicrecon"
  exit 1
fi
