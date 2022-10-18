#!/bin/bash

source $(dirname $0)/common.sh $*

# Reconstruct
JANA_HOME=/usr/local/lib/EICrecon run_eicrecon_reco_flags.py ${JUGGLER_SIM_FILE} $(basename ${JUGGLER_REC_FILE} .root)

if [ "$?" -ne "0" ] ; then
  echo "ERROR running eicrecon"
  exit 1
fi
