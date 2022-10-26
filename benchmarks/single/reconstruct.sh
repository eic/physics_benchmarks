#!/bin/bash

source $(dirname $0)/common.sh $*

# Reconstruct
run_eicrecon_reco_flags.py ${JUGGLER_SIM_FILE} ${JUGGLER_REC_FILE_BASE}
if [ "$?" -ne "0" ] ; then
  echo "ERROR running eicrecon"
  exit 1
fi

rootls -t ${JUGGLER_REC_FILE_BASE}.tree.edm4eic.root
