#!/bin/bash
source strict-mode.sh

source $(dirname $0)/common.sh $*

# Reconstruct
if [ ${RECO} == "eicrecon" ] ; then
  eicrecon ${JUGGLER_SIM_FILE} -Ppodio:output_file=${JUGGLER_REC_FILE}
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running eicrecon"
    exit 1
  fi
fi

if [[ ${RECO} == "juggler" ]] ; then
  gaudirun.py options/reconstruction.py || [ $? -eq 4 ]
  if [ "$?" -ne "0" ] ; then
    echo "ERROR running juggler"
    exit 1
  fi
fi

if [ -f jana.dot ] ; then cp jana.dot ${JUGGLER_REC_FILE_BASE}.dot ; fi

rootls -t ${JUGGLER_REC_FILE_BASE}.tree.edm4eic.root
