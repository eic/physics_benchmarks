#!/bin/bash
source strict-mode.sh
source benchmarks/your_benchmark/setup.config $*

# Reconstruct
if [ ${RECO} == "eicrecon" ] ; then
  eicrecon ${OUTPUT_FILE} -Ppodio:output_file=${REC_FILE}
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

if [ -f jana.dot ] ; then cp jana.dot ${REC_FILE_BASE}.dot ; fi

#rootls -t ${REC_FILE_BASE}.tree.edm4eic.root
rootls -t ${REC_FILE}
