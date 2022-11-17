#!/bin/bash

source $(dirname $0)/common.sh $*

# Reconstruct
for rec in options/*.py ; do
  unset tag
  [[ $(basename ${rec} .py) =~ (.*)\.(.*) ]] && tag=".${BASH_REMATCH[2]}"
  JUGGLER_REC_FILE=${JUGGLER_REC_FILE/.root/${tag:-}.root} \
    /usr/bin/time -v \
    gaudirun.py ${JUGGLER_GAUDI_OPTIONS:-} ${rec}
  if [[ "$?" -ne "0" ]] ; then
    echo "ERROR running juggler"
    exit 1
  fi
done

/usr/bin/time -v run_eicrecon_reco_flags.py ${JUGGLER_SIM_FILE} ${JUGGLER_REC_FILE_BASE}
if [ "$?" -ne "0" ] ; then
  echo "ERROR running eicrecon"
  exit 1
fi

rootls -t ${JUGGLER_REC_FILE_BASE}.tree.edm4eic.root
