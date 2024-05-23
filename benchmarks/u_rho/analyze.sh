#!/bin/bash
source strict-mode.sh


source benchmarks/u_rho/setup.config $*

OUTPUT_PLOTS_DIR=sim_output/nocampaign
mkdir -p ${OUTPUT_PLOTS_DIR}
# Analyze
/usr/bin/time -v \
root -l -b -q benchmarks/u_rho/analysis/uchannelrho.cxx+("${REC_FILE}","${OUTPUT_PLOTS_DIR}/plots.root")
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR analysis failed"
  exit 1
fi

if [ ! -d "${OUTPUT_PLOTS_DIR}_figures" ]; then
    mkdir "${OUTPUT_PLOTS_DIR}_figures"
    echo "${OUTPUT_PLOTS_DIR}_figures directory created successfully."
else
    echo "${OUTPUT_PLOTS_DIR}_figures directory already exists."
fi
root -l -b -q benchmarks/u_rho/macros/plot_rho_physics_benchmark.C("${OUTPUT_PLOTS_DIR}/plots.root")
cat sim_output/*.json
