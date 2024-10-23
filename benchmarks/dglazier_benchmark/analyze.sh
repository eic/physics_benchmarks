#!/bin/bash
source strict-mode.sh
source benchmarks/dglazier_benchmark/setup.config $*

OUTPUT_PLOTS_DIR=sim_output/nocampaign
mkdir -p ${OUTPUT_PLOTS_DIR}
# Analyze
command time -v \
root -l -b -q "benchmarks/dglazier_benchmark/analysis/uchannelrho.cxx(\"${REC_FILE}\",\"${OUTPUT_PLOTS_DIR}/plots.root\")"
if [[ "$?" -ne "0" ]] ; then
  echo "ERROR analysis failed"
  exit 1
fi

if [ ! -d "${OUTPUT_PLOTS_DIR}/plots_figures" ]; then
    mkdir "${OUTPUT_PLOTS_DIR}/plots_figures"
    echo "${OUTPUT_PLOTS_DIR}/plots_figures directory created successfully."
else
    echo "${OUTPUT_PLOTS_DIR}/plots_figures directory already exists."
fi
root -l -b -q "benchmarks/dglazier_benchmark/macros/plot_rho_physics_benchmark.C(\"${OUTPUT_PLOTS_DIR}/plots.root\")"
cat benchmark_output/*.json
