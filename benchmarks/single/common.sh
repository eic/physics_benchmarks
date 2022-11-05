if [[ ! -n  "${JUGGLER_N_EVENTS}" ]] ; then 
  export JUGGLER_N_EVENTS=100
fi

export JUGGLER_FILE_NAME_TAG="${1:-e-_1GeV_45to135deg}"
export JUGGLER_GEN_FILE="benchmarks/single/${JUGGLER_FILE_NAME_TAG}.steer"
export JUGGLER_SIM_FILE="sim_output/sim_${JUGGLER_FILE_NAME_TAG}.edm4hep.root"
export JUGGLER_REC_FILE="sim_output/rec_${JUGGLER_FILE_NAME_TAG}.root"

export BENCHMARK_TAG="single"
echo "Setting up the local environment for the ${BENCHMARK_TAG^^} benchmarks"

RESULTS_PATH="results/${BENCHMARK_TAG}"
mkdir -p ${RESULTS_PATH}
mkdir -p "${RESULTS_PATH}/truth_reconstruction"
export RESULTS_PATH=`realpath ${RESULTS_PATH}`
echo "RESULTS_PATH:           ${RESULTS_PATH}"
