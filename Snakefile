import functools
import os
import shlex
from snakemake.logging import logger

configfile: "snakemake.yml"

# Keep the YAML as the canonical source of benchmark defaults; the fallbacks below
# preserve local CLI use when the config file is absent or intentionally overridden.
DETECTOR_PREFIX = config.get("DETECTOR_PREFIX", "/opt/detector/epic-main")
DETECTOR_CONFIG = config.get("DETECTOR_CONFIG", "epic_craterlake")
BENCHMARK_N_EVENTS = config.get("BENCHMARK_N_EVENTS", 100)
BENCHMARK_N_THREADS = config.get("BENCHMARK_N_THREADS", 10)
BENCHMARK_RNG_SEED = config.get("BENCHMARK_RNG_SEED", 1)

DETECTOR_PREFIX_SHELL = shlex.quote(DETECTOR_PREFIX)
THIS_EPIC_SHELL = shlex.quote(os.path.join(DETECTOR_PREFIX, "bin", "thisepic.sh"))
DETECTOR_CONFIG_SHELL = shlex.quote(DETECTOR_CONFIG)
BENCHMARK_N_EVENTS_SHELL = shlex.quote(str(BENCHMARK_N_EVENTS))
BENCHMARK_N_THREADS_SHELL = shlex.quote(str(BENCHMARK_N_THREADS))
BENCHMARK_RNG_SEED_SHELL = shlex.quote(str(BENCHMARK_RNG_SEED))

shell.prefix(
    "set -e; "
    f"if [ -f {THIS_EPIC_SHELL} ]; then source {THIS_EPIC_SHELL}; fi; "
    f"export DETECTOR_PREFIX={DETECTOR_PREFIX_SHELL}; "
    f"export DETECTOR_CONFIG={DETECTOR_CONFIG_SHELL}; "
    f"export BENCHMARK_N_EVENTS={BENCHMARK_N_EVENTS_SHELL}; "
    f"export BENCHMARK_N_THREADS={BENCHMARK_N_THREADS_SHELL}; "
    f"export BENCHMARK_RNG_SEED={BENCHMARK_RNG_SEED_SHELL}; "
    f"export JUGGLER_N_EVENTS={BENCHMARK_N_EVENTS_SHELL}; "
    f"export ROOT_MAX_THREADS={BENCHMARK_N_THREADS_SHELL}; "
    f"export JUGGLER_RNG_SEED={BENCHMARK_RNG_SEED_SHELL}; "
)


@functools.cache
def get_spack_package_hash(package_name):
    import json
    import subprocess
    try:
        ver_info = json.loads(subprocess.check_output(["spack", "find", "--json", package_name]))
        return ver_info[0]["hash"]
    except FileNotFoundError as e:
        logger.warning("Spack is not installed")
        return ""
    except subprocess.CalledProcessError as e:
        print(e)
        return ""


@functools.cache
def find_epic_libraries():
    import ctypes.util
    # if library is not found (not avaliable) we return an empty list to let DAG still evaluate
    libs = []
    lib = ctypes.util.find_library("epic")
    if lib is not None:
        libs.append(os.environ["DETECTOR_PATH"] + "/../../lib/" + lib)
    return libs


ROOT_BUILD_DIR = os.getenv("ROOT_BUILD_DIR", None)

if ROOT_BUILD_DIR is not None:
    ROOT_BUILD_DIR_PREFIX = f"{ROOT_BUILD_DIR.rstrip('/')}/{os.getcwd().lstrip('/')}/"
else:
    ROOT_BUILD_DIR_PREFIX = ""


rule compile_analysis:
    input:
        "{path}/{filename}.cxx",
    output:
        ROOT_BUILD_DIR_PREFIX + "{path}/{filename}_cxx.d",
        ROOT_BUILD_DIR_PREFIX + "{path}/{filename}_cxx.so",
        ROOT_BUILD_DIR_PREFIX + "{path}/{filename}_cxx_ACLiC_dict_rdict.pcm",
    shell:
        """
root -l -b -q -e '.L {input}+'
"""


rule fetch_epic:
    output:
        filepath="EPIC/{PATH}"
    cache: True
    shell: """
xrdcp root://dtn-eic.jlab.org//volatile/eic/{output.filepath} {output.filepath}
"""


rule warmup_run:
    output:
        "warmup/{DETECTOR_CONFIG}.edm4hep.rnt.root",
    message: "Ensuring that calibrations/fieldmaps are available for {wildcards.DETECTOR_CONFIG}"
    shell: """
ddsim \
  --runType batch \
  --numberOfEvents 1 \
  --compactFile "$DETECTOR_PATH/{wildcards.DETECTOR_CONFIG}.xml" \
  --outputConfig.useRNTuple true \
  --outputFile "{output}" \
  --enableGun
"""

include: "benchmarks/Exclusive-Diffraction-Tagging/demp/Snakefile"
include: "benchmarks/Exclusive-Diffraction-Tagging/diffractive_vm/Snakefile"
include: "benchmarks/Exclusive-Diffraction-Tagging/dvmp/Snakefile"
include: "benchmarks/Exclusive-Diffraction-Tagging/semi_coherent/Snakefile"
include: "benchmarks/Jets-HF/jets/Snakefile"
include: "benchmarks/Inclusive/dis/Snakefile"
