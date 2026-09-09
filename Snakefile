configfile: "snakemake.yml"

import functools
import os
import subprocess
from snakemake.logging import logger


## The detector install prefix is the single source of truth: the geometry XML
## directory, the geometry library and thisepic.sh are all derived from it.
## DETECTOR_PATH is published back into the configuration so that rules can
## refer to it as config["DETECTOR_PATH"] rather than reaching for a
## module-level name, and so that a detector with a non-standard layout can
## still be pointed at directly with --config DETECTOR_PATH=...
DETECTOR_PREFIX = config["DETECTOR_PREFIX"]
config.setdefault("DETECTOR_PATH", f"{DETECTOR_PREFIX}/share/{config['DETECTOR']}")

## ROOT mirrors the absolute source path underneath its build directory, so
## ROOT_BUILD_DIR has to be absolute for the declared outputs to match where
## ACLiC actually writes. An empty value means build in-tree.
ROOT_BUILD_DIR = config.get("ROOT_BUILD_DIR") or None

if ROOT_BUILD_DIR is not None:
    ROOT_BUILD_DIR = os.path.abspath(ROOT_BUILD_DIR)
    ROOT_BUILD_DIR_PREFIX = f"{ROOT_BUILD_DIR}/{os.getcwd().lstrip('/')}/"
else:
    ROOT_BUILD_DIR_PREFIX = ""


## Every shell: block gets the detector environment and ROOT settings from the
## configuration, so the pipeline does not depend on the CI job having set them
## up. Five things to keep in mind when editing this:
##  - No stray curly braces: shell.prefix() runs its argument through
##    snakemake's own format(), which would interpret them as format fields.
##    That rules out ${VAR:-} style defaults, hence the ordering below.
##  - "set -euo pipefail" comes last, not first. Snakemake normally prepends it
##    to every shell: block, but only while no shell.prefix() is set, so we have
##    to reinstate it ourselves or silently lose it. It has to come *after* the
##    environment setup: thisepic.sh reads $LD_LIBRARY_PATH unguarded and the
##    ROOT_INCLUDE_PATH append reads its own previous value, both of which abort
##    under "set -u" when unset. The rule body still runs fully strict.
##  - thisepic.sh takes the detector configuration as $1. That relies on
##    "source" rather than POSIX ".", which does not give the sourced script
##    its own positional parameters; snakemake runs shell: blocks under bash
##    (it calls shell.executable("bash") on import), so this is fine.
##  - DETECTOR_PATH is re-exported afterwards, because thisepic.sh sets it from
##    its own install tree, which would otherwise disagree with
##    config["DETECTOR_PATH"] for the rules that read $DETECTOR_PATH in a
##    shell: block.
##  - ";" and not "&&", so that the environment setup does not abort a rule
##    where thisepic.sh is absent (e.g. the lager image used by dvmp:generate).
shell.prefix(
    f"source {DETECTOR_PREFIX}/bin/thisepic.sh {config['DETECTOR_CONFIG']}; "
    f"export DETECTOR_PATH={config['DETECTOR_PATH']}; "
    f"export ROOT_MAX_THREADS={config['BENCHMARK_N_THREADS']}; "
    f"export ROOT_INCLUDE_PATH={os.path.abspath(workflow.basedir)}/include:$ROOT_INCLUDE_PATH; "
    + (f"export ROOT_BUILD_DIR={ROOT_BUILD_DIR}; " if ROOT_BUILD_DIR else "")
    + "set -euo pipefail; "
)


@functools.cache
def get_spack_package_hash(package_name):
    import json
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
        libs.append(f"{DETECTOR_PREFIX}/lib/{lib}")
    return libs


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
