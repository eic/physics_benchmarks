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
    shell: """
xrdcp root://dtn-eic.jlab.org//work/eic2/{output.filepath} {output.filepath}
"""


rule warmup_run:
    output:
        "warmup/{DETECTOR_CONFIG}.edm4hep.root",
    message: "Ensuring that calibrations/fieldmaps are available for {wildcards.DETECTOR_CONFIG}"
    shell: """
ddsim \
  --runType batch \
  --numberOfEvents 1 \
  --compactFile "$DETECTOR_PATH/{wildcards.DETECTOR_CONFIG}.xml" \
  --outputFile "{output}" \
  --enableGun
"""

include: "benchmarks/diffractive_vm/Snakefile"
include: "benchmarks/dis/Snakefile"
include: "benchmarks/demp/Snakefile"
include: "benchmarks/zdc_lyso/Snakefile"
