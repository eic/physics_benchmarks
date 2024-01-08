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

include: "benchmarks/diffractive_vm/Snakefile"
