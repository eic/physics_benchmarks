rule compile_analysis:
    input:
        "{path}/{filename}.cxx",
    output:
        "{path}/{filename}_cxx.d",
        "{path}/{filename}_cxx.so",
        "{path}/{filename}_cxx_ACLiC_dict_rdict.pcm",
    shell:
        """
root -l -b -q -e '.L {input}+'
"""

include: "benchmarks/diffractive_vm/Snakefile"
