import glob
import os

#SUBSETS = ['fungi','plant','invertebrate','vertebrate_mammalian', 'vertebrate_other']
SUBSETS = ['fungi'] 
INPUTS = {}
NAMES ={}

DATADIR = "data"

rule all:
    input: 
        # dynamic(expand("data/genbank/{subset}/{{id}}/{{name}}_rm.out.gz", subset=SUBSETS)),
        dynamic(expand("data/sigs/{subset}/{{id}}/{{name}}.sig", subset = SUBSETS))

rule download_ncbi_genomes:
    params: 
        sub = SUBSETS,
        out = "data"
    output: dynamic("data/genbank/{subset}/{id}/{name}_rm.out.gz")
    conda: "envs/env.yml"
    shell: """
           touch {output}
           #ncbi-genome-download -o {params.out} -s genbank -v -F rm fungi
           """

rule compute_sigs:
    input: "data/genbank/{subset}/{id}/{name}_rm.out.gz"
    output: "data/sigs/{subset}/{id}/{name}.sig"
    conda: "envs/env.yml"
    shell:'''
        sourmash compute -k 21,31,51 --scaled 1000 --track-abundance --name-from-first -o {output} {input}
    '''
 
rule index_sigs:
    output: 'trees/{subset}/{subset}-k{ksize}.sbt.json'
    input: dynamic("data/sigs/{subset}/{name}.sig") 
    params:
        ksize="{ksize}",
        subset="{subset}"
    conda: "envs/env.yml"
    shell: """
        sourmash index -k {params.ksize} \
                       -x 1e6 \
                       --traverse-directory \
                       {output} \
                       `dirname {input[0]}`
    """


#rule compute_sigs:
#    input: os.path.join(DATADIR,'{subset}/{id}/{id}_{name}_genomic.fna.gz')
#    output: os.path.join(OUTDIR, 'sigs/{subset}/{id}_{name}.sig')
#    wildcard_constraints:
#        id="\w+_\d+\.\d+"
#    conda: "envs/env.yml"
#    shell: """
#        sourmash compute -k 21,31,51 \
#                         --scaled 1000 \
#                         --track-abundance \
#                         --name-from-first \
#                         -o {output} \
#                         {input}
#    """

#def sigs_for_tree(w):
#    return [os.path.join(OUTDIR,"sigs/{subset}/{name}.sig").format(subset=w.subset, name=name)
#            for name in NAMES[w.subset]] 

#rule index_sigs:
#    output: os.path.join(OUTDIR,'trees/{subset}/{subset}-k{ksize}.sbt.json')
#    input: sigs_for_tree 
#    params:
#        ksize="{ksize}",
#        subset="{subset}"
#    conda: "envs/env.yml"
#    shell: """
#        sourmash index -k {params.ksize} \
#                       -x 1e6 \
#                       --traverse-directory \
#                       {output} \
#                       `dirname {input[0]}`
#    """
