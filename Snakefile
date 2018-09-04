import glob
import os


SUBSETS = ['fungi','plant','invertebrate','vertebrate_mammalian', 'vertebrate_other']
#SUBSETS = 'fungi' 
INPUTS = {}
NAMES ={}

DATADIR = "data"
OUTDIR = "euk_sbts"

#for subset in SUBSETS:
#    INPUTS[subset] = glob.glob(os.path.join(DATADIR, subset,"*/*.fna.gz"))
#    NAMES[subset] = [os.path.basename(f).split('_genomic.fna.gz')[0] for f in INPUTS[subset] if '_rna_from_genomic' not in f]

rule all:
    input:
       # expand(os.path.join(OUTDIR,'trees/{subset}/{subset}-k{ksize}.sbt.json'),
       #        subset=SUBSETS,
       #        ksize=[21, 31, 51])
       dynamic(expand("data/genbank/{subset}/{{id}}/{{id}}_{{name}}_rm.out.gz", subset=SUBSETS))


rule download_ncbi_genomes:
    params: sub=SUBSETS, out='data'
    output: dynamic("data/genbank/{subset}/{id}/{id}_{name}_rm.out.gz")
    conda: "envs/env.yml"
    shell: """
           #ncbi-genome-download --section genbank {params}
		   ncbi-genome-download -o {params.out} -s genbank -v -F rm fungi
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
