import glob
import os
import pandas as pd


SUBSETS = ['fungi'] # 'archaea', 'protozoa', 'plant','invertebrate','vertebrate_mammalian', 'vertebrate_other', 'bacteria', 'viral']
INPUTS = {}
NAMES ={}

DATADIR = "genbank"
OUTDIR = ""
METADATA = {}

for subset in SUBSETS:
    INPUTS[subset] = glob.glob(os.path.join(DATADIR, subset,"*/*_rna_from_genomic.fna.gz"))
    NAMES[subset] = [os.path.basename(f).split('_rna_from_genomic.fna.gz')[0] for f in INPUTS[subset]] #if '_rna_from_genomic' not in f]

for subset in SUBSETS:
    METADATA[subset] = (pd.read_table('{}-metadata.tab'.format(subset))
                          .set_index('assembly_accession'))

rule all:
    input:
        expand(os.path.join(OUTDIR,'trees/{subset}/{subset}-k{ksize}.sbt.json'),
               subset=SUBSETS,
               ksize=[21, 31, 51])

def generate_sig_name(w):
    metadata = METADATA[subset].loc[w.id]
    return "{} {}".format(metadata[w.id], metadata['organism_name']) 

rule compute_sigs:
    input: os.path.join(DATADIR,'{subset}/{id}/{id}_{name}_rna_from_genomic.fna.gz')
    output: os.path.join(OUTDIR, 'sigs/{subset}/{id}_{name}.sig')
    params:
        id="{id}",
        name=generate_sig_name
    wildcard_constraints:
        id="\w+_\d+\.\d+"
    conda: "envs/env.yml"
    shell: """
        sourmash compute -k 21,31,51 \
                         --scaled 1000 \
                         --track-abundance \
                         --merge {params.name:q} \
                         -o {output} \
                         {input}
    """

def sigs_for_tree(w):
    return [os.path.join(OUTDIR,"sigs/{subset}/{name}.sig").format(subset=w.subset, name=name)
            for name in NAMES[w.subset]] 

rule index_sigs:
    output: os.path.join(OUTDIR,'trees/{subset}/{subset}-k{ksize}.sbt.json')
    input: sigs_for_tree 
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
