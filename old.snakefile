import glob
import os
import re


SUBSETS = ['archaea', 'fungi', 'protozoa', 'plant','invertebrate','vertebrate_mammalian', 'vertebrate_other'] #, 'bacteria', 'viral']
INPUTS = {}
NAMES ={}

DATADIR = "genbank"
OUTDIR = ""

for subset in SUBSETS:
    INPUTS[subset] = glob.glob(os.path.join(DATADIR, subset,"*/*_rna_from_genomic.fna.gz"))
    NAMES[subset] = [os.path.basename(f).split('_rna_from_genomic.fna.gz')[0] for f in INPUTS[subset]] #if '_rna_from_genomic' not in f]

rule all:
    input:
        expand(os.path.join(OUTDIR,'trees/{subset}/{subset}-k{ksize}.sbt.json'),
               subset=SUBSETS,
               ksize=[21, 31, 51])

def generate_sig_name(w):
    file_name = os.path.join(DATADIR, f'{w.subset}/{w.id}/{w.id}_{w.name}_assembly_report.txt')
    #import pdb;pdb.set_trace()
    if os.path.isfile(file_name):
        with open(file_name, 'r') as f:
            for i, line in enumerate(f):
                if i == 0:
                    accession = line.strip()
                    accession = re.sub("# Assembly name: ", "", accession)
                if i == 1:
                    organism = line.strip()
                    organism = re.sub("# Organism name: ", "", organism)
        return "{} {}".format(accession, organism)
    else:
        path = os.path.join(DATADIR, f'{w.subset}/{w.id}/*assembly_report.txt')
        for file_name in glob.glob(path):
            with open(file_name, 'r') as f:
                for i, line in enumerate(f):
                    if i == 0:
                        accession = line.strip()
                        accession = re.sub("# Assembly name: ", "", accession)
                    if i == 1:
                        organism = line.strip()
                        organism = re.sub("# Organism name: ", "", organism)
            return "{} {}".format(accession, organism)


rule compute_sigs:
    input: os.path.join(DATADIR,'{subset}/{id}/{id}_{name}_rna_from_genomic.fna.gz')
    output: os.path.join(OUTDIR, 'sigs/{subset}/{id}_{name}.sig')
    params:
        id="{id}",
        subset = "{subset}",
        name="{name}",
        sig_name=generate_sig_name
    wildcard_constraints:
        id="\w+_\d+\.\d+"
    conda: "envs/env.yml"
    shell: """
        sourmash compute -k 21,31,51 \
                         --scaled 1000 \
                         --track-abundance \
                         --merge {params.sig_name:q} \
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
