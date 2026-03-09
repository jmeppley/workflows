
# **workflows** (aloha_3_1)
This repository contains collection of Snakemake workflows for metagenomics. The aloha_3_1 branch has been pruned down to (mostly) the workflows used in building version 3.1 of the ALOHA gene catalog.

## Prerequisites

For the most part, all you need is git and conda.

## Installation

Just clone this repository and set up you conda environment. Each of the primary workflows (the ones in the root folder) have an associated conda
configuration file in tests/conda. 

You can create the necessary conda environments with thses files. For example, to run the clean.illumina.snake workflow, create a conda environment with:

```
conda env create -n qc -f test/conda/illumina.qc.yml
```

## tests

The tests are run via `bats`. Install bats and run 

```
bats test/bats
```

The tests will create the necessary conda environments and run some example workflows

## Workflows

### Assembly

`assembly.metagenomic.snake` will take raw illumina reads, clean them, and assemble them into contigs

Conda environment definition file:

 * test/conda/assembly.yml
 
Configuration:

```yaml
trimmomatic:
    threads: 9
    minlen: 100
    leading: 10
    trailing: 10
    sliding_window: "4:20"
bbduk:
    threads: 9
cmsearch:  
    threads: 9
bfc:
    threads: 9
    params: "-k 21"
assembler: megahit
megahit:
    threads: 20

assembly_name: HSDXXX
sample_data:
   HSDXXX_RUN1:
      raw:
       - /path/to/R1.fastq
       - /path/to/R2.fastq
```

Run the workflow (once in each assembly directory) with:

```bash
 $ snakemake -s /path/to/workflows/assembly.metagenomic.snake --configfile config.yaml -j {threads}
```

### Gene catalog

`annotation.gene_catalog.snake` will take genes from multiple assemblies and cluster them into one catalog of non-redundant genes.

Conda environment definition file:

 * test/conda/catalog.yml
 
Configuration:
```yaml
assembly_list_file: {file listing assembly directories}
gene_file_root: contigs.all.annotations
cross_tab: True
clade_ranks: 
 - order
 - genus
output_style: long
hmmer:
    threads: 4
lastal:
    threads: 20
prodigal:
    threads: 20
aa_conversion: prodigal
clustering_method: mmseqs2
mmseqs2:
    tmp_dir_root: {path to fast temporary filesystem}
    threads: 80
dbs:
    GTDB: 
        path: {path to last-formated GTDB genes}
        format: lastp
        type: tax
    RefSeq: 
        path: {path to last-formated RefSeq genes}
        format: lastp
        type: tax
    KEGG: 
        path: {path to last-formatted KEGG genes}
        format: lastp
        assign_type: kegg
    PFAM:
        path: {path to PFAM HMMs}
        frags: 10
    COG:
        path: {path to COG HMMs}
        frags: 10
    TIGRFAM:
        path: {path to TIGR HMMs}
        frags: 10
```

Run the workflow with:

```bash
snakemake -s /path/to/workflows/annotation.gene_catalog.snake --configfile=config.yaml -p -k -j {threads} --notemp
```

## gene abundances

`scripts/gene_abundances.snake` will calculate abundances of gene sequences by mapping reads from each assembly to the final gene catalog.

Conda environment definition file:

 * scripts/gene_abundances.yaml
 
Run the workflow with:

```bash
snakemake -s /path/to/workflows/scripts/gene_abundances.snake \
    --config gene_cluster_reps=/path/to/all_genes.ffn.mmseqs.cluster.setcover.95.ffn \
             assembly_root_dir=/path/to/assemblies \
    -p -k -j 20 --notemp \
 > abundance.log 2>&1
```
