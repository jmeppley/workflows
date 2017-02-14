# workflows
A collection of Snakemake workflows for metagenomics

## Prerequisites

For the most part, all you need is (ana)conda. The only tools used that aren't available in conda are pandaseq (though, you can use pear or flash) and megahit (for OSX).

## Installation

Just clone this repository and set up you conda environment. Each of the primary workflows (the ones in the root folder) have an associated conda
configuration file in tests/conda. 

You can create the necessary conda environments with thses files. For example, to run the clean.illumina.snake workflow, create a conda environment with:

```
conda env create -n qc -f test/conda/illumina.qc.yml
```

## tests

Currently only bats tests are available. Install bats and run 

```
bats test/bats
```

The tests will create the necessary conda environments and run some example workflows

# Workflows

## Illumina Prep

clean.illumina.snake will take pars of raw illumina reads and create joined, cleaned reads using:

 * trimmomatic (to remove adapters)
 * bbduk (to remove more adapters)
 * pear, pandaseq, or flash to assemble pairs
 * custom code to connected unjoined pairs with N's
 * trimmomatic to end trim low quality bases from everything

Conda env file:

 * test/conda/illumina.qc.yml

## anvio metagenomics

anvio.metagenomic.snake will take reads from multiple samples and create the files necessary to load into anvio metagenomic mode.

 * reads cleaned (optionally) with bbduk, bfc, and trimmomatic
 * reads assembled (or provide your own contigs) with megahit
 * reads mapped to contigs with bwa
 * everything imported into anvio (including sample clustering)

Conda env files

 * test/conda/assembly.yml (for running snakemake)
 * test/conda/anvi2.yml (for anvio steps)

## read annotation

annotate.reads.snake will compare reads to reference DBs to cross compile hits to taxonomic referecens and gene faimly DBs. Lastal can be used for protein sequence dbs like RefSeq or KEGG, and hmmer for profile searches against the likes of COG and PFAM.
