# workflows
A collection of my Snakemake workflows for meetagenomics including assembly, annotation, and anvio preprocessing. 

In this repositories are makefiles to:

 * QC illumina reads (with bbmap and bfc), assemble with SPAdes, and idenitfy and annotate genes
 * QC illumina reads (with trimmamatic and pandaseq), compare to public databases, and tabulate counts by taxon and gene family.
 * Assemble multiple samples with megahit, map reads with BWA, and import everything into Anvi'o for visualization

## A work in progress

This repository is in the process of being refactored, documented, and tested. Please have bear with me, but don't hesitate to contact me if you have requests.

## anvio.metagenomic.snake

The multiple sample assembly and visualization workflow is the only one remotely ready for public consumption at this point. Instructions for running can be found in the initial comments.

## installation

Just clone this repository and install snakemake. I recommend using (ana)conda for this. You will also need to install software used by the workflow. The workflow comments should list the requirements.


