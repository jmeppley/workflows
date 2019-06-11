setup() {
    eval "$(conda shell.bash hook)"
    mkdir -p test/conda/envs
    ENV=illumina.qc
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi
    conda activate $ENV_DIR
}

@test "Prep a pair of scriptseq files with pear" {
    rm -rf test/scratch/qc.pear
    mkdir -p test/scratch/qc.pear
    cd test/scratch/qc.pear

    ln -s ../../data/raw_reads/2014_ALOHA_XVII_1-1B_S1_R1_001.head40k.fastq aloha1b.R1.fastq
    ln -s ../../data/raw_reads/2014_ALOHA_XVII_1-1B_S1_R2_001.head40k.fastq aloha1b.R2.fastq
    run bash -c "snakemake -j 10 -s ../../../meta.snake -p --config workflow=qc/pear.snake --verbose aloha1b.scripseq.ATCACG.trim_adapt.joined.fastq > pear.join.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Prep two pair of scriptseq files with flash" {
    rm -rf test/scratch/qc.flash
    mkdir -p test/scratch/qc.flash
    cd test/scratch/qc.flash

    run bash -c "snakemake -j 10 -s ../../../clean.illumina.snake -p --configfile ../../data/configs/illumina.yaml --verbose --config joining_program=flash discover_fastx_for_stats=True > flash.join.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Prep two pair of scriptseq files with pandaseq" {
    which pandaseq || skip "Pandaseq is not installed"

    rm -rf test/scratch/qc.pandaseq
    mkdir -p test/scratch/qc.pandaseq
    cd test/scratch/qc.pandaseq

    run bash -c "snakemake -j 10 -s ../../../clean.illumina.snake -p --configfile ../../data/configs/illumina.yaml --verbose --config joining_program=pandaseq > pandaseq.join.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Join and clean reads, and count spiked standards" {
    rm -rf test/scratch/qc.spike
    mkdir -p test/scratch/qc.spike
    cd test/scratch/qc.spike

    run bash -c "snakemake -j 10 -s ../../../normalize_by_standards.snake -p --configfile ../../data/configs/illumina.yaml --verbose --config joining_program=flash discover_fastx_for_stats=True cleaning_protocol=join_and_standards standards_fasta=../../data/other/standards.fasta spike_amounts_table=../../data/other/spike_amounts.tsv > spikes.log 2>&1"
    [ "$status" -eq 0 ]
    [ -e stats/fortyk.scripseq.ATCACG.trim_adapt.solo.paired.noadaptp.trimmed.fastq.hist ]
}

@test "Split reads into rRNA or not" {
    rm -rf test/scratch/qc.rna
    mkdir -p test/scratch/qc.rna
    cd test/scratch/qc.rna

    ln -s ../../data/raw_reads/2014_ALOHA_XVII_1-1B_S1_R2_001.head40k.fastq aloha1b.R2.fastq
    run bash -c "snakemake -j 10 -s ../../../meta.snake -p --config workflow=qc/sort.rna.snake file_root=aloha1b.R2 -k sort_rna_default_all > sort.rna.log 2>&1"
    [ "$status" -eq 0 ]
}

