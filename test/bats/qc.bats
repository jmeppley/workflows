setup() {
    ENV=illumina.qc
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet
    fi
    source activate $ENV_DIR

    rm -rf test/scratch/qc
    mkdir -p test/scratch/qc
    cd test/scratch/qc
}

@test "Prep a pair of scriptseq files with pear" {
    ln -s ../../data/raw_reads/2014_ALOHA_XVII_1-1B_S1_R1_001.head40k.fastq aloha1b.R1.fastq
    ln -s ../../data/raw_reads/2014_ALOHA_XVII_1-1B_S1_R2_001.head40k.fastq aloha1b.R2.fastq
    run bash -c "snakemake -j 10 -s ../../../test.snake -p --config workflow=qc/pear.snake --verbose aloha1b.scripseq.ATCACG.trim_adapt.joined.fastq > pear.join.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Prep two pair of scriptseq files with flash" {
    run bash -c "snakemake -j 10 -s ../../../clean.illumina.snake -p --configfile ../../data/configs/illumina.yaml --verbose --config joining_program=flash run_stats=True > flash.join.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Prep two pair of scriptseq files with pandaseq" {
    which pandaseq || skip "Pandaseq is not installed"
    run bash -c "snakemake -j 10 -s ../../../clean.illumina.snake -p --configfile ../../data/configs/illumina.yaml --verbose --config joining_program=pandaseq > pandaseq.join.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Prep two pair of files for assembly" {
    run bash -c "snakemake -j 10 -s ../../../clean.illumina.snake -p --configfile ../../data/configs/illumina.yaml --verbose --config cleaning_protocol=assembly > assembly.join.log 2>&1"
    [ "$status" -eq 0 ]
}

