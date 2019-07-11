setup() {
    eval "$(conda shell.bash hook)"
    mkdir -p test/conda/envs
    ENV=annotate
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > /dev/null 2>&1
    fi
    conda activate $ENV_DIR
}

@test "Annotate reads from a bam file with prodigal" {
    rm -rf test/scratch/reads.prod
    mkdir -p test/scratch/reads.prod
    cd test/scratch/reads.prod
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.prod.yaml -p -k -j 2 > reads.annot.prod.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.prod.yaml -p -k -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "Annotate reads from a bam file with sixframe translation" {
    rm -rf test/scratch/reads.six
    mkdir -p test/scratch/reads.six
    cd test/scratch/reads.six
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.six.yaml -p -k -j 2 > reads.annot.six.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.six.yaml -p -k -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "Annotate genes from multiple faa files" {
    rm -rf test/scratch/reads.faa
    mkdir -p test/scratch/reads.faa
    cd test/scratch/reads.faa
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.faa.yaml -p -k -j 2 > reads.annot.faa.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.faa.yaml -p -k -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "Annotate reads from multiple fastq files" {
    rm -rf test/scratch/reads.glob
    mkdir -p test/scratch/reads.glob
    cd test/scratch/reads.glob
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.glob.yaml -p -k -j 2 > reads.annot.glob.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.glob.yaml -p -k -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

