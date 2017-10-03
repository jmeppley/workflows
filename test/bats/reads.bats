setup() {
    mkdir -p test/conda/envs
    ENV=annotate
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > /dev/null 2>&1
    fi
    source activate $ENV_DIR
}

@test "Annotate reads from a bam file with prodigal" {
    rm -rf test/scratch/reads.prod
    mkdir -p test/scratch/reads.prod
    cd test/scratch/reads.prod
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.prod.yaml -p -k -j 20 > reads.annot.prod.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Annotate reads from a bam file with sixframe translation" {
    rm -rf test/scratch/reads.six
    mkdir -p test/scratch/reads.six
    cd test/scratch/reads.six
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.six.yaml -p -k -j 20 > reads.annot.six.log 2>&1"
    [ "$status" -eq 0 ]
}

