setup() {
    ENV=anvi2
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > /dev/null 2>&1
    fi
    ANVIENV=$ENV_DIR
    ENV=assembly
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > /dev/null 2>&1
    fi
    source activate $ENV_DIR
    rm -rf test/scratch/anvio.contig
    mkdir -p test/scratch/anvio.contig
    cd test/scratch/anvio.contig
}

@test "From reaads and contigs into anvio" {
    echo $ANVIENV > test.anvi.env.txt
    echo $PATH >> test.anvi.env.txt
    run bash -c "snakemake -s ../../../anvio.metagenomic.snake --configfile ../../data/configs/anvio.contigs.all.yaml -p -j 20 --config anvio_env=$ANVIENV --verbose > anvio.contigs.all.log 2>&1"
    [ "$status" -eq 0 ]
}
