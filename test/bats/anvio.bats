setup() {
    ENV=anvi2
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    ANVIENV=$ENV_DIR
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > /dev/null 2>&1
        source activate $ANVIENV
        anvi-setup-ncbi-cogs --just-do-it > $ANVIENV/.cog.log 2>&1
        source deactivate
    fi

    ENV=assembly
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > /dev/null 2>&1
    fi
    source activate $ENV_DIR

    rm -rf test/scratch/anvio
    mkdir -p test/scratch/anvio
    cd test/scratch/anvio
}

@test "From reaads through megahit to anvio" {
    run bash -c "snakemake -s ../../../anvio.metagenomic.snake --configfile ../../data/configs/anvio.megahit.all.yaml -p -j 20 > anvio.megahit.all.log 2>&1"
    [ "$status" -eq 0 ]
}
