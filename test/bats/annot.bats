setup() {
    eval "$(conda shell.bash hook)"
    mkdir -p test/conda/envs
    ENV=annotate
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi

    conda activate $ENV_DIR
}

@test "Annotate genes" {
    rm -rf test/scratch/annot
    mkdir -p test/scratch/annot
    cd test/scratch/annot
    run bash -c "snakemake -s ../../../annotation.genes.snake --configfile ../../data/configs/genes.annot.yaml -p -k -j 2 > genes.annot.prod.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../annotation.genes.snake --configfile ../../data/configs/genes.annot.yaml -p -k -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "Annotate existing gene catalog" {
    rm -rf test/scratch/annot_cat
    mkdir -p test/scratch/annot_cat
    cd test/scratch/annot_cat
    run bash -c "snakemake -s ../../../annotation.gene_catalog.snake --configfile ../../data/configs/genes.annot.yaml -p -k -j 2 --notemp > genes.annot.prod.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../annotation.gene_catalog.snake --configfile ../../data/configs/genes.annot.yaml -p -k -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}
