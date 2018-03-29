setup() {
    mkdir -p test/conda/envs
    ENV=catalog
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi
    # suppress warnings
    source activate $ENV_DIR 2>/dev/null || true
    # quit if activation failed
    if [ "$CONDA_PREFIX" != "$ENV_DIR" ]; then
        echo "ERROR activating cond env $ENV_DIR"
        exit 1
    fi
}


@test "catalog: create empty folders" {
    rm -rf test/scratch/catalog

    mkdir -p test/scratch/catalog/ALOHA_04
    mkdir -p test/scratch/catalog/ALOHA_15
    mkdir -p test/scratch/catalog/ALOHA_21
}

@test "assemble sample ALOHA_04" {
    cd test/scratch/catalog/ALOHA_04
    run bash -c "snakemake -j 10 -s ../../../../assembly.metagenomic.snake -p --configfile ../../../data/configs/megahit.A04.yaml > assembly.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "assemble samples ALOHA_15" {
    cd test/scratch/catalog/ALOHA_15
    run bash -c "snakemake -j 10 -s ../../../../assembly.metagenomic.snake -p --configfile ../../../data/configs/megahit.A15.yaml > assembly.log 2>&1"
    [ "$status" -eq 0 ]
}    

@test "assemble samples ALOHA_21" {
    cd test/scratch/catalog/ALOHA_21
    run bash -c "snakemake -j 10 -s ../../../../assembly.metagenomic.snake -p --configfile ../../../data/configs/megahit.A21.yaml > assembly.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Compile gene catalog" {
    cd test/scratch/catalog
    run bash -c "snakemake -s ../../../annotation.gene_catalog.snake --configfile ../../data/configs/gene_catalog.yaml -p -k -j 20 --notemp > gene_catalog.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../annotation.gene_catalog.snake --configfile ../../data/configs/gene_catalog.yaml -p -k -j 20 -n 2>&1 | grep '^Nothing to be done'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}
