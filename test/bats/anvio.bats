setup() {
    ENV=anvio
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    ANVIENV=$ENV_DIR
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf "$ENV_DIR"
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
        source activate $ANVIENV 2>/dev/null || true
        if [ "$CONDA_PREFIX" != "$ENV_DIR" ]; then
            echo "ERROR activating cond env $ENV_DIR"
            exit 1
        fi
        anvi-setup-ncbi-cogs --just-do-it > $ANVIENV/.cog.log 2>&1
        source deactivate
    fi

    # activate for testing below
    source activate $ANVIENV 2>/dev/null || true
}

@test "From reads through megahit to anvio" {
    which megahit > /dev/null 2>&1 || skip "Megahit is not installed" 

    rm -rf test/scratch/anvio
    mkdir -p test/scratch/anvio
    cd test/scratch/anvio

    run bash -c "snakemake -s ../../../anvio.metagenomic.snake --configfile ../../data/configs/anvio.megahit.all.yaml -p -j 20 > anvio.megahit.all.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "From reads and contigs into anvio" {
    rm -rf test/scratch/anvio.contig
    mkdir -p test/scratch/anvio.contig
    cd test/scratch/anvio.contig

    run bash -c "snakemake -s ../../../anvio.metagenomic.snake --configfile ../../data/configs/anvio.contigs.all.yaml -p -j 20 --verbose > anvio.contigs.all.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "From subset of reads and contigs into anvio" {
    cd test/scratch/anvio.contig

    run bash -c "snakemake -s ../../../anvio.metagenomic.snake --configfile ../../data/configs/anvio.contigs.teens.yaml -p -j 20 --verbose > anvio.contigs.teens.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "From explicit list of and contigs into anvio" {
    cd test/scratch/anvio.contig
    rm contigs.raw.fasta
    cp ../../data/contigs/contigs.aloha.fa contigs.raw.fasta

    run bash -c "snakemake -s ../../../anvio.metagenomic.snake --configfile ../../data/configs/anvio.contigs.explicit.yaml -p -j 20 --verbose > anvio.contigs.explicit.log 2>&1"
    [ "$status" -eq 0 ]
}



@test "Pangenome" {
    rm -rf test/scratch/pang
    mkdir -p test/scratch/pang
    cd test/scratch/pang

    run bash -c "snakemake -s ../../../anvio.pangenome.snake --configfile ../../data/configs/anvio.pangenome.yaml -p -j 20 --verbose > anvio.pangenome.log 2>&1"
    [ "$status" -eq 0 ]
}


