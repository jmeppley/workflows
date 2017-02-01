setup() {
    rm -rf test/scratch/anvio.contig
    mkdir -p test/scratch/anvio.contig
    cd test/scratch/anvio.contig
}

@test "From reaads and contigs into anvio" {
    run bash -c "snakemake -s ../../../anvio.metagenomic.snake --configfile ../../data/configs/anvio.contigs.all.yaml -p -j 20 > anvio.contigs.all.log 2>&1"
    [ "$status" -eq 0 ]
}
