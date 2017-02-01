setup() {
    rm -rf test/scratch/anvio
    mkdir -p test/scratch/anvio
    cd test/scratch/anvio
}

@test "From reaads through megahit to anvio" {
    run bash -c "snakemake -s ../../../anvio.metagenomic.snake --configfile ../../data/configs/anvio.megahit.all.yaml -p -j 20 > anvio.megahit.all.log 2>&1"
    [ "$status" -eq 0 ]
}
