setup() {
    rm -rf test/scratch/anvio
    mkdir -p test/scratch/anvio
    cd test/scratch/anvio
}

@test "From reaads through megahit to anvio" {
    run snakemake -s ../../../anvio.metagenomic.snake --configfile ../../data/configs/anvio.megahit.all.yaml -p -j 20
    [ "$status" -eq 0 ]
}
