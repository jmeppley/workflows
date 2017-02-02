setup() {
    rm -rf test/scratch/reads
    mkdir -p test/scratch/reads
    cd test/scratch/reads
}

@test "Annotate reads from a bam file" {
    run bash -c "snakemake -s ../../../annotation.reads.snake --configfile ../../data/configs/reads.annot.yaml -p -k -j 20 > reads.annot.log 2>&1"
    [ "$status" -eq 0 ]
}
