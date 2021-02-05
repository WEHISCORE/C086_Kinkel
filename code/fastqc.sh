#!/usr/bin/env bash
# Run FastQC (and MultiQC on the output).
# Peter Hickey
# 2021-02-05

module load fastqc

OUTDIR="../output/FastQC"
mkdir -p ${OUTDIR}

fastqc -o ${OUTDIR} \
       --threads 8 \
       ../extdata/NN206/Fastq/*Kinkel*fastq.gz

multiqc --title C086_Kinkel \
        --outdir ${OUTDIR} \
        --no-megaqc-upload \
        ${OUTDIR}
