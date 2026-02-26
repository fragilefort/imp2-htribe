#!/usr/bin/env bash

BAM_DIR="/home/emadeldin/HTRIBE_ncRNA/notrim_mm10/dedupped_reads"
SAM_DIR="/home/emadeldin/HTRIBE_ncRNA/notrim_mm10/dedupped_reads"


for bam in "$BAM_DIR"/*.bam; do
    filename=$(basename -- "$bam")
    # remove the suffix
    samplename=${filename%.bam}
    samtools view -h -o "$SAM_DIR"/"$samplename".sam "$bam"
done
