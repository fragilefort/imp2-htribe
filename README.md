# Editing sites summary, code, and multiqc report


### Pre-HTRIBE

Reads were aligned to the mouse reference genome (GRCm39 - GENCODE edition) using the `STAR` aligner (`./code/star_align.py`), and the resulting BAM files were converted to SAM format using `samtools` (`./code/bam2sam.py`). Then, the aligned reads were dedupped using the script `./code/dedup.py`. The multiqc of theses reads are in `./multiqc_report.html`


### Running HyperTRIBE

The SAM files were converted into matrix files and uploaded to a MariaDB database under the table `lastHope` using `./code/sendto_maria.py`. The software then queries this database to compare sample groups and identify editing sites.

The minimum coverage threshold was set to 9 by default. All pairwise comparisons were run using `./code/runall.py`. The primary output files of interest are the bedgraphs. Finally, summarizing the results was done using `./code/annotate_A2G_site.py`. The results for this are in `./htribe_results/`.

### Enrichment analysis
Enrichment analysis was done on the three comparsions: wt-mcherry (background), wt-imp2, mcherry-imp2. This is done usine `./code/holy_enrichment.r` which uses the utils from `./code/enrichment_utils.r`

