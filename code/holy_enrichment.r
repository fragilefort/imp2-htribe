# This script should do the following:
# Remove the genes in the wt-imp2, mcherry=imp2 that are in the background
# Perform enrichment analysis for those genes:universe is all mm10 genes
# It also perform enrichment for background (wt-mcheryy) using the same universe
# I will also try to use the standard universe on cluster profiler
# Maybe we can also try using the mm39 genes to see any difference
# environment name is clusterprofiler


library(here)
library(data.table)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(stringr)
library(scales)
library(grid)

setwd("/home/emadeldin/work/hiwi/HTRIBE_ncRNA/grcm39_Feb26/")
data_dir <- file.path(here::here(), "data")
out_dir  <- file.path(here::here(), "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(here::here(), "code", "enrichment_utils.r"))

read_file <- function(file_path) {
    fread(file_path,
          sep = "\t",
          header = TRUE,
          encoding = "UTF-8",
          showProgress = FALSE,
          data.table = TRUE)
}

wt_imp2_df      <- read_file(file.path(data_dir, "wt_imp2_A2G_1pct.xls"))
mcherry_imp2_df <- read_file(file.path(data_dir, "mcherry_imp2_A2G_1pct.xls"))
background_df   <- read_file(file.path(data_dir, "wt_mcherry_A2G_1pct.xls"))

wt_imp2_genes      <- unique(wt_imp2_df$Gene_name)
mcherry_imp2_genes <- unique(mcherry_imp2_df$Gene_name)
background_clean   <- unique(background_df$Gene_name)

remove_bkg <- function(genes, bkg) {
    if (length(genes) == 0) return(character(0))
    setdiff(unique(genes), unique(bkg))
}

wt_imp2_clean      <- remove_bkg(wt_imp2_genes, background_clean)
mcherry_imp2_clean <- remove_bkg(mcherry_imp2_genes, background_clean)

cat("wt_imp2 genes after background removal:     ", length(wt_imp2_clean), "\n")
cat("mcherry_imp2 genes after background removal:", length(mcherry_imp2_clean), "\n")
cat("background genes:                           ", length(background_clean), "\n")

wt_imp2_res      <- use_enrichment(wt_imp2_clean,      
                                   gene_universe = "default",
                                   out_prefix = file.path(out_dir, "wt_imp2_1pct"))

mcherry_imp2_res <- use_enrichment(mcherry_imp2_clean, 
                                   gene_universe = "default",
                                   out_prefix = file.path(out_dir, "mcherry_imp2_1pct"))

wt_mcherry_res   <- use_enrichment(background_clean,   
                                   gene_universe = "default",
                                   out_prefix = file.path(out_dir, "wt_mcherry_1pct"))

save_plot <- function(res, filename, width = 10, height = 10) {
    if (!is.null(res$plot)) {
        png(file.path(out_dir, filename),
            width = width, height = height, units = "in", res = 300)
        print(res$plot)
        dev.off()
        message("Saved: ", filename)
    } else {
        message("No significant enrichment to plot for: ", filename)
    }
}

save_plot(wt_imp2_res,      "wt_imp2_1pct.png")
save_plot(mcherry_imp2_res, "mcherry_imp2_1pct.png")
save_plot(wt_mcherry_res,   "wt_mcherry_1pct.png")

plots_available <- list()

if (!is.null(wt_imp2_res$plot))      plots_available[["WT vs. IMP2"]]      <- wt_imp2_res$plot
if (!is.null(mcherry_imp2_res$plot)) plots_available[["mCherry vs. IMP2"]] <- mcherry_imp2_res$plot
if (!is.null(wt_mcherry_res$plot))   plots_available[["WT vs. mCherry"]]   <- wt_mcherry_res$plot

if (length(plots_available) > 0) {
    png(file.path(out_dir, "GO_annots_combined.png"),
        width = 18, height = 10, units = "in", res = 300)
    print(
        ggarrange(plotlist = plots_available,
                  ncol = 2,
                  nrow = ceiling(length(plots_available) / 2),
                  common.legend = TRUE,
                  legend = "bottom")
    )
    dev.off()
    message("Saved: GO_annots_combined.png")
} else {
    message("No plots available for combined figure — no significant enrichment found in any comparison.")
}
