# enrichment_utils.r
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(stringr)
library(scales)

symbols_to_entrez <- function(gene_symbols, orgdb = org.Mm.eg.db) {
    map <- AnnotationDbi::select(orgdb,
                                 keys    = unique(as.character(gene_symbols)),
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL")
    unique(na.omit(map$ENTREZID))
}

get_enrichment_data <- function(gene_list, gene_universe = "default", orgdb = org.Mm.eg.db, 
                                ont = "BP", pvalueCutoff = 0.05) {
    if (length(gene_list) == 0) return(NULL)

    gene_entrez <- symbols_to_entrez(gene_list, orgdb)
    if (length(gene_entrez) == 0) {
        warning("No Entrez IDs found for provided gene symbols.")
        return(NULL)
    }
    message(length(gene_entrez), " / ", length(gene_list), " gene symbols mapped to Entrez IDs")

    universe_entrez <- AnnotationDbi::keys(orgdb, keytype = "ENTREZID")

    clusterProfiler::enrichGO(gene          = gene_entrez,
                              universe      = universe_entrez,
                              OrgDb         = orgdb,
                              ont           = ont,
                              pAdjustMethod = "fdr",
                              pvalueCutoff  = pvalueCutoff,
                              qvalueCutoff  = 1,
                              readable      = TRUE)
}

create_enrichment_plot <- function(ego, title = NULL, showCategory = 8) {
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)

    p <- enrichplot::dotplot(ego, showCategory = showCategory) +
        scale_color_gradient(low = "#1a9850", high = "#FFD966", name = "Adjusted p-value") +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) +
        scale_x_continuous(breaks = pretty_breaks(n = 2),
                           expand = expansion(add = c(0.005, 0.005)),
                           name   = "Gene Ratio") +
        theme_light() +
        theme(legend.text  = element_text(size = 12),
              legend.title = element_text(size = 13),
              axis.title   = element_text(size = 13),
              axis.text    = element_text(colour = "black", size = 11),
              plot.title   = element_text(size = 14, face = "bold"),
              plot.margin  = margin(10, 10, 10, 10))

    if (!is.null(title)) p <- p + ggtitle(title)
    return(p)
}

use_enrichment <- function(gene_list, gene_universe = "default", out_prefix = "enrichment_result",
                           ont = "BP", showCategory = 8, pvalueCutoff = 0.05) {
    if (length(gene_list) == 0) {
        message("No genes supplied for: ", out_prefix)
        return(list(enrich = NULL, plot = NULL))
    }

    ego      <- get_enrichment_data(gene_list, gene_universe, ont = ont, pvalueCutoff = pvalueCutoff)
    plot_obj <- create_enrichment_plot(ego, title = basename(out_prefix), showCategory = showCategory)

    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        message("No enriched terms found for: ", out_prefix)
    } else {
        message(nrow(as.data.frame(ego)), " enriched terms found for: ", out_prefix)
        saveRDS(ego, paste0(out_prefix, "_enrich.rds"))
        if (!is.null(plot_obj))
            ggsave(paste0(out_prefix, "_dotplot.png"), plot_obj, width = 10, height = 7, dpi = 200)
    }

    return(list(enrich = ego, plot = plot_obj))
}
