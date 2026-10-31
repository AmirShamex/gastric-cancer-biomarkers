#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: Paired Differential Expression Analysis & Visualization (limma)
# Associated Study: Biomarker Discovery in Gastric Adenocarcinoma (GSE79973)
# Congress: Second International Congress of Cancer Genomics (CGC2024)
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
})

run_deg_pipeline <- function(
    input_file = "01_transcriptomics_microarray/results/GSE79973_normalized_annotated.tsv",
    output_dir = "01_transcriptomics_microarray/results",
    log2fc_cutoff = 1.5,
    adj_p_cutoff = 0.05
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  message(sprintf("[INFO] Loading normalized dataset from: %s", input_file))
  
  data <- read_tsv(input_file, show_col_types = FALSE)
  sample_cols <- grep("^GSM", colnames(data), value = TRUE)
  
  # Separate cancer (C) and normal (N) paired samples
  cancer_cols <- sort(grep("C_HG", sample_cols, value = TRUE))
  normal_cols <- sort(grep("N_HG", sample_cols, value = TRUE))
  
  cancer_ids <- sub(".*_(13[0-9]+)C_.*", "\\1", cancer_cols)
  normal_ids <- sub(".*_(13[0-9]+)N_.*", "\\1", normal_cols)
  
  stopifnot(identical(cancer_ids, normal_ids))
  n_pairs <- length(cancer_ids)
  message(sprintf("[INFO] Verified %d matched patient pairs.", n_pairs))
  
  # Construct paired design matrix
  ordered_samples <- c(normal_cols, cancer_cols)
  expr_matrix <- as.matrix(data[, ordered_samples])
  rownames(expr_matrix) <- data$ID
  
  patient <- factor(rep(cancer_ids, 2))
  tissue <- factor(rep(c("Normal", "Tumor"), each = n_pairs), levels = c("Normal", "Tumor"))
  
  design <- model.matrix(~ patient + tissue)
  colnames(design)[ncol(design)] <- "Tumor_vs_Normal"
  
  # Fit linear model with empirical Bayes moderation
  message("[INFO] Fitting paired linear model and applying eBayes moderation...")
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  
  deg_table <- topTable(fit, coef = "Tumor_vs_Normal", number = Inf, sort.by = "P") %>%
    rownames_to_column(var = "ID") %>%
    left_join(data[, c("ID", "Gene Symbol")], by = "ID") %>%
    rename(Gene_Symbol = `Gene Symbol`)
  
  # Assign significance classification
  deg_table <- deg_table %>%
    mutate(
      Significance = case_when(
        logFC >= log2fc_cutoff & adj.P.Val < adj_p_cutoff ~ "Upregulated (Tumor)",
        logFC <= -log2fc_cutoff & adj.P.Val < adj_p_cutoff ~ "Downregulated (Tumor)",
        TRUE ~ "Not Significant"
      )
    )
  
  # Export full statistics table
  deg_out_file <- file.path(output_dir, "GSE79973_DEG_limma_results.tsv")
  write_tsv(deg_table, deg_out_file)
  message(sprintf("[SUCCESS] DEG statistical results exported to: %s", deg_out_file))
  
  # --------------------------------------------------------------------------
  # 1. Publication-Quality Volcano Plot
  # --------------------------------------------------------------------------
  message("[INFO] Generating publication-grade Volcano Plot...")
  
  key_biomarkers <- c("INHBA", "COL10A1", "FNDC1", "COL11A1", "SFRP4",
                      "GKN2", "GKN1", "ATP4B", "GIF", "KIF26B", "PIWIL2")
  
  volcano_labels <- deg_table %>%
    filter(Gene_Symbol %in% key_biomarkers) %>%
    group_by(Gene_Symbol) %>%
    slice_min(P.Value, n = 1) %>%
    ungroup()
  
  p_volcano <- ggplot(deg_table, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
    geom_point(alpha = 0.45, size = 1.4) +
    scale_color_manual(values = c(
      "Upregulated (Tumor)" = "#D9381E",
      "Downregulated (Tumor)" = "#2B5B84",
      "Not Significant" = "#CCCCCC"
    )) +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", linewidth = 0.4) +
    geom_point(data = volcano_labels, color = "black", shape = 1, size = 3, stroke = 1.1) +
    geom_text_repel(
      data = volcano_labels,
      aes(label = Gene_Symbol),
      color = "black",
      fontface = "bold",
      size = 3.6,
      max.overlaps = 20,
      box.padding = 0.4
    ) +
    labs(
      title = "Volcano Plot: Gastric Adenocarcinoma vs. Normal Mucosa (GSE79973)",
      subtitle = sprintf("Paired eBayes Moderated Analysis (FDR < %.2f, |log2FC| >= %.1f)", adj_p_cutoff, log2fc_cutoff),
      x = expression(bold(log[2]~"Fold Change (Tumor / Normal)")),
      y = expression(bold(-log[10]~italic(p)*"-value"))
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      legend.position = "top",
      legend.title = element_blank()
    )
  
  volcano_png <- file.path(output_dir, "volcano_plot_GSE79973.png")
  ggsave(volcano_png, plot = p_volcano, width = 8.5, height = 6.2, dpi = 300)
  message(sprintf("[SUCCESS] Volcano plot saved to: %s", volcano_png))
  
  # --------------------------------------------------------------------------
  # 2. Clustered Heatmap of Top 30 Significant Biomarkers
  # --------------------------------------------------------------------------
  message("[INFO] Generating Clustered Heatmap...")
  
  top_up <- deg_table %>%
    filter(Significance == "Upregulated (Tumor)", !is.na(Gene_Symbol)) %>%
    arrange(desc(logFC)) %>%
    distinct(Gene_Symbol, .keep_all = TRUE) %>%
    head(15)
    
  top_down <- deg_table %>%
    filter(Significance == "Downregulated (Tumor)", !is.na(Gene_Symbol)) %>%
    arrange(logFC) %>%
    distinct(Gene_Symbol, .keep_all = TRUE) %>%
    head(15)
    
  top_candidates <- bind_rows(top_up, top_down)
  
  # Interleave pairs side-by-side: P1_N, P1_T, P2_N, P2_T, ...
  interleaved_samples <- as.vector(rbind(normal_cols, cancer_cols))
  interleaved_labels <- as.vector(rbind(paste0("P", cancer_ids, "_N"), paste0("P", cancer_ids, "_T")))
  
  heatmap_mat <- expr_matrix[top_candidates$ID, interleaved_samples]
  rownames(heatmap_mat) <- top_candidates$Gene_Symbol
  colnames(heatmap_mat) <- interleaved_labels
  
  # Row-wise Z-score transformation
  z_matrix <- t(scale(t(heatmap_mat)))
  
  annotation_col <- data.frame(
    Tissue = rep(c("Normal", "Tumor"), n_pairs),
    row.names = interleaved_labels
  )
  ann_colors <- list(Tissue = c(Normal = "#2B5B84", Tumor = "#D9381E"))
  
  heatmap_png <- file.path(output_dir, "heatmap_top_DEGs_GSE79973.png")
  pheatmap(
    z_matrix,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    scale = "none",
    color = colorRampPalette(c("#2B5B84", "#F7F7F7", "#D9381E"))(100),
    breaks = seq(-2.5, 2.5, length.out = 101),
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    fontsize_row = 8,
    fontsize_col = 8,
    main = "Expression Profiles of Top Gastric Cancer Biomarkers (10 Matched Pairs)",
    filename = heatmap_png,
    width = 9.5,
    height = 7.5
  )
  message(sprintf("[SUCCESS] Heatmap saved to: %s", heatmap_png))
}

# Run
run_deg_pipeline()