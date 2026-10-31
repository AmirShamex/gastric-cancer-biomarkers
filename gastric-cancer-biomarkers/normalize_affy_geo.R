#!/usr/bin/env Rscript
# ==============================================================================
# Microarray Expression Preprocessing & RMA Normalization Pipeline
# Associated with: Second International Congress of Cancer Genomics (CGC2024)
# Reference Study: Unveiling Potential Biomarkers for Gastric Cancer Through
#                  Bioinformatics Approaches (Certificate No: 40321102)
# ==============================================================================

suppressPackageStartupMessages({
  library(affy)
  library(GEOquery)
  library(tidyverse)
})

run_normalization <- function(geo_id = "GSE79973",
                              data_dir = "01_transcriptomics_microarray/data",
                              output_dir = "01_transcriptomics_microarray/results") {
  
  message(sprintf("[INFO] Initializing preprocessing for GEO dataset: %s", geo_id))
  
  raw_tar <- file.path(data_dir, paste0(geo_id, "_RAW.tar"))
  extracted_dir <- file.path(data_dir, paste0(geo_id, "_extracted"))
  
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Step 1: Extract CEL archive if archive exists and hasn't been extracted
  if (file.exists(raw_tar) && !dir.exists(extracted_dir)) {
    message("[INFO] Extracting raw CEL archive...")
    untar(tarfile = raw_tar, exdir = extracted_dir)
  }
  
  # Step 2: Read Affymetrix arrays and execute RMA normalization
  if (!dir.exists(extracted_dir)) {
    stop(sprintf("[ERROR] Neither archive '%s' nor directory '%s' was found.", raw_tar, extracted_dir))
  }
  
  message("[INFO] Reading raw probe intensities and performing RMA normalization...")
  raw_data <- ReadAffy(celfile.path = extracted_dir)
  normalized_data <- rma(raw_data)
  expr_matrix <- as.data.frame(exprs(normalized_data))
  message(sprintf("[INFO] RMA normalization complete: %d probes x %d samples.", 
                  nrow(expr_matrix), ncol(expr_matrix)))
  
  # Step 3: Retrieve platform feature annotations via GEOquery
  message(sprintf("[INFO] Querying GEO platform metadata for %s...", geo_id))
  gse <- getGEO(geo_id, GSEMatrix = TRUE)
  feature_data <- gse[[1]]@featureData@data
  
  # Extract Probe ID and primary Gene Symbol columns
  annot_subset <- feature_data[, c(1, min(11, ncol(feature_data)))]
  
  # Step 4: Merge annotation with normalized expression values
  message("[INFO] Assembling annotated expression matrix...")
  combined_data <- cbind(annot_subset, expr_matrix)
  
  # Step 5: Export results
  output_path <- file.path(output_dir, paste0(geo_id, "_normalized_annotated.tsv"))
  write_tsv(combined_data, file = output_path)
  message(sprintf("[SUCCESS] Annotated matrix exported to: %s", output_path))
}

# Command-line execution entry point
args <- commandArgs(trailingOnly = TRUE)
target_geo <- ifelse(length(args) > 0, args[1], "GSE79973")
run_normalization(geo_id = target_geo)