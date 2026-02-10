#!/usr/bin/env Rscript

# Required libraries
lapply(c("dplyr", "readr", "stringr", "argparse"), library, character.only = TRUE)

# Argument parser
parser <- ArgumentParser(description = "Process RNA-seq data and calculate enrichment")
parser$add_argument("output_up", type = "character", help = "Upregulated output CSV file name")
parser$add_argument("output_down", type = "character", help = "Downregulated output CSV file name")
args <- parser$parse_args()

path_dir <- getwd()
csv_files <- list.files(path = path_dir, pattern = "*_T4metadata.csv", full.names = TRUE)

# Your desired order
desired_order <- c(
  "Antibiotic response/sensing",
  "Cell division",
  "Cell growth and death",
  "Cellular community",
  "Structural proteins",
  "Membrane transport",
  "Signal transduction",
  "DNA repair and editing",
  "Folding, sorting, degradation",
  "Mobile element/pseudogene",
  "Replication",
  "Transcription",
  "Translation",
  "Unknown",
  "Amino acid metabolism",
  "Capsule metabolism",
  "Carbohydrate metabolism",
  "Cell wall metabolism",
  "Cofactor and vitamin metabolism",
  "Energy metabolism",
  "Lipid metabolism",
  "Nucleotide metabolism",
  "Various metabolism"
)

process_file <- function(file_path, sig_type) {
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  data$Category1[is.na(data$Category1)] <- "Unknown"
  
  categories <- unique(data$Category1)
  
  # Choose significance criteria
  if (sig_type == "up") {
    datasig <- dplyr::filter(data, log2FoldChange > 1 & padj < 0.05)
  } else {
    datasig <- dplyr::filter(data, log2FoldChange < -1 & padj < 0.05)
  }
  
  total <- nrow(data)
  group1 <- nrow(datasig)
  results <- data.frame(Category = character(),
                        Overlap = numeric(),
                        CategorySize = numeric(),
                        PValue = numeric(),
                        stringsAsFactors = FALSE)
  
  for (cat in categories) {
    categoryGenes <- subset(data, Category1 == cat)
    group2 <- nrow(categoryGenes)
    overlap <- length(intersect(datasig$gene_id, categoryGenes$gene_id))
    
    pval <- if (group1 > 0) {
      phyper(overlap - 1, group2, total - group2, group1, lower.tail = FALSE)
    } else {
      1
    }
    
    results <- rbind(results, data.frame(Category = cat,
                                         Overlap = overlap,
                                         CategorySize = group2,
                                         PValue = pval,
                                         stringsAsFactors = FALSE))
  }
  
  results$AdjustedPValue <- p.adjust(results$PValue, method = "BH")
  results$Filename <- basename(file_path)
  
  # Order by your desired order, if present
  results$Category <- factor(results$Category, levels = desired_order)
  results <- results[order(results$Category), ]
  
  return(results)
}

all_results_up <- data.frame()
all_results_down <- data.frame()

for (file_path in csv_files) {
  # Upregulated
  enrichment_up <- process_file(file_path, "up")
  all_results_up <- rbind(all_results_up, enrichment_up)
  
  # Downregulated
  enrichment_down <- process_file(file_path, "down")
  all_results_down <- rbind(all_results_down, enrichment_down)
}

# Output files
output_file_up <- args$output_up
output_file_down <- args$output_down

write.csv(all_results_up, file = file.path(path_dir, output_file_up), row.names = FALSE)
write.csv(all_results_down, file = file.path(path_dir, output_file_down), row.names = FALSE)
cat("Upregulated enrichment done ->", output_file_up, "\n")
cat("Downregulated enrichment done ->", output_file_down, "\n")
