#!/usr/bin/env Rscript 

# Load required libraries
lapply(c("dplyr", "readr", "stringr", "argparse"), library, character.only = TRUE)

# Set up argument parser
parser <- ArgumentParser(description = "Process RNA-seq data and calculate enrichment")
parser$add_argument("output_file", type = "character", help = "Output CSV file name")

# Parse arguments
args <- parser$parse_args()

# Get the current working directory
path_dir <- getwd()

# List all CSV files in the directory with the pattern *_metadata.csv
csv_files <- list.files(path = path_dir, pattern = "*_metadata.csv", full.names = TRUE)

# Function to process each file
process_file <- function(file_path) {
  # Read the .csv file
  data <- read.csv(file_path)
  
  # Assign NA of Category1 as "Unknown"
  data["Category1"][is.na(data["Category1"])] <- "Unknown"
  
  # Get unique categories
  categories <- unique(data$Category1)
  
  # Assign Significance
  datasig <- data %>% filter(
    abs(log2FoldChange) > 1 & padj < 0.05
  )
  
  total <- as.numeric(nrow(data))
  group1 <- as.numeric(nrow(datasig))
  
  # Initialize results dataframe
  results <- data.frame(Category = character(), 
                        Overlap = numeric(), 
                        CategorySize = numeric(), 
                        PValue = numeric(), 
                        stringsAsFactors = FALSE)
  
  # Loop through each category
  for (cat in categories) {
    # Subset data for current category
    categoryGenes <- subset(data, Category1 == cat)
    
    # Calculate overlap
    overlap <- nrow(as.data.frame(intersect(datasig$gene_id, categoryGenes$gene_id)))
    
    # Category size
    group2 <- nrow(categoryGenes)
    
    # Calculate p-value
    pval <- phyper(overlap - 1, group2, total - group2, group1, lower.tail = FALSE)
    
    # Add to results
    results <- rbind(results, data.frame(Category = cat, 
                                         Overlap = overlap, 
                                         CategorySize = group2, 
                                         PValue = pval))
  }
  
  # Adjust p-values for multiple testing
  results$AdjustedPValue <- p.adjust(results$PValue, method = "BH")
  
  # Sort results by p-value
  results <- results[order(results$PValue), ]
  
  # Add filename column
  results$Filename <- basename(file_path)
  
  return(results)
}

# Initialize an empty dataframe to store all results
all_results <- data.frame()

# Process all CSV files and combine results
for (file_path in csv_files) {
  # Process the file
  enrichment_results <- process_file(file_path)
  
  # Append results to all_results
  all_results <- rbind(all_results, enrichment_results)
}

# Write all results to a single CSV file
output_file <- args$output_file
write.csv(all_results, file = file.path(path_dir, output_file), row.names = FALSE)

# Print a message indicating completion
cat("Enrichment analysis complete. Results saved to:", output_file, "\n")