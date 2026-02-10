#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

# Get the current working directory
dir_path <- getwd()

# List all CSV files in the directory
csv_files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)

# Function to process each CSV file
process_file <- function(file_path) {
  # Read the CSV file
  data <- read_csv(file_path)
  
  # Ensure the first column is named "gene_id"
  names(data)[1] <- "gene_id"
  
  # Remove "gene_" or "gene-" prefix from gene_id
  data <- data %>%
    mutate(gene_id = str_remove(gene_id, "^gene[_-]"))
  
  # Assign significance based on log2FoldChange and adjusted p-value
  data <- data %>%
    mutate(significance = case_when(
      is.na(log2FoldChange) ~ "insignificant",
      !is.numeric(log2FoldChange) ~ "insignificant",
      abs(log2FoldChange) > 1 & padj < 0.05 ~ "significant",
      TRUE ~ "insignificant"
    ))
  
  # Arrange  significant first, then by absolute log2FoldChange in descending order
  data <- data %>%
    arrange(desc(significance == "significant"),
            desc(abs(log2FoldChange)))
  
  # Create new filename with "processed" added
  new_filename <- str_replace(basename(file_path), "\\.csv$", "_processed.csv")
  new_filepath <- file.path(dir_path, new_filename)
  
  # Write processed data to new CSV file
  write_csv(data, new_filepath)
  
  cat("Processed:", basename(file_path), "->", new_filename, "\n")
  
  # Print summary of significance categories
  cat("Significance summary:\n")
  print(table(data$significance))
  cat("\n")
}

# Process all CSV files in the directory
lapply(csv_files, process_file)

cat("All files have been processed.\n")
