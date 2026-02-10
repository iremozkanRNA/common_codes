#!/usr/bin/env Rscript

# Load necessary libraries
lapply(c("stringr", "dplyr", "readr"), library, character.only = TRUE)

# Set directory path
dir_path <- getwd()

# Collect CSV files
csv_files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)

# Read metadata
metadata <- read.csv("T4metadata.csv")

# Function to process each CSV file
process_file <- function(file_path) {
  
  # Read the CSV file
  data <- read.csv(file_path)
  
  # Remove "gene_" or "gene-" prefix from gene_id
  data <- data %>%
    mutate(gene_id = str_remove(Gene, "^gene[_-]"))
  
  # Merge data with metadata based on gene_id and newlocus
  data <- data %>%
    left_join(metadata, by = c("gene_id" = "Gene"))
  
  # Create new file name with "metadata" added
  new_filename <- str_replace(basename(file_path), "\\.csv$", "_metadata.csv")
  new_filepath <- file.path(dir_path, new_filename)
  
  # Write the processed data to a new CSV file
  write_csv(data, new_filepath)
  
  # Print a message indicating the file has been processed
  cat("Processed:", basename(file_path), "->", new_filename, "\n")
}

# Process all CSV files
lapply(csv_files, process_file)
