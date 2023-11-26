

# Load required libraries
library(rentrez)
library(readxl)

# Read the Excel file
file_path <- "TRPM8_orthologs.xlsx" # Make sure the file is in your working directory or provide the full path
orthologs_data <- read_excel(file_path)

# Extract the RefSeq accession numbers and organism names
refseq_accessions <- orthologs_data$`RefSeq Transcript accessions`
common_names <- orthologs_data$`Common name`

# Function to read the last accessed index from a file
read_last_index <- function(file_name) {
  if (file.exists(file_name)) {
    return(as.numeric(readLines(file_name)))
  } else {
    return(0)
  }
}

# Function to write the last accessed index to a file
write_last_index <- function(file_name, index) {
  writeLines(as.character(index), file_name)
}

# Initialize variables
last_index_file <- "last_index.txt"
last_index <- read_last_index(last_index_file)
fasta_file <- "TRPM8_orthologs_sequences.fasta"

# Check if the fasta file exists and append to it if it does
if(file.exists(fasta_file)) {
  fasta_sequences <- readLines(fasta_file)
} else {
  fasta_sequences <- c()
}

# Loop through the accession numbers and fetch the sequences starting from the last index
for (i in (last_index + 1):length(refseq_accessions)) {
  accession <- refseq_accessions[i]
  organism <- common_names[i]
  
  # Attempt to fetch the sequence from NCBI
  tryCatch({
    fasta_data <- entrez_fetch(db = "nucleotide", id = accession, rettype = "fasta", retmode = "text")
    
    # Add the organism name as a descriptor to the FASTA header
    fasta_header <- sub(">", paste0(">", organism, " | "), fasta_data)
    
    # Append the modified FASTA data to our list
    fasta_sequences <- c(fasta_sequences, fasta_header)
    
    # Write the current sequence to the FASTA file
    writeLines(fasta_header, fasta_file, append = TRUE)
    
    # Optional: print a message to track progress
    cat("Fetched sequence for", organism, "with accession", accession, "\n")
    
    # Update the last index file
    write_last_index(last_index_file, i)
  }, error = function(e) {
    cat("Error fetching sequence for accession", accession, ":", e$message, "\n")
  })
}

writeLines(fasta_sequences, "TRPM8_nucleotide_sequences.fasta")


# Print completion message
cat("FASTA file updated up to the last successful fetch.\n")
