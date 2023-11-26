library(rentrez)
library(readxl)

# Read the Excel file
file_path <- "TRPM8_orthologs.xlsx" # Make sure the file is in your working directory or provide the full path
orthologs_data <- read_excel(file_path)

# Extract the RefSeq protein accession numbers and organism names
refseq_accessions <- orthologs_data$`RefSeq Protein accessions`
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
fasta_file <- "TRPM8_protein_sequences.fasta"
fasta_sequences <- "" # Variable to accumulate the FASTA sequences

# Loop through the accession numbers and fetch the sequences starting from the last index
for (i in (last_index + 1):length(refseq_accessions)) {
  accession <- refseq_accessions[i]
  organism <- common_names[i]
  
  # Attempt to fetch the protein sequence from NCBI
  tryCatch({
    fasta_data <- entrez_fetch(db = "protein", id = accession, rettype = "fasta", retmode = "text")
    
    # Add the organism name as a descriptor to the FASTA header
    fasta_header <- sub(">", paste0(">", organism, " | "), fasta_data)
    fasta_sequences <- paste0(fasta_sequences, fasta_header, "\n") # Append to the accumulator
    
    # Write the current sequence to the FASTA file
    base::writeLines(fasta_header, fasta_file, append = TRUE)
    
    # Optional: print a message to track progress
    cat("Fetched protein sequence for", organism, "with accession", accession, "\n")
    
    # Update the last index file
    write_last_index(last_index_file, i)
  }, error = function(e) {
    cat("Error fetching protein sequence for accession", accession, ":", e$message, "\n")
  })
}

# Write the accumulated sequences to a file at the end
writeLines(fasta_sequences, "TRPM8_AminoAcid_sequences.fasta")
