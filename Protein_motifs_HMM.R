install.packages("readtext")
library(readtext)
library(HMM)

# Function to read a Word document and extract the text
read_fasta_from_word <- function(file_path) {
  # Read the Word document using the readtext package
  doc_text <- readtext(file_path)$text
  
  # Assuming each sequence in the Word document is separated by two newline characters
  # and each sequence starts with a '>'
  sequences <- unlist(strsplit(doc_text, "\\n\\n"))
  sequences <- sapply(sequences, function(x) gsub("\\r|\\n", "", x)) # Remove all newlines
  
  return(sequences)
}

# Function to find 'QNE' motif and encode the sequence as states
encode_sequence <- function(sequence) {
  states <- rep("OTHER", nchar(sequence))
  motif <- "QNE"
  motif_length <- nchar(motif)
  
  for (i in 1:(nchar(sequence) - motif_length + 1)) {
    window <- substr(sequence, i, i + motif_length - 1)
    if (window == motif) {
      states[i:(i + motif_length - 1)] <- "QNE"
    }
  }
  return(states)
}

# Read the Word documents containing the training data
mammal_sequences <- read_fasta_from_word("mammals_training_1.docx")
non_mammal_sequences <- read_fasta_from_word("nonmammals_training_1.docx")


# Encode all your sequences
mammal_encoded_sequences <- sapply(mammal_sequences, encode_sequence)
non_mammal_encoded_sequences <- sapply(non_mammal_sequences, encode_sequence)

# Concatenate all states from all sequences to create a 'super sequence' for HMM training
mammal_super_sequence <- unlist(mammal_encoded_sequences)
non_mammal_super_sequence <- unlist(non_mammal_encoded_sequences)

# Define the possible states and observations
states <- c("QNE", "OTHER")
amino_acids <- c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

# Initialize the HMM parameters. This is a simplified example where
# the transition and emission probabilities are set arbitrarily.
# In practice, these should be informed by your data.

transition_probs <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)
emission_probs <- matrix(runif(40), nrow = 2)

# Initialize the HMMs for mammals and non-mammals
hmm_mammal <- initHMM(states, amino_acids, transProbs = transition_probs, emissionProbs = emission_probs)
hmm_non_mammal <- initHMM(states, amino_acids, transProbs = transition_probs, emissionProbs = emission_probs)


# Function to perform Viterbi training
viterbi_training <- function(sequences, hmm) {
  for (sequence in sequences) {
    # Perform Viterbi algorithm to find the most likely state sequence for the given observation sequence
    viterbi_path <- viterbi(hmm, sequence)
    
    
  }
  
  
  return(hmm)
}


# Perform Viterbi training on the mammal and non-mammal sequences
hmm_mammal <- viterbi_training(mammal_encoded_sequences, hmm_mammal)
hmm_non_mammal <- viterbi_training(non_mammal_encoded_sequences, hmm_non_mammal)


encode_test_sequence <- function(test_sequence) {
  states <- rep("OTHER", nchar(test_sequence))
  motif <- "QNE"
  motif_length <- nchar(motif)
  
  for (i in 1:(nchar(test_sequence) - motif_length + 1)) {
    window <- substr(test_sequence, i, i + motif_length - 1)
    if (window == motif) {
      states[i:(i + motif_length - 1)] <- "QNE"
    }
  }
  return(states)
}


# This function takes a sequence and an HMM, predicts the state sequence, and returns the results.
predict_with_hmm <- function(encoded_sequence, hmm_model) {
  viterbi_path <- viterbi(hmm_model, encoded_sequence)
  return(viterbi_path$states) # or return the full viterbi_path if you need more details
}

# Function to process and predict states for new test sequences
test_new_data <- function(new_sequences, hmm_mammal, hmm_non_mammal) {
  encoded_sequences <- sapply(new_sequences, encode_test_sequence)
  predictions_mammal <- lapply(encoded_sequences, predict_with_hmm, hmm_model = hmm_mammal)
  predictions_non_mammal <- lapply(encoded_sequences, predict_with_hmm, hmm_model = hmm_non_mammal)
  
  # Added print statements for debugging
  print("Encoded Sequences:")
  print(encoded_sequences)
  print("Mammal Predictions:")
  print(predictions_mammal)
  print("Non-Mammal Predictions:")
  print(predictions_non_mammal)
  
  return(list(mammal = predictions_mammal, non_mammal = predictions_non_mammal))
}

# Example usage with new test sequences
new_test_sequences <- c("MLFQVSMGTMRHRRNGNFESSRLLYSSMSRSIDVACSDADLANFIQENFKKRECVFFTKDTKSMGNLCKC
GYPENQHIEGTQVNTSEKWNYKKHTKELPTDAFGDIQFENLGKRGKYIRLSCDTDSETLYDLMTQHWHLK
TPNLVISVTGGAKNFALKPRMRKIFSRLIYIAQSKGAWIFTGGTHYGLMKYIGEVVRDNTISRSSEENVV
AIGIAAWGMISNRETLIRTADSDGNYLAHYIMDDLKRDPLYCLDNNHTHLLLVDNGTHGHPTIEAKVRTQ
LEKYISERVIPESNYGGKISINVAIKSKIPCVVVEGSGRIADVIASLMEAEGTLASSCVKESLLRYLPRT
ISRLSEEETESWIKWIKEVLENPHLLTVIKIEEAGDEIVSNAISFALYKAFSTNEHDRDNWNGQLKLLLE
WNQLELASDEIFTNDRNWE"
                        , "MSFRAARLSMRNRRNDTLDSTRTLYSSASRSTDLSYSESDLVNFIQANFKKRECVFFTKDSKATENVCKC
GYAQSQHMEGTQINQSEKWNYKKHTKEFPTDAFGDIQFETLGKKGKYIRLSCDTDAEILYELLTQHWHLK
TPNLVISVTGGAKNFALKPRMRKIFSRLIYIAQSKGAWILTGGTHYGLMKYIGEVVRDNTISRSSEENIV
AIGIAAWGMVSNRDTLIRNCDAEGYFLAQYLMDDFTRDPLYILDNNHTHLLLVDNGCHGHPTVEAKLRNQ
LEKYISERTIQDSNYGGKIPIVCFAQGGGKETLKAINTSIKNKIPCVVVEGSGQIADVIASLVEVEDALT
SSAVKEKLVRFLPRTVSRLPEEETESWIKWLKEILECSHLLTVIKMEEAGDEIVSNAISYALYKAFSTNE
QDKDNWNGQLKLLLEWNQLDLANDEIFTNDRRWESADLQEVMFTALIKDRPKFVRLFLENGLNLRKFLTH
DVLTELFSNHFSTLVYRNLQIAKNSYNDALLTFVWKLVANFRRGFRKEDRNGRDEMDIELHDVSPITRHP
LQALFIWAILQNKKELSKVIWEQTRGCTLAALGASKLLKTLAKVKNDINAAGESEELANEYETRAVELFT
ECYSSDEDLAEQLLVYSCEAWGGSNCLELAVEATDQHFIAQPGVQNFLSKQWYGEISRDTKNWKIILCLF
IIPLVGCGFVSFRKKPVDKHKKLLWYYVAFFTSPFVVFSWNVVFYIAFLLLFAYVLLMDFHSVPHPPELV
LYSLVFVLFCDEVRQWYVNGVNYFTDLWNVMDTLGLFYFIAGIVFRLHSSNKSSLYSGRVIFCLDYIIFT
LRLIHIFTVSRNLGPKIIMLQRMLIDVFFFLFLFAVWMVAFGVARQGILRQNEQRWRWIFRSVIYEPYLA
MFGQVPSDVDGTTYDFAHCTFTGNESKPLCVELDEHNLPRFPEWITIPLVCIYMLSTNILLVNLLVAMFG
YTVGTVQENNDQVWKFQRYFLVQEYCSRLNIPFPFIVFAYFYMVVKKCFKCCCKEKNMESSVCCFKNEDN
ETLAWEGVMKENYLVKINTKANDTSEEMRHRFRQLDTKLNDLKGLLKEIANKIK"
                        , "MSLQGTRLSMRSRRNCTLGSTRTLYSSASRSTDVSYSESDLVNFIQANFKKRECVFFTKDSKATENVCKC
GYAQSQHIEGTQINQSEKWNYKKHTKEFPTDAFGDIQFETLGKKGKYIRLSCDTDAETLYELLTQHWHLK
TPNLVISVTGGAKNFALKPRMRKIFSRLIYIAQSKGAWILTGGTHYGLMKYIGEVVRDNTISRSSEENVV
AIGIVAWGMVSNRDALIRNCDVEGYFSAQYIMDDFKRDPLYILDNNHTHLLLVDSGCHGHPTVEAKLRNQ
LEKYISERTIQDSNYGGKIPIVCFAQGGGKETLKAINTSIKSKIPCVVVEGSGQIADVIASLVEMEDALT
SSIIKEKLVRFLPRTVSRLPEEETESWIRWLKEILESSHLLTVIKMEEAGDEIVSNAISYALYKAFSTND
QDKDNWNGQLKLLLEWNQLDLANDEIFTNDRRWESADLQEVMFTALIKDRPKFVRLFLENGLNLRKFLTS
DVLTELFSNHFSSLVYQNLQIAKNSYNDALLTFVWKLVANFRRGFRKEDRNSKDEMDVELHDVSPITRHP
LQALFIWAILQNKKELSKVIWEQTKGCTLAALGASKLLKTLAKVKNDINAAGESEELANEYETRAVELFT
ECYSSDEDLAEQLLVYSCEAWGGSNCLELAVEATDQHFIAQPGVQNFLSKQWYGEISRDTKNWKIILCLF
IIPLAGCGFISFRKKPLDKHRKLLWSYVAFFTSPFVVFSWNVVFYIAFLLLFAYVLLMDFHSVPHPPELV
LYALVFVLFCDEVRQWYVNGVSYFTDLWNVMDTLGLFYFIAGIVFRLHSSHKTSLYSGRVIFCLDYIIFT
LRLVHIFTVSRNLGPKIIMLQRMLIDVFFFLFLFAVWMVAFGVARQGILRQNEHRWRWIFRSVIYEPYLA
MFGQVPSDVDGTTYDFSHCTFTGNESKPLCVELDEHNLPRFPEWITIPLVCIYMLSTNILLVNLLVAMFG
YTVGTVQENNDQVWKFQRYFLVKEYCSRLNVPFPFVVLAYFYMVVKKCFGCCCQDRGVESSACCFKNEDN
MTLAWEGVMKENYLVKINTKANDTSEEMRHRFRQLDTKLNDLKGLLKEIANKIK
")  



# Output the test results

test_results <- test_new_data(new_test_sequences, hmm_mammal, hmm_non_mammal)
print(test_results)

