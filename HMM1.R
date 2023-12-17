

# Install BiocManager if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Biostrings
BiocManager::install("Biostrings")
install.packages("ggplot2")

install.packages("HMM")
install.packages("tidyr")
install.packages("pheatmap")


# Required libraries
library(Biostrings)
library(HMM)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(dplyr)

read_fasta <- function(file_path) {
  valid_amino_acids <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
                         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  
  fasta_data <- readLines(file_path)
  sequences <- c()
  current_seq <- ""
  
  for (line in fasta_data) {
    if (startsWith(line, ">")) {
      if (current_seq != "") {
        sequences <- c(sequences, current_seq)
        current_seq <- ""
      }
    } else {
      current_seq <- paste0(current_seq, line)
    }
  }
  
  if (current_seq != "") {
    sequences <- c(sequences, current_seq)
  }
  
  sequences <- lapply(sequences, function(seq) {
    paste(sapply(strsplit(seq, "")[[1]], function(x) 
      if (x %in% valid_amino_acids) x else NA_character_), collapse = "")
  })
  
  return(sequences)
}

classify_amino_acid <- function(amino_acid) {
  polar <- c('S', 'T', 'C', 'Y', 'N', 'Q')
  non_polar <- c('A', 'V', 'L', 'I', 'P', 'W', 'F', 'M', 'G')
  negative <- c('D', 'E')
  positive <- c('R', 'K', 'H')
  
  if (amino_acid %in% polar) {
    return('polar')
  } else if (amino_acid %in% non_polar) {
    return('non_polar')
  } else if (amino_acid %in% negative) {
    return('negative')
  } else if (amino_acid %in% positive) {
    return('positive')
  } else {
    return(NA)  # for amino acids not in the list
  }
}

# Function to encode sequences based on amino acid states
encode_sequence <- function(sequence) {
  states <- sapply(strsplit(sequence, "")[[1]], classify_amino_acid)
  amino_acids <- strsplit(sequence, "")[[1]]
  encoded_pairs <- paste(amino_acids, states, sep=":")
  return(encoded_pairs)
}


states <- c("polar", "non_polar", "negative", "positive")

train_hmm <- function(encoded_sequences) {
  valid_amino_acids <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
                         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  states <- c("polar", "non_polar", "negative", "positive")
  
  start_states <- sapply(encoded_sequences, function(x) strsplit(x[1], ":")[[1]][2])
  startProbs <- table(factor(start_states, levels = states)) / length(start_states)
  startProbs <- as.numeric(startProbs)  # Convert to numeric
  
  state_indices <- setNames(1:length(states), states)
  
  
  transitionProbs <- matrix(0, nrow = length(states), ncol = length(states))
  colnames(transitionProbs) <- rownames(transitionProbs) <- states
  
  for (seq in encoded_sequences) {
    for (i in 1:(length(seq) - 1)) {
      current_state <- state_indices[strsplit(seq[i], ":")[[1]][2]]
      next_state <- state_indices[strsplit(seq[i + 1], ":")[[1]][2]]
      transitionProbs[current_state, next_state] <- transitionProbs[current_state, next_state] + 1
    }
  }
  transitionProbs <- sweep(transitionProbs, 1, rowSums(transitionProbs), FUN="/")
  
  emissionProbs <- matrix(0, nrow = length(states), ncol = length(valid_amino_acids))
  rownames(emissionProbs) <- states
  colnames(emissionProbs) <- valid_amino_acids
  
  for (seq in encoded_sequences) {
    for (state_pair in seq) {
      state <- state_indices[strsplit(state_pair, ":")[[1]][2]]
      amino_acid <- strsplit(state_pair, ":")[[1]][1]
      emissionProbs[state, amino_acid] <- emissionProbs[state, amino_acid] + 1
    }
  }
  emissionProbs <- sweep(emissionProbs, 1, rowSums(emissionProbs), FUN="/")
  
  hmm <- initHMM(states, valid_amino_acids, startProbs, transitionProbs, emissionProbs)
  
  
  #attr(hmm, "symbols") <- valid_amino_acids
  hmm$symbols <- valid_amino_acids
  
  hmm$transitionProbs <- transitionProbs
  hmm$emissionProbs <- emissionProbs
  
  return(hmm)
}

mammal_sequences <- read_fasta("mammal1_train.fasta")
non_mammal_sequences <- read_fasta("bird1_train.fasta")

mammal_encoded_sequences <- lapply(mammal_sequences, encode_sequence)
non_mammal_encoded_sequences <- lapply(non_mammal_sequences, encode_sequence)

hmm_mammal <- train_hmm(mammal_encoded_sequences)
hmm_non_mammal <- train_hmm(non_mammal_encoded_sequences)

predict_with_hmm <- function(encoded_sequence, hmm_model) {
  amino_acid_sequence <- sapply(strsplit(encoded_sequence, ":"), function(pair) pair[1])
  
  
  valid_observation <- amino_acid_sequence[amino_acid_sequence %in% hmm_model$symbols]
  if (length(valid_observation) == 0) {
    print("No valid observations found. Returning NA.")
    return(NA)
  }
  viterbi_path <- viterbi(hmm_model, observation = valid_observation)
  
 
  if (is.character(viterbi_path)) {
    return(viterbi_path)
  } else {
   
    warning("Unexpected return type from viterbi function")
    return(NA)
  }
}

test_new_data <- function(new_sequences, hmm_mammal, hmm_non_mammal) {
  encoded_sequences <- lapply(new_sequences, encode_sequence)
  
  
  predictions_mammal <- lapply(encoded_sequences, predict_with_hmm, hmm_model = hmm_mammal)
  predictions_non_mammal <- lapply(encoded_sequences, predict_with_hmm, hmm_model = hmm_non_mammal)
  
  list(mammal = predictions_mammal, non_mammal = predictions_non_mammal)
}

new_test_sequences <- read_fasta("bird_testing.fasta")

print(hmm_mammal$transitionProbs)
print(hmm_non_mammal$transitionProbs)

# Testing the HMM models on new data
test_results1 <- test_new_data(new_test_sequences, hmm_mammal, hmm_non_mammal)
print(test_results1)

guess_file_type <- function(test_) {
  lines <- readLines(file_path, n = 5, warn = FALSE)
  
  if (any(grepl("^>", lines))) {
    return("FASTA")
  } else if (any(grepl(",", lines))) {
    return("CSV")
  } else {
    return("Unknown")
  }
}

write_output_to_csv <- function(test_results, file_path) {
  output_data <- data.frame(Model = character(), SequenceID = integer(), States = character(), stringsAsFactors = FALSE)
  
  for (model in names(test_results)) {
    sequences <- test_results[[model]]
    for (i in seq_along(sequences)) {
      states_str <- paste(sequences[[i]], collapse = ", ")
      output_data <- rbind(output_data, data.frame(Model = model, SequenceID = i, States = states_str))
    }
  }
  
  write.csv(output_data, file_path, row.names = FALSE)
}

# Save to a CSV file
write_output_to_csv(test_results1, "test_results1_output.csv")


# Read the CSV file
data <- read.csv("test_results1_output.csv", stringsAsFactors = FALSE)


str(data)


data_long <- data %>%
  separate_rows(States, sep = ",\\s*")
head(data_long)

data_long_percent <- data_long %>%
  count(States) %>%
  mutate(FreqPercent = n / sum(n) * 100)


ggplot(data_long_percent, aes(x = States, y = FreqPercent, fill = States)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +  # Using a color scale suitable for discrete variables
  theme_minimal() +
  labs(title = "Combined State Frequency Distribution Across All Models (Percentages)", 
       x = "State", 
       y = "Frequency (%)")

ggplot(data_long, aes(x = States, fill = Model)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "State Frequency Distribution", x = "State", y = "Frequency") +
  facet_wrap(~ Model)

data$States <- strsplit(data$States, ",\\s*")

calculate_transition_matrix <- function(sequences) {
  # Initialize matrix with state names and zeros
  unique_states <- unique(unlist(sequences))
  transition_matrix <- matrix(0, nrow = length(unique_states), ncol = length(unique_states),
                              dimnames = list(unique_states, unique_states))
  
  
  for (seq in sequences) {
    for (i in 1:(length(seq) - 1)) {
      transition_matrix[seq[i], seq[i + 1]] <- transition_matrix[seq[i], seq[i + 1]] + 1
    }
  }
  
  transition_matrix <- sweep(transition_matrix, 1, rowSums(transition_matrix), FUN="/")
  return(transition_matrix)
}

transition_matrix_mammal <- calculate_transition_matrix(data[data$Model == "mammal", "States"])
transition_matrix_non_mammal <- calculate_transition_matrix(data[data$Model == "non_mammal", "States"])

pheatmap(transition_matrix_mammal, main = "Transition Probabilities - Mammal")
pheatmap(transition_matrix_non_mammal, main = "Transition Probabilities - Non-Mammal")

identify_extremes <- function(transition_matrix) {
  max_transition <- which(transition_matrix == max(transition_matrix), arr.ind = TRUE)
  min_transition <- which(transition_matrix == min(transition_matrix), arr.ind = TRUE)
  list(max_transition = max_transition, min_transition = min_transition)
}

extremes_mammal <- identify_extremes(transition_matrix_mammal)
extremes_non_mammal <- identify_extremes(transition_matrix_non_mammal)

print(extremes_mammal)
print(extremes_non_mammal)

mammal_sequences <- read_fasta("mammal1_train.fasta")

first_mammal_encoded_sequence <- encode_sequence(mammal_sequences[[1]])


calculate_emission_probabilities_single_sequence <- function(encoded_sequence) {
  valid_amino_acids <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
                         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  states <- c("polar", "non_polar", "negative", "positive")
  
  
  emission_matrix <- matrix(0, nrow = length(valid_amino_acids), ncol = length(states),
                            dimnames = list(valid_amino_acids, states))
  
  for (state_pair in encoded_sequence) {
    components <- strsplit(state_pair, ":")[[1]]
    state <- components[2]  # State is the second element
    amino_acid <- components[1]  # Amino acid is the first element
    if (amino_acid %in% valid_amino_acids && state %in% states) {
      emission_matrix[amino_acid, state] <- emission_matrix[amino_acid, state] + 1
    }
  }
  
  emission_matrix <- sweep(emission_matrix, 2, colSums(emission_matrix), FUN="/")
  return(emission_matrix)
}

emission_probabilities_first_sequence <- calculate_emission_probabilities_single_sequence(first_mammal_encoded_sequence)
pheatmap(emission_probabilities_first_sequence, main = "Emission Probabilities - First Mammalian Sequence")

calculate_aa_transition_matrix <- function(sequences) {
  valid_amino_acids <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
                         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  
  
  transition_matrix <- matrix(0, nrow = length(valid_amino_acids), 
                              ncol = length(valid_amino_acids), 
                              dimnames = list(valid_amino_acids, valid_amino_acids))
  
  for (sequence in sequences) {
    amino_acids <- unlist(strsplit(sequence, ""))
    for (i in 1:(length(amino_acids) - 1)) {
      if (amino_acids[i] %in% valid_amino_acids && amino_acids[i + 1] %in% valid_amino_acids) {
        transition_matrix[amino_acids[i], amino_acids[i + 1]] <- transition_matrix[amino_acids[i], amino_acids[i + 1]] + 1
      }
    }
  }
  
  transition_matrix <- sweep(transition_matrix, 1, rowSums(transition_matrix), FUN = "/")
  return(transition_matrix)
}

mammal_sequences <- read_fasta("mammal1_train.fasta")
aa_transition_matrix <- calculate_aa_transition_matrix(mammal_sequences)

pheatmap(aa_transition_matrix, main = "Amino Acid Transition Probabilities - Mammalian Sequences")
