# This script is designed to calculate a normalized index of hydrophobic clustering,
# similar to Jason's Python script. It adjusts average run length by the expected
# run length from amino acid composition.

####################################################################
# IMPORTANT: CURRENTLY DOES NOT ACCOUNT FOR UNCERTAIN AMINO ACIDS! #
####################################################################

# Required packages.
library(stringr)

# Define the hydrophobic residues here. FILMV are the most hydrophobic and are listed
# below as the default, but the list may be edited for other combinations.
hydrophobic.residues <- c("F", "I", "L", "M", "V") # Five most hydrophobic residues.

# Function for counting the run lengths of hydrophobic residues in an amino acid sequence.
# Requires an amino acid sequence as an input, and returns a list of hydrophobic residues.
# Needs to reference the "hydrophobic.residues" list above to determine which amino acids
# are hydrophobic. Also returns mean run length.
# NOTE: Requires the "stringr" package!
count.hydrophobic.run.lengths <- function(aa.sequence){
  counts.vector <- c(rep(0, str_length(aa.sequence)))
  aa.vector <- str_split(aa.sequence, "")
  aa.length <- str_length(aa.sequence)
  hydro.run.length <- 0
  for (i in 1:aa.length) {
    if (aa.vector[[1]][i] %in% hydrophobic.residues){
      hydro.run.length <- hydro.run.length + 1
    } else {
      if (hydro.run.length > 0) {
        counts.vector[hydro.run.length] <- counts.vector[hydro.run.length] + 1
        hydro.run.length <- 0
      }
    }
  }
  if (hydro.run.length > 0){
    counts.vector[hydro.run.length] <- counts.vector[hydro.run.length] + 1
    hydro.run.length <- 0
  }
  return(counts.vector)
}

# Function for calculating the normalized index of dispersion.
normalized.hydro.run.length <- function(aa.seq){
  run.lengths <- count.hydrophobic.run.lengths(aa.sequence = aa.seq)
  aa.vector <- str_split(aa.seq, "")[[1]]
  aa.table <- table(aa.vector)
  total.hydro <- sum(aa.table[names(aa.table) %in% hydrophobic.residues])
  mean.run.length <- total.hydro / sum(run.lengths)
  hydro.percent <- total.hydro / str_length(aa.seq)
  normalization <- 1 / (1 - hydro.percent)
  return(mean.run.length / normalization)
}