myFA <-             readFASTA("data/RAB39B_HSa_coding.fa")
myFA <- rbind(myFA, readFASTA("data/PTPN5_HSa_coding.fa"))
myFA <- rbind(myFA, readFASTA("data/PTPN11_HSa_coding.fa"))
myFA <- rbind(myFA, readFASTA("data/KRAS_HSa_coding.fa"))
rownames(myFA)<-c("RAB39B", "PTPN5", "PTPN11", "KRAS") # Assign row names

gen_mutations <- function(seq, N) {
  stats <- c()
  stats <- cbind(stats, c(0, 0, 0))
  rownames(stats) <- c("silent", "missense", "nonsense")
  colnames(stats) <- c("occurrences")
  # Actual function
  for (i in 1:217) {
    # select index for mutation
    working_seq <- Biostrings::DNAString(seq)
    aa_seq <- Biostrings::translate(working_seq, no.init.codon = TRUE)
    mut_action <- sample(c("ins", "del", "sub"), 1, TRUE)
    mut_seq <- Biostrings::DNAString(seq)
    if (mut_action == "sub") {
      mut_index <- sample(1:length(working_seq), 1, replace = TRUE)
      possible_mutations <- Biostrings::DNA_BASES
      possible_mutations <- possible_mutations[possible_mutations != as.character(unlist(working_seq[mut_index]))]
      mut_change <- sample(possible_mutations, 1, replace = TRUE)
      mut_seq <- Biostrings::replaceLetterAt(mut_seq, mut_index, mut_change)
    } else if (mut_action == "ins") {
      mut_index <- sample(1:length(working_seq) - 2, 1, replace = TRUE)
      possible_mutations <- Biostrings::DNA_BASES
      mut_seq <- Biostrings::DNAString(paste(substring(working_seq, 1, mut_index - 1), sample(possible_mutations, 1), substring(working_seq, mut_index), sep = ""))
    } else {
      mut_index <- sample(1:length(working_seq), 1, replace = TRUE)
      mut_seq <- mut_seq[-mut_index]
    }
    mut_seq <- Biostrings::DNAString(substring(mut_seq, 1, length(mut_seq) - (length(mut_seq) %% 3)))
    mut_aa <- Biostrings::translate(mut_seq, no.init.codon = TRUE)

    # Note: we need silent, nonsense, and missense
    mut_aa_stop <- match("*", Biostrings::as.matrix(mut_aa))
    aa_seq_stop <- match("*", Biostrings::as.matrix(aa_seq))
    if (!is.na(mut_aa_stop) & (is.na(aa_seq_stop) | mut_aa_stop < aa_seq_stop)) {
      stats["nonsense", "occurrences"] <- 1 + stats["nonsense", "occurrences"]
    } else if (mut_aa == aa_seq) {
      stats["silent", "occurrences"] <- 1 + stats["silent", "occurrences"]
    } else {
      stats["missense", "occurrences"] <- 1 + stats["missense", "occurrences"]
    }
  }
  return(stats)
}
N_test <- 1200
gen_mutations("ATGATGATGATGATGATG", N_test)
gen_mutations("CCCCCCCCCCCCCCCCCC", N_test)
gen_mutations("TATTACTATTACTATTAC", N_test)
gen_mutations("TGGTGGTGGTGGTGGTGGTGGTGG", N_test)
gen_mutations("TGTTGTTGTTGTTGTTGTTGTTGT", N_test)
