gen_mutations <- function(seq, N) {
  sealKey() # See: http://steipe.biochemistry.utoronto.ca/abc/index.php/BCH441_Code_submisson_instructions
  stats <- c()
  stats <- cbind(stats, c(0, 0, 0))
  rownames(stats) <- c("silent", "missense", "nonsense")
  colnames(stats) <- c("occurrences")
  # Actual function
  for (i in 1:N) {
    original_seq <- Biostrings::DNAString(seq)
    aa_seq <- Biostrings::translate(original_seq, no.init.codon = TRUE)

    mut_seq <- Biostrings::DNAString(seq)
    mut_index <- sample(1:length(original_seq), 1, replace = TRUE)
    possible_mutations <- Biostrings::DNA_BASES
    possible_mutations <- possible_mutations[possible_mutations != as.character(unlist(original_seq[mut_index]))]
    mut_seq <- Biostrings::replaceLetterAt(mut_seq, mut_index, sample(possible_mutations, 1, replace = TRUE))
    mut_aa <- Biostrings::translate(mut_seq, no.init.codon = TRUE)


    term_aa <- regexpr(pattern = "\\*", aa_seq)
    term_mut_aa <- as.integer(regexpr(pattern = "\\*", mut_aa))
    if ((term_aa == -1 && term_mut_aa != -1) || (term_mut_aa != -1 && term_mut_aa < term_aa)) {
      stats["nonsense", "occurrences"] <- 1 + stats["nonsense", "occurrences"]
    } else if (mut_aa == aa_seq) {
      stats["silent", "occurrences"] <- 1 + stats["silent", "occurrences"]
    } else {
      stats["missense", "occurrences"] <- 1 + stats["missense", "occurrences"]
    }
  }
  sealKey()
  return(stats)
}

gen_mutations("ATGATGATGATGATGATG", 1000)
gen_mutations("CCCCCCCCCCCCCCCCCC", 500)
gen_mutations("TATTACTATTACTATTAC", 500)
gen_mutations("TGGTGGTGGTGGTGGTGGTGGTGG", 500)
gen_mutations("TGTTGTTGTTGTTGTTGTTGTTGT", 500)
gen_mutations("TGTTGTTGTTGTTGTTGTTGTTGA", 500)


myFA <-             readFASTA("data/RAB39B_HSa_coding.fa")
myFA <- rbind(myFA, readFASTA("data/PTPN5_HSa_coding.fa"))
myFA <- rbind(myFA, readFASTA("data/PTPN11_HSa_coding.fa"))
myFA <- rbind(myFA, readFASTA("data/KRAS_HSa_coding.fa"))
rownames(myFA)<-c("RAB39B", "PTPN5", "PTPN11", "KRAS") # Assign row names

gen_mutations(myFA["RAB39B", 2], 10000)
gen_mutations(myFA["PTPN5", 2], 10000)
gen_mutations(myFA["PTPN11", 2], 10000)
gen_mutations(myFA["KRAS", 2], 10000)
