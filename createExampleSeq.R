
library(gtools)
#install.packages("gtools")
# generate a sequence that has one of every codon
# this could be useful for debugging my dn/ds script
seqMat = permutations(4, 3, v = c("A", "C", "T", "G"), repeats.allowed = T)
seqVec = as.vector(t(seqMat)) 
seqStr = paste(seqVec, collapse = "")

# Generate one seq for every codon, used to make sure I am 
for(i in 1:(length(seqVec)/3)){
  codonVec = seqVec[((i-1)*3+1):(i*3)]
  codonStr = paste(codonVec, collapse = "")
  cat(paste(">", codonStr, sep = ""))
  cat("\n")
  cat(paste(codonStr, collapse = ""))
  cat("\n")
}

