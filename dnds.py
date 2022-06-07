from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
import sys
import numpy as np

# input for overall script, just fasta file of codon-aligned sequences
fastaInput = sys.argv[1]
analysisMode = sys.argv[2]

# DNA alphabet, could be changed but probably not
alphabet = ["A", "T", "G", "C"]

#############
# FUNCTIONS #
#############

### Determine expected number of nonsynonymous sites in a codon
# input: codon
# output: expected number of nonsynonymous sites for codon
def nonsynonymousSites(codon):
  codon = MutableSeq(codon)
  first = codon[0]
  second = codon[1]
  third = codon[2]
  mutatedCodons = list()

  # Generate a list of codons mutated at first position
  alphabetMinusOne = [x for x in alphabet if x not in first]
  # copy codon sequence, need [:] to avoid changing both original and copy
  mutatedCodon = codon[:]
  for letter in alphabetMinusOne:
    mutatedCodon[0] = letter
    mutatedCodons.append(mutatedCodon[:])
    #print(mutatedCodon)

  alphabetMinusOne = [x for x in alphabet if x not in second]
  mutatedCodon = codon[:]
  for letter in alphabetMinusOne:
    mutatedCodon[1] = letter
    mutatedCodons.append(mutatedCodon[:])
    #print(mutatedCodon)

  alphabetMinusOne = [x for x in alphabet if x not in third]
  mutatedCodon = codon[:]
  for letter in alphabetMinusOne:
    mutatedCodon[2] = letter
    mutatedCodons.append(mutatedCodon[:])
    #print(mutatedCodon)

  # Translate seq to identify nonsynonymous mutations. For every nonsynonymous mutation, add 1/3 to number of nonsynonymous sites
  n = 0
  for mutatedCodon in mutatedCodons:
    if mutatedCodon.translate() != codon.translate():
      n = n + (1/3)
  return(n)

### Count substitutions in a step-wise mutational pathway
# input: list of codons from a mutational pathway
# output: number of non-synonymous and synonymous differences in pathway
def countSubs(pathway):
  aaList = [x.translate() for x in pathway]
  nd = 0
  sd = 0
  for i in range(len(aaList) - 1):
    if aaList[i] != aaList[i+1]:
      nd = nd + 1
    else:
      sd = sd + 1
  return([nd, sd])

### Determine distance between reference and sample codon
# input: two codons of same length
# output: distance between two codons (i.e. number of pairwise differences)
def codonDistance(codon1, codon2):
  # make sequences mutable
  codon1 = MutableSeq(codon1)
  codon2 = MutableSeq(codon2)

  # vectors to store outputs
  distance = 0
  nd = 0
  sd = 0
  mismatchPositions = list() # mismatches between codons

  #calculate differences between two codons
  for i in range(len(codon1)):
    letter1 = codon1[i]
    letter2 = codon2[i]
    if letter1 != letter2:
      distance = distance + 1
      mismatchPositions.append(i)

  # If the codons are identical, then there are no differences
  if distance == 0:
    return([nd, sd])
  # If codons differ by one basepair, their is either one synonymous or nonsynonymous site
  if distance == 1:
    if codon1.translate() == codon2.translate():
      sd = sd + 1
      return([nd, sd])
    else:
      nd = nd + 1
      return([nd, sd])

  # If codons differ by two basepairs, there are two possible step-wise mutational pathways
  if distance == 2:
    # construct pathway one
    pathway = list()
    interCodon = codon1[:] # intermediate codon
    pathway.append(interCodon[:])
    interCodon[mismatchPositions[0]] = codon2[mismatchPositions[0]]
    pathway.append(interCodon[:])
    interCodon[mismatchPositions[1]] = codon2[mismatchPositions[1]]
    pathway.append(interCodon[:])
    pathway1Count = countSubs(pathway) # count subs
    #print(pathway)

    # construct pathway two
    pathway = list()
    interCodon = codon1[:]
    pathway.append(interCodon[:])
    interCodon[mismatchPositions[1]] = codon2[mismatchPositions[1]]
    pathway.append(interCodon[:])
    interCodon[mismatchPositions[0]] = codon2[mismatchPositions[0]]
    pathway.append(interCodon[:])
    pathway2Count = countSubs(pathway) # count substitutions
    #print(pathway)

    # assume pathways equally probable, so average result
    sum_list = [(a + b)/2 for a, b in zip(pathway1Count, pathway2Count)]
    return(sum_list)

  if distance == 3:
    mutPos = [[0,1,2],[0,2,1],[1,2,0],[1,0,2],[2,1,0],[2,0,1]]
    subCounts = list() # output counts of substitutions per pathway
    for i in mutPos:
      pathway = list()
      interCodon = codon1[:] # intermediate codon
      pathway.append(interCodon[:])
      for j in i:
        interCodon[j] = codon2[j]
        pathway.append(interCodon[:])
      subCounts.append(countSubs(pathway)) # count subs
      #print(pathway)
      #print(countSubs(pathway))
    #print(subCounts)
    #print(list(np.sum(subCounts, axis=0)/6))
    #sum_list = [sum(x) for x in zip(*subCounts)]
    # assume pathways are equally probable, so average result
    return(list(np.sum(subCounts, axis=0)/6))

##############
# LOAD FASTA #
##############

# Load sequences
records = list(SeqIO.parse(fastaInput, "fasta"))

# Number of sequences
recNum = len(records)

##########################
# SEQUENCE PAIR ANALYSIS #
##########################

if analysisMode == "pair":
	# print header of output
	print("\t".join(["reference", "sample", "N", "S", "Nd", "Sd", "pN", "pS", "dN", "dS", "dNdS"]))

	# Loop over all pairs of sequences
	for j in range((recNum - 1)):
		for k in range(j+1, recNum):
			refRecord = records[j]
			samRecord = records[k]
			refSeq = refRecord.seq
			samSeq = samRecord.seq
			samName = samRecord.name
			refName = refRecord.name

    			# count number of codons in sequence
			codonCount = len(refSeq)/3

    			# counts of nonsynonymous sites
			Nref = list()
			Nsam = list()

    			# loop over codons in reference sequence
			for i in range(int(codonCount)):
				codonRefSeq = refSeq[i*3:i*3+3]
				codonSamSeq = samSeq[i*3:i*3+3]
				Nref.append(nonsynonymousSites(codonRefSeq))
				Nsam.append(nonsynonymousSites(codonSamSeq))

    			# sum over nonsynonymous sites of individual codons to get total
			Nref = sum(Nref)
			Nsam = sum(Nsam)

    			# average the number of nonsynonymous sites in both sequences
			N = (Nref + Nsam)/2

    			# calculate number of synonymous sites: total sites - nonsynonymous sites
			Sref = 3*codonCount - Nref
			Ssam = 3*codonCount - Nsam

			# average synonymous sites
			S = (Sref + Ssam)/2

    			# calculate number of differences between sequence pair
			NdSd = list()
			for i in range(int(codonCount)):
				codonRefSeq = refSeq[i*3:i*3+3]
				codonSamSeq = samSeq[i*3:i*3+3]
				NdSd.append(codonDistance(codonRefSeq, codonSamSeq))

    			# Sum differences across codons to get total differences
			NdSd = np.sum(NdSd, axis=0)
			Nd = NdSd[0]
			Sd = NdSd[1]

    			# Calculate proportion of differences
			pN = Nd/N
			pS = Sd/S

    			# Calculate substitution rates
			# correct for multiple substitutions per site using Jukes and Cantor 1969 formula
			dN = (-3/4)*np.log(1-(4/3)*pN)
			dS = (-3/4)*np.log(1-(4/3)*pS)

    			# output final result
			print("\t".join([refName, samName, str(N), str(S), str(Nd), str(Sd), str(pN), str(pS), str(dN), str(dS), str(dN/dS)]))

############################
# SINGLE SEQUENCE ANALYSIS #
############################

if analysisMode == "sing":
	print("\t".join(["sequence", "N", "S", "L"]))
	# Loop over each sequence, counting nonsynonymous sites
	for j in range(recNum):
		oneRecord = records[j]
		oneSeq = oneRecord.seq
		oneName = oneRecord.name

		# Get sequence length
		L = len(oneSeq)

		# count number of codons in sequence
		codonCount = L/3

    		# counts of nonsynonymous sites per codon
		N = list()

    		# loop over codons, counting nonsynonymous sites
		for i in range(int(codonCount)):
			codonSeq = oneSeq[i*3:i*3+3]
			N.append(nonsynonymousSites(codonSeq))

		# sum over nonsynonymous sites of individual codons to get total
		N = sum(N)

		# calculate number of synonymous sites: total sites - nonsynonymous sites
		S = 3*codonCount - N

		# output final result
		print("\t".join([oneName, str(N), str(S), str(L)]))

# Test countSubs function
#print(countSubs([Seq("CCG"), Seq("CTG"), Seq("CTA")]))
#print(countSubs([Seq("CCG"), Seq("CCA"), Seq("CTA")]))
#print(countSubs([Seq("AAA"), Seq("TAA"), Seq("TTA"), Seq("TTG")]))

# Test codonDistance
#print(codonDistance(Seq("AAA"), Seq("AAA"))) # no difference
#print(codonDistance(Seq("AAA"), Seq("TAA"))) # one nonsynonymous difference
#print(codonDistance(Seq("AAA"), Seq("AAG"))) # one synonymous difference
#print(codonDistance(Seq("TTT"), Seq("GTA"))) # two differences
#print(codonDistance(Seq("AAA"), Seq("TTG"))) # three differences
