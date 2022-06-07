# Purpose
This script will perform pairwise comparisons of all sequences in a fasta file to calculate dN/dS, Nd (also called Dn), and Sd (also called Ds). dN/dS is a common statistic for measuring the stregth of direction on protein sequences.

# Quick start
Download this repository. Install python3 along with the Biopython and numpy modules.

The general usage for this code is:

`python3 dnds.py <fasta> <pair/sing>`

## Sequence pair analysis mode

In this mode, the program will loop through all pairs of sequences in the fasta file and calculate dN/dS for each pair.

Apply script to pair of sequences in a fasta file
`python3 dnds.py data/examplePair.fasta pair`

Apply script to all pairs of sequences in a multifasta file
`python3 dnds.py data/exampleMulti.fasta pair`

To save output to a file:
`python3 dnds.py data/examplePair.fasta pair > exampleOutput.txt`

Apply script to a folder of fasta files
```
for FILE in data/*; do python3 dnds.py $FILE pair >> output.txt; done
```

## Single sequence analysis mode

In this mode, the program will loop through each sequence in the fasta file, individually, and calculate N and S for each.

Apply script to pair of sequences in a fasta file
`python3 dnds.py data/examplePair.fasta sing`

Apply script to all pairs of sequences in a multifasta file
`python3 dnds.py data/exampleMulti.fasta sing`

To save script output to a file:
`python3 dnds.py data/examplePair.fasta sing > exampleOutput.txt`

Apply script to a folder of fasta files
```
for FILE in data/*; do python3 dnds.py $FILE sing >> output.txt; done
```

# Understanding the input

## Sequence pair analysis mode
This script takes a single fasta-formatted file that MUST meet all of these criteria:
* There must be at least two nucleotide sequences
* All sequences must be codon-aligned
* All sequences must be the same length
* All sequences must have no gap characters
* All sequences must have no ambiguous characters
* All sequences must have no stop characters
* All sequences must not be super short (<20 bp) otherwise the dN/dS calculation could be undefined and return weird errors

See the data folder for some simplified example inputs.

Don't know how to generate input files that meet these criteria? See below

## Single sequence analysis mode

# Ideas for how to get input files
* Apply orthofinder to a folder containing protein sequences from two or more species
* Align sequences within each orthogroup of interest (usually the single-copy orthogroups) using MAFFT or a similar alignment tool
* Use seqkit grep to build fasta files containing coding sequences for each orthogroup
* Use pal2nal to convert amino acid alignments to codon-based alignments

# Understanding the output

## Sequence pair analysis mode
The script outputs tab-separated text with 11 columns to standard output.

reference - the name of the "reference" sequence in the dN/dS calculation (it shouldn't matter too much which sequence is treated as the reference)

sample - the name of the "sample" sequence in the dN/dS calculation 

N - number of nonsynonymous sites in the reference sequence

S - number of synonymous sites in the reference sequence

Nd - number of nonsynonymous differences between reference and sample

Sd - number of synonymous differences between reference and sample

pN - proportion of nonsynonymous sites that differ between reference and sample

pS - proportion of synonymous sites that differ between reference and sample

dN - the number of nonsynonymous differences per nonsynonymous site (i.e. the nonsynonymous substitution rate)

dS - the number of synonymous differences per synonymous site (i.e. the synonymous substitution rate)

dNdS - dN/dS (the ratio of the nonsynonymous and synonymous substitution rates)

Notes on output:
* N + S gives the total length of the reference sequence
* Nd + Sd gives the total number of differences between reference and sample
* pN and pS do not refer to the numbers of polymorphisms segregating in a population (although the notation is the same)

## Single sequence analysis mode
The script outputs a tab-separated text file with four columns

sequence - the name of the sequence

N - number of nonsynonymous sites in sequence

S - number of synonymous sites in sequence

L - length of sequence

# References
This script employs the Nei and Gojobori method for calculating dN/dS. See here for the details:
Nei M, Gojobori T (1986) Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. Molecular Biology and Evolution 3: 418–426
