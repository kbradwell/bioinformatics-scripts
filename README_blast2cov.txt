PURPOSE
Script to indicate whether genes are covered by a threshold level of sequencing reads across their entire length: provides one type of evidence for genomic presence of a gene. 

USAGE
script <tBLASTn results file> <outfile for gene list>

INPUT
BLAST file of genes of interest vs. reads

OUTPUT
List of genes that are legitimately "FOUND" based on the reads and the given criteria (currently set to genes where at least 4 BLAST-aligned reads cover over 60% of the gene length). 
