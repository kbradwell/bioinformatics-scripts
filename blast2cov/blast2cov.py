#! /usr/bin/python

'''
Script to indicate whether genes are covered by a threshold level of sequencing reads across their entire length: provides one type of evidence for genomic presence of a gene. 
Input: 
BLAST file of genes of interest vs. reads
Output:
List of genes that are legitimately "FOUND" based on the reads and the given criteria (currently set to genes where at least 4 BLAST-aligned reads cover over 60% of the gene length). 
'''

import string
from sys import argv
import re
import os

# take input files from command line
script, blastOut, outGenes = argv

# open and read the BLAST file of genes vs. sequencing reads
allBLAST = open(blastOut)
blastLines = allBLAST.readlines()

# open the outfile to write to
of=open(outGenes,'w') 

# create dicts for the genes (queries in BLAST file)
queryDict={}
queryLenDict={}

# read the BLAST file
for line in blastLines:
	lineStr = line.strip()
	lineSplit = lineStr.split()
	genename = lineSplit[0]
	readname = lineSplit[1]
	qgene=genename.strip()
	# save the start and end coords of where the 
	gsStr = lineSplit[2].strip()
	geStr = lineSplit[3].strip()
	gs=int(gsStr)
	ge=int(geStr)
	# output from BLAST+ had the 10th column specifying the query gene length
	glen = lineSplit[9].strip()
	# save the gene name and its length
	if qgene not in queryLenDict:
		queryLenDict[qgene]=int(glen)
	# use a counter to indicate how many times each position in the query gene is covered by a read
	if qgene in queryDict:
		for pos in range(gs,ge+1):
			if pos in queryDict[qgene]: 
				queryDict[qgene][pos]+=1
			else:
				queryDict[qgene][pos]=1
	else:
		queryDict[qgene]={}
		for pos in range(gs,ge+1):
                	if pos in queryDict[qgene]:
                        	queryDict[qgene][pos]+=1
                	else:
                        	queryDict[qgene][pos]=1

# loop through each query gene
for g in queryLenDict:
	posCounter=0
	#print g
	for position in range(1,queryLenDict[g]):
		if position in queryDict[g]:
			# for each position in the gene check whether at least 4 reads cover it
			if queryDict[g][position]>3:
				posCounter+=1
				#print position,"\t", queryDict[g][position]
			else:
				pass
				#print position, "\t", queryDict[g][position]
		else:
			queryDict[g][position]=0
			#print position, "\t",queryDict[g][position]
	# check how many positions in the gene are covered by at least 4 reads
	percOver3Reads = ((float(posCounter)/float(queryLenDict[g]))*100)
	print g, "\t", percOver3Reads
	# save gene name if over 60% of the query gene is covered by at least 4 reads
	if (percOver3Reads>60):
		of.write(g)
		of.write("\n")
	
allBLAST.close()	
of.close()
