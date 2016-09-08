#! /usr/bin/python

'''
Script to obtain percent identity of each pairwise combination of sequences in an aligned fasta file
Input: 
Sequence alignment
Output:
Percent identities for each pairwise sequence comparison
'''

import sys
from Bio import AlignIO
# alignment file is first command line argument
file2read=sys.argv[1]
align = AlignIO.read(file2read, "fasta")

# get the number of sequences in the alignment
numSeqs=len(align)

# get the length of the alignment
A=list(align[0])

for i in range(0, numSeqs):
	for j in range(0, numSeqs):
		if i != j:
			count=0
			gaps=0
			A1=list(align[i])
			A2=list(align[j])
			for n in range(0, len(A)):
				f=A1[n]
				s=A2[n]
				if f.upper()==s.upper():
					if f!="-":
						count=count+1
					else:
						gaps=gaps+1
			#print count
			#print float((len(A)-gaps))
			print align[i].id
			print align[j].id
			print "Percent identity", 100*(count/float((len(A)-gaps)))
