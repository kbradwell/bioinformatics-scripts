#! /usr/bin/python

'''
Script takes maf and blastab format files from lastal (LAST, http://last.cbrc.jp/) search and makes a gff of potential pseudogenes (without filtering for
called gene overlaps).
Input: 
MAF and BLASTTAB format files from lastal (query file, i.e. genome assembly sequences, must contain "contig" in query names). 
Output:
GFF file of putative pseudogenes. This can be used to subsequently filter out any coordinates overlapping called genes.
'''

import sys
from sys import argv
import string
import re
import os
import itertools
import subprocess

script, mafInfile, blasttabInfile, outGFF = argv

blastdbcmdLoc = None # add in path as double quoted string here (blasdbcmd program location e.g. "/PATH/TO/PROGRAM/blastdbcmd")
dbLoc = None # add in path as double quoted string here (nr database location e.g. "/PATH/TO/DATABASE/nr")

if blastdbcmdLoc is None or dbLoc is None: 
	sys.exit("please add full path of blastdbcmd program and database (e.g. nr) to beginning of script.")

# read in each alignment from the MAF file, and ignore any that don't contain the req gb or ref subj hit or > 40% subject aln

f=open(mafInfile,"r")

dbstrings = ("|ref|","|gb|")
firstFlag=0
arrayC=[]
for line in f:
	if (line.startswith( 'a' ) and firstFlag==0):
		firstFlag=1
		scoreStr=line
		scoreList=scoreStr.split()
		evalStr=scoreList[3]
		evalList=evalStr.split("=")
		eval=evalList[1]
	elif (line.startswith( 's' ) and "gi|" in line and any (s in line for s in dbstrings)):
		subjline=line
	elif (line.startswith( 's' ) and "contig" in line.lower()):
		queryline=line
	elif (line.startswith( 's' ) and "gi|" in line and not any (s in line for s in dbstrings)):
		subjline="NOTCORRECTDB"
	elif (line.startswith( 'a' ) and firstFlag==1):
		if "NOTCORRECTDB" not in subjline:
			# check that >40% of the subject is aligned
			subjLineList=subjline.split()
			subjLen=subjLineList[5]
			subjAlnLen=subjLineList[3]
			fracSubj=float(subjAlnLen)/float(subjLen)
			if fracSubj>0.4:
				queryLineList=queryline.split()
				# save to arrayC
				contig=queryLineList[1]
				alnStart=queryLineList[2]
				alnLen=queryLineList[3]
				alnEnd=int(alnStart)+int(alnLen)
				strand=queryLineList[4]
				queryLen=queryLineList[5]
				alnSeq=queryLineList[6]
				if "-" in strand:
					realAlnStart=int(queryLen)-(int(alnStart)+int(alnLen))
					realAlnEnd=int(queryLen)-int(alnStart)
				else:
					realAlnStart=alnStart
					realAlnEnd=alnEnd
				subjHit=subjLineList[1]
				percSubj=fracSubj*100
				coordInfo=[contig,realAlnStart,realAlnEnd,queryLen,strand,subjHit,eval,subjLen,percSubj,alnSeq]
				#print "subj line", subjline
				#print "query line", queryline
				#print "eval", eval
				#print coordInfo
				arrayC.append(coordInfo)
		scoreStr=line
		scoreList=scoreStr.split()
		evalStr=scoreList[3]
		evalList=evalStr.split("=")
		eval=evalList[1]


f.close()

#remove overlapping coords from arrayC, keeping the longest one

sArrayC = sorted(arrayC, key = lambda x: (x[0], int(x[1])))
#print "initial coords list length", len(sArrayC)
toRem=[]
possRem=[]
tied={}
for i in sArrayC:
	if i in toRem:
		continue
	#print "assessing i",i
	for j in sArrayC:
		if j in toRem:
			continue
		#print "assessing j",j
		elif i==j:
			continue
		elif i[0] != j[0]:
			continue
		else:
			x = range(int(i[1]),int(i[2]))
			y = range(int(j[1]),int(j[2]))
			xs = set(x)
			overlap=xs.intersection(y)
			#print "overlap length",len(overlap)
			if len(overlap)>0:
				iInterval = int(i[2])-int(i[1])
				jInterval = int(j[2])-int(j[1])
				if iInterval>jInterval:
					toRem.append(j)
					#print "removed",j,"kept",i
				elif jInterval>iInterval:
					toRem.append(i)
					#print "removed",i,"kept",j
				elif jInterval==iInterval:
					# choose the one with the best eval or if evals same choose to remove j
					if float(i[6])<float(j[6]):
						toRem.append(j)
					elif float(i[6])==float(j[6]):
						possRem.append(j)
						possRem.append(i)
						tiedKey=i[0]+str(i[1])
						if tiedKey not in tied:
							tied[tiedKey]=[]
							tied[tiedKey].append(i)
						else:
							tied[tiedKey].append(i)
					else:
						toRem.append(i)

longestCoords = [x for x in sArrayC if x not in toRem]
nonDupsLongestCoords = [x for x in longestCoords if x not in possRem]

for tk in tied:
	if tied[tk][0] not in toRem:
		if tied[tk][0] not in nonDupsLongestCoords:
			#print "tied key", tied[tk][0]
			nonDupsLongestCoords.append(tied[tk][0])

#print "final longest coords list length", len(nonDupsLongestCoords)

# get only the pseudogene coords
pseudoChars = ["*", "\\","/"]

toRem=[]
for c in nonDupsLongestCoords:
	alnseq = c[9]
	if any(char in alnseq for char in pseudoChars):
		pass
		#print "pseudo char", alnseq
	else:
		toRem.append(c)
		#print "not pseudo char", alnseq

pseudoCoords = [x for x in nonDupsLongestCoords if x not in toRem]

#print "only pseudos list length", len(pseudoCoords)
#for c in pseudoCoords:
#	print "pseudo",c 
#	print "\n"

# putative pseudogene coordinates are now saved - percID, descriptions of hits, and code for type of pseudogene feature need to be added

# first read the blasttab file to get the percID, and add the pseudogene feature hallmarks and blastdbcmd descriptions

frameshifts=["/","\\"]
premStop=["*"]
f2=open(blasttabInfile,"r")
for line in f2:
	lineArr=line.split()
	contig=lineArr[0]
	subj=lineArr[1]
	st=lineArr[6]
	end=lineArr[7]
	bCoords=[int(st),int(end)]
	coords=sorted(bCoords)
	for c in pseudoCoords:
		if ((c[0] in contig) and (c[5] in subj) and (abs(int(coords[0]) - int(c[1])) < 3)):
			pseudoDescArr=[]
			percID=lineArr[2]
			#print "obtained perc ID"
			c.append(percID)
			if any(char in c[9] for char in frameshifts):
				pseudoDescArr.append("F")
			if any(char in c[9] for char in premStop):
				pseudoDescArr.append("S")
			pseudoDescStr=",".join(pseudoDescArr)
			c.append(pseudoDescStr)
			id=c[5]
			parts=id.split("|")
			gi=parts[1]
			dbt='prot'
			try:
				proc = subprocess.Popen([blastdbcmdLoc,"-db", dbLoc, "-dbtype", dbt, "-entry",gi], stdout=subprocess.PIPE)
				tmp = proc.stdout.read()
				tlList=tmp.split("\n")
				descList=tlList[0].split(">")
				firstDesc=descList[1]
				firstDescParts=firstDesc.split("|")
				nonIDdesc=firstDescParts[4].strip()
				c.append(nonIDdesc)
			except:
				c.append("no description")
				continue
					
f2.close()

# all info for coords are now saved and are ready to make into gff to rem the called gene intersection coords

sortedPseudos=sorted(pseudoCoords, key = lambda x: (x[0], int(x[1])))

of=open(outGFF,"w")
# print the gff info to the GFF file 
for c in sortedPseudos:
	# parts of the GFF 
	# c[0]=Q_ID,c[1]=Q_BEG, c[2]=Q_END, c[3]=Q_LEN, c[4]=STRAND, c[5]=S_ID, c[6]=E_VAL, c[7]S_LEN, c[8]=PERC_S_LEN, c[10]=IDENT, c[11]=PSEUDO_FEATURE(F=frameshift,S=premature stop codon), c[12]=S_DESCRIPTION
	perc_s_len = "{0:.2f}".format(c[8])
	lastgffcolList=[str(c[3]),str(c[7]),str(perc_s_len),str(c[10]),str(c[11]),str(c[5]),str(c[12])]
	lastgffcol="_".join(lastgffcolList)
	gffList=[c[0],"LASTsearch","DNA",str(c[1]),str(c[2]),str(c[6]),str(c[4]),".",lastgffcol]
	gffString = "\t".join(gffList)
	of.write(gffString)
	of.write("\n")

of.close()

