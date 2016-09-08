#! /usr/bin/python

from sys import argv
import string
import re
import os
import itertools

'''
Script makes the .tbl file, which is necessary for tbl2asn script input (GenBank data deposition).
~Input~ 
Various files containing different types of annotation information
Example files that all need to be present for this version of the script to work:
tRNA locations:
EXAMPLE_INPUT_annotation2tbl_tRNA.gff3
rRNA locations:
EXAMPLE_INPUT_annotation2tbl_rna.gff
protein-coding gene locations:
EXAMPLE_INPUT_annotation2tbl_genes.gff
Genes fasta:
EXAMPLE_INPUT_annotation2tbl_genes.fasta
EC #s and best nr BLAST hits:
EXAMPLE_INPUT_annotation2tbl_genes_supp_tbl
Intergenic repeats:
EXAMPLE_INPUT_annotation2tbl_ig-reps.outtable
Pseudogenes (all):
EXAMPLE_INPUT_annotation2tbl_pseudogenes.gff
EC number lookup file:
EXAMPLE_INPUT_annotation2tbl_genes.fasta.paths.detail
File with necessary prefixes:
EXAMPLE_INPUT_annotation2tbl_prefixes2use.txt, which contains Locus tag prefix, ProteinID prefix, Transcript ID prefix
e.g.
locus_tag       miscORG
protein_id      gnl|xyzVCU|
transcript_id   gnl|xyzVCU|mrna.
~Output~ 
.tbl file
Lookup table for gene names to locus tags and coordinates of other features to locus tags
'''

def getUniqECs(ECs):

	#make into a list of lists
	LoL=[]
	for e in ECs:
		es=e.strip(".")
		eList=es.split(".")
		#remove all the dashes
		eL=[ x for x in eList if "-" not in x ]
		#print eL
		LoL.append(eL)
	#print LoL
	
	for i in LoL:
		for j in LoL:
			if i==j:
				pass
			else:
				if len(i)<len(j):
					comp=j[0:len(i)]
					if i==comp:
						if i in LoL:
							LoL.remove(i)
				elif len(j)<len(i):
					comp=i[0:len(j)]
					if j==comp:
						if j in LoL:
							LoL.remove(j)
				else:
					pass
	#print LoL	
	LoL.sort()
	newEClist=list(LoL for LoL,_ in itertools.groupby(LoL))
	#add "-" to fill spaces up to 4 positions, where necessary
	newECsStr=[]
	for enz in newEClist:
		if len(enz)<4:
			dashes2add=4-len(enz)
			for i in range(dashes2add):
				enz.append("-")
		#turn list into a string again
		ecStr=".".join(enz)
		newECsStr.append(ecStr)
	return newECsStr

script, outfile, mappingFile = argv
# get the relevant file names that contain the annotation information
for file in os.listdir("./"):
	if file.endswith("tRNA.gff3"):
		print "tRNA gff3 file:", file
		trnaGFF=file
	if file.endswith("rna.gff"):
		print "rRNA gff file:", file
		rrnaGFF=file
	if file.endswith("_genes.gff"):
		print "genes gff file:", file
		genesGFF=file
	if file.endswith("genes.fasta"):
		print "genes nt fasta file:", file
		genesFasta=file
	if file.endswith("genes_supp_tbl"):
		print "EC no. and nr BLAST best hit info file - supp genes tbl:", file
		suppTbl=file
	if file.endswith("ig-reps.outtable"):
		print "Intergenic repeats file:", file
		igReps=file
	if file.endswith("pseudogenes.gff"):
		print "pseudogenes gff file:", file
		pseudosGFF=file
	if file.endswith("genes.fasta.paths.detail"):
		print "EC number: description lookup file:", file
		ecDescs=file
	if file.endswith("prefixes2use.txt"):
		print "prefixes file:", file
		prefixNames=file


print "............................"

# get the prefix names
prefixDict={}
prefixesList = open(prefixNames)
lines = prefixesList.readlines()
for line in lines:
	line = line.strip()
	cols = line.split("\t")
	prefixType=cols[0]
	pname=cols[1]
	#print prefixType,pname
	if prefixType not in prefixDict:
		prefixDict[prefixType]=pname
		
prefixesList.close()
	
# save EC: enzyme description info. from paths.detail file in case nr BLAST hit is "hypothetical" or "unknown" yet an EC number is present

ecDescDict={}
ecpathsFile = open(ecDescs)
lines = ecpathsFile.readlines()
for line in lines:
	if line.startswith(">"):
		pass
	elif line.startswith("#"):
		pass
	elif line.startswith("K"):
		# get the EC number from the description
		if "[EC:" in line:
			line = line.strip()
			cols = line.split("\t")
			ecDescription=cols[1].strip()
			kNumLineParts=ecDescription.split("[EC:")
			ecDescription=kNumLineParts[0].strip()
			ecbracket=kNumLineParts[1]
			ecIDparts=ecbracket.split("]")
			ecID=ecIDparts[0].strip()
			if ecID not in ecDescDict:
				ecDescDict[ecID]=ecDescription
	else:
		line = line.strip()
		cols = line.split("\t")
		ecID=str(cols[0]).strip()
		ecDescription=cols[1]
		if ecID not in ecDescDict:
			ecDescDict[ecID]=ecDescription

ecpathsFile.close()

# save any nr BLAST best hit and EC number information from the supplementary tables for all the protein coding genes
# remove any EC numbers when the product description is "unknown" or "hypothetical protein" and an EC number description is not present
supp = open(suppTbl)
lines = supp.readlines()

gInfo={} # holds unique EC numbers and nr (gb and ref) best hit info for genes
for line in lines:
	line = line.strip()
	cols = line.split("\t")
	genename = cols[0].strip()
	dbAnnotation = cols[13]
	dbID=cols[10]
	nrTup=("--","--")
	if "--" not in dbID:
		dbIDparts=dbID.split("|")
		dbShortName=dbIDparts[2]
		dbAccession=dbIDparts[3]
		# make sure GB is changed to INSD, and ref is changed to RefSeq, as per NCBI instructions
		if "GB" in dbShortName.upper():
			dbShortName="INSD"
		elif "REF" in dbShortName.upper():
			dbShortName="RefSeq"
		else:
			print "not valid database hit - has to be gb or ref hit"
		similarity=dbShortName+":"+dbAccession
		# remove the species in brackets from the annotation
		nrBLASThit1=re.sub('\[.*?\]','', dbAnnotation)
		# also make sure to remove any appearent EC number, inference, or info on partial or putative nature in protein description 
		nrBLASThit2=re.sub('\(EC\s{0,1}[0-9]{1,3}\..*?\s','',nrBLASThit1)
		nrBLASThit3=re.sub('\(inferred from BLASTP similarity\)','',nrBLASThit2)
		nrBLASThit4=re.sub(', putative, partial','',nrBLASThit3)
		nrBLASThit=re.sub(', partial','',nrBLASThit4)
		nrHit=nrBLASThit.strip()
		if "hypothetical protein" in nrHit:
			nrHit="hypothetical protein"
		nrTup=(nrHit,similarity)
	ecNumsInfo=cols[15].strip()
	ecNums=ecNumsInfo.split(",")
	#print "ecNums",ecNums
	#check all EC number are unique (no subsets e.g. enzyme classes together with lower level EC numbers belonging to them)
	if "--" not in ecNums:
		uniqECs=getUniqECs(ecNums)
		# check that product description allows for an EC number i.e. not "hypothetical" or "unknown", if not, use the EC number lookup description
		if "hypothetical" in nrHit.lower() or "unknown" in nrHit.lower() or "predicted protein" in nrHit.lower() or "related to" in nrHit.lower() or "--" in nrTup[0]:
			newEnzDesc=[]
			for enz in uniqECs:
				if enz in ecDescDict:
					enzymedesc1=str(ecDescDict[enz]).strip(".")
					enzymedesc=enzymedesc1.rstrip("s").lower()
					newEnzDesc.append(enzymedesc)	
			if len(newEnzDesc)>0:
				nrHit="/".join(newEnzDesc)
				nrTup=(nrHit,"NONE")
			else: # if there is no database hit or only hypothetical, unknown or predicted protein as desc, and no alt enzyme desc from paths.detail
				uniqECs=["--"]
	else:
		uniqECs=ecNums
	if genename not in gInfo:
		gInfo[genename]={}
		gInfo[genename]["EC"]=uniqECs
		gInfo[genename]["bestHit"]=nrTup
	else:
		print "ERROR: duplicated gene info in supplementary table!"
	#print ecNums
	#print similarity
	
supp.close()

# reading the fasta file to get the sequence info in a dict (used later to determine partial genes)
# code for gene ID : gene seq mapping modified from https://github.com/alunem/bioman/blob/master/fasta/bmn-SubsetFasta.py
handle = open(genesFasta, "r")

storedSeqs={}
seqIDmap={}
seq_id = handle.next()
while (seq_id[0]!=">"):
    seq_id = handle.next()
while True:
    try:
        seq = handle.next()
        line = handle.next()
        while (line[0]!=">"):
            seq = seq+line
            line = handle.next()
        storedSeqs[seq_id]=seq
        shortID=string.split(seq_id, " ")[0][1:].strip()
        seqIDmap[shortID]=seq_id
        seq_id = line 
    except StopIteration:
        break
# store last line
storedSeqs[seq_id]=seq
seqIDmap[string.split(seq_id, " ")[0][1:].strip()]=seq_id

# seqIDmap[g] will be used to get the gene name
# storedSeqs[seqIDmap[g]].strip() will be used to get the sequence

handle.close()

# go through each of the annotation files and get the coordinates and orientation of each feature, then save in data structure ("tblDict") that mirrors the .tbl file format
# use of lists in the tblDict is to get the feature keys and qualifier keys to print in the right order
# for protein coding genes check the fasta sequence to determine partial genes, add ">", "<", or notes as necessary

tblDict={}

#tRNA
trnaInfo = open(trnaGFF)
lines = trnaInfo.readlines()[1:]
for line in lines:
	line = line.strip()
	cols = line.split("\t")
	#print cols
	col1keyTup=()
	if "gene" in cols[2]:
		contig = cols[0]
		start=cols[3]
		end=cols[4]
		strand=cols[6]
		infoLine=cols[8]
		infoLineParts=infoLine.split(";")
		codonParts=infoLineParts[0].split("=")
		codonNumbered=codonParts[1]
		codon=codonNumbered.split("-")[0]
		trnaProduct=codon.replace(":","-")
		inferenceStr="COORDINATES:profile:tRNAscan-SE:2.19"
		col1keyTup=(start,end,strand)
		if contig not in tblDict:
			#print "new key for contig: ", col1keyTup
			tblDict[contig]={}
		if col1keyTup in tblDict[contig]:
			print "ERROR: feature start, end and orientation have exact overlap - check for duplication", col1keyTup, contig, trnaProduct
		else:
			#print "new key for contig: ", col1keyTup
			tblDict[contig][col1keyTup]=[{},{}] # tRNA coordinates contain two feature keys: "gene" and "tRNA"
			tblDict[contig][col1keyTup][0]['gene']=[{}] # gene for tRNA has one qualifier key: "locus_tag"
			tblDict[contig][col1keyTup][0]['gene'][0]['locus_tag']="X"  # give the locus_tag a dummy value for now - the locus tags will be assigned at the end of the script
			tblDict[contig][col1keyTup][1]['tRNA']=[{},{}]
			tblDict[contig][col1keyTup][1]['tRNA'][0]['product']=trnaProduct
			tblDict[contig][col1keyTup][1]['tRNA'][1]['inference']=inferenceStr

print "tRNA features have been added"				
trnaInfo.close()

#rRNA
rrnaInfo = open(rrnaGFF)
lines = rrnaInfo.readlines()
rrnaCoordsList = [] # this will save the rRNA coordinates so that any pseudogene overlapping them can be removed from the annotation
for line in lines:
	line = line.strip()
	cols = line.split("\t")
	col1keyTup=()
	if "#" not in line:
		contig = cols[0]
		start=cols[3]
		end=cols[4]
		# save to the coords list so that any pseudogenes overlapping them can be removed
		rcoords=[contig,start,end]
		rrnaCoordsList.append(rcoords)
		strand=cols[6]
		rRNAtype=cols[8]
		typeParts=rRNAtype.split("_")
		rRNA=typeParts[0].strip()
		rRNAname=rRNA.upper()
		if rRNAname=="8S":
			rRNAname="5S"
		rrnaProduct=rRNAname+" ribosomal RNA"
		inferenceStr="COORDINATES:profile:RNAmmer:1.2"
		#print start,end,strand, contig,rrnaProduct
		col1keyTup=(start,end,strand)
		if contig not in tblDict:
			tblDict[contig]={}
		if col1keyTup in tblDict[contig]:
			print "ERROR: feature start, end and orientation have exact overlap - check for duplication", col1keyTup, contig, rrnaProduct
		else:
			tblDict[contig][col1keyTup]=[{},{}]
			tblDict[contig][col1keyTup][0]['gene']=[{}]
			tblDict[contig][col1keyTup][0]['gene'][0]['locus_tag']="X"
			tblDict[contig][col1keyTup][1]['rRNA']=[{},{}]
			tblDict[contig][col1keyTup][1]['rRNA'][0]['product']=rrnaProduct
			tblDict[contig][col1keyTup][1]['rRNA'][1]['inference']=inferenceStr

print "rRNA features have been added"
rrnaInfo.close()

#protein-coding genes
geneCoordDict={} # used for mapping locus tags back to genenames later on, in the mapping file
stopCodons=["TAA","TGA","TAG"]
genesInfo = open(genesGFF)
lines = genesInfo.readlines()
for line in lines:
	partial5prime=""
	partial3prime=""
	line = line.strip()
	cols = line.split("\t")
	col1keyTup=()
	contig = cols[0]
	#print contig
	start=cols[3]
	end=cols[4]
	strand=cols[6]
	genename=cols[8]
	#check if the gene is partial or not
	sequence=storedSeqs[seqIDmap[genename]].strip()
	ntSeq = sequence.replace('\n', '').replace('\r', '')
	initiationCodon=ntSeq[0:3].upper()
	if initiationCodon != "ATG":
		partial5prime="<"
	terminationCodon=ntSeq[-3:].upper()
	if not any(x in terminationCodon for x in stopCodons):
		partial3prime=">"
	#print genename, initiationCodon, partial5prime, terminationCodon, partial3prime, "sequence"
	if strand=="+":
		s=partial5prime+start
		e=partial3prime+end
	else:
		s=partial3prime+start
		e=partial5prime+end
	col1keyTup=(s,e,strand)
	if genename not in geneCoordDict:
		geneCoordDict[col1keyTup]=genename
	else:
		print genename, " gene seems to be duplicated in the genes GFF file"
	#print col1keyTup
	if contig not in tblDict:
		tblDict[contig]={}
	if col1keyTup in tblDict[contig]:
		print "ERROR: feature start, end and orientation have exact overlap - check for duplication", col1keyTup, contig, genename
	else:	
		tblDict[contig][col1keyTup]=[{},{},{}] # three dicts because each gene has "gene","mRNA" and "CDS" feature keys
		tblDict[contig][col1keyTup][0]['gene']=[{}]
		tblDict[contig][col1keyTup][0]['gene'][0]['locus_tag']="X"
		tblDict[contig][col1keyTup][1]['mRNA']=[{},{},{}]
		if "--" not in gInfo[genename]["bestHit"]:
			tblDict[contig][col1keyTup][1]['mRNA'][0]['product']=gInfo[genename]["bestHit"][0]
		else:
			tblDict[contig][col1keyTup][1]['mRNA'][0]['product']="hypothetical protein"
		tblDict[contig][col1keyTup][1]['mRNA'][1]['protein_id']="X"
		tblDict[contig][col1keyTup][1]['mRNA'][2]['transcript_id']="X"
		tblDict[contig][col1keyTup][2]['CDS']=[{}]
		if "--" not in gInfo[genename]["bestHit"]:
			tblDict[contig][col1keyTup][2]['CDS'][0]['product']=gInfo[genename]["bestHit"][0]
		else:
			tblDict[contig][col1keyTup][2]['CDS'][0]['product']="hypothetical protein"
		qualifierCounter=1
		if "--" not in gInfo[genename]["EC"]:
			qualifierCounter=0
			for i in xrange(len(gInfo[genename]["EC"])):
				qualifierCounter+=1
				tblDict[contig][col1keyTup][2]['CDS'].append({})
				tblDict[contig][col1keyTup][2]['CDS'][i+1]["EC_number"]=gInfo[genename]["EC"][i]
		tblDict[contig][col1keyTup][2]['CDS'].append({})
		tblDict[contig][col1keyTup][2]['CDS'][qualifierCounter]['protein_id']="X"
		tblDict[contig][col1keyTup][2]['CDS'].append({})
		tblDict[contig][col1keyTup][2]['CDS'][qualifierCounter+1]['transcript_id']="X"
		tblDict[contig][col1keyTup][2]['CDS'].append({})
		tblDict[contig][col1keyTup][2]['CDS'][qualifierCounter+2]['inference']="ab initio prediction:GeneMarkS:4.7b"
		if "--" not in gInfo[genename]["bestHit"]:
			tblDict[contig][col1keyTup][2]['CDS'].append({})
			if "INSD" in str(gInfo[genename]["bestHit"][1]):
				tblDict[contig][col1keyTup][2]['CDS'][qualifierCounter+3]['inference']="similar to AA sequence:"+str(gInfo[genename]["bestHit"][1])
				if qualifierCounter>1:
					tblDict[contig][col1keyTup][2]['CDS'].append({})
					tblDict[contig][col1keyTup][2]['CDS'][qualifierCounter+4]['note']="bifunctional"	
			if "RefSeq" in str(gInfo[genename]["bestHit"][1]):
				tblDict[contig][col1keyTup][2]['CDS'][qualifierCounter+3]['inference']="similar to RNA sequence, mRNA:"+str(gInfo[genename]["bestHit"][1])
				if qualifierCounter>1:
					tblDict[contig][col1keyTup][2]['CDS'].append({})
					tblDict[contig][col1keyTup][2]['CDS'][qualifierCounter+4]['note']="bifunctional"
			if "NONE" in str(gInfo[genename]["bestHit"][1]):
				# leave out the similarity part - the description will come from the EC number and not the nr results
				if qualifierCounter>1:
					tblDict[contig][col1keyTup][2]['CDS'].append({})
					tblDict[contig][col1keyTup][2]['CDS'][qualifierCounter+3]['note']="bifunctional"
		else: 
			if qualifierCounter>1:
				tblDict[contig][col1keyTup][2]['CDS'].append({})
				tblDict[contig][col1keyTup][2]['CDS'][qualifierCounter+3]['note']="bifunctional"
						
print "protein-coding features have been added"
genesInfo.close()

#intergenic repeats
repCats2save=["LINE","LTR"]
reps2save=["VIPER","SIRE","CZAR","L1Tc"]
repsInfo = open(igReps)
lines = repsInfo.readlines()
for line in lines:
	line = line.strip()
	cols = line.split("\t")
	repCat=cols[7]
	# get the mobile elements
	if any(x in repCat for x in repCats2save):
		repType=cols[8]
		#print repType
		if "TcVIPER" in repType:
			repType="VIPER"
		if "SIRE" in repType:
			repType="SIRE"
		if "CZAR" in repType:
			repType="CZAR"
		if "L1_TC" in repType:
			repType="L1Tc"
		if any(x in repType for x in reps2save):
			#print line
			contig=cols[0]
			start=cols[1]
			end=cols[2]
			strand=cols[3]
			if "C" in strand:
				strand="-"
			repName="retrotransposon:"+repType
			col1keyTup=(start,end,strand)
			if contig not in tblDict:
				tblDict[contig]={}
			if col1keyTup in tblDict[contig]:
				print "ERROR: feature start, end and orientation have exact overlap - check for duplication", col1keyTup, contig, repName
			else:
				tblDict[contig][col1keyTup]=[{}] # each repeat element has just one feature key: "repeat_region"
				tblDict[contig][col1keyTup][0]['mobile_element']=[{}] # just one qualifier key: "mobile_element"
				tblDict[contig][col1keyTup][0]['mobile_element'][0]['mobile_element_type']=repName
	# get the repeat_region features - simple repeats ("microsatellites") and satellites
	elif repCat=="Simple_repeat":
		contig=cols[0]
		start=cols[1]
		end=cols[2]
		strand=cols[3]
		if "C" in strand:
			strand="-"
		repName="microsatellite"
		repType=cols[8].strip()
		rptSeq=repType[repType.find("(")+1:repType.find(")")]
		rptUnitSeq=rptSeq.lower()
		#print rptUnitSeq
		col1keyTup=(start,end,strand)
		if contig not in tblDict:
			tblDict[contig]={}
		if col1keyTup in tblDict[contig]:
			print "ERROR: feature start, end and orientation have exact overlap - check for duplication", col1keyTup, contig, repName
		else:
			tblDict[contig][col1keyTup]=[{}] # each repeat element has just one feature key: "repeat_region"
			tblDict[contig][col1keyTup][0]['repeat_region']=[{},{},{}] # qualifier keys: "rpt_type", "rpt_unit_seq" and "satellite"
			tblDict[contig][col1keyTup][0]['repeat_region'][0]['rpt_type']="tandem"
			tblDict[contig][col1keyTup][0]['repeat_region'][1]['rpt_unit_seq']=rptUnitSeq
			tblDict[contig][col1keyTup][0]['repeat_region'][2]['satellite']=repName
	elif repCat=="Satellite":
		contig=cols[0]
		start=cols[1]
		end=cols[2]
		strand=cols[3]
		if "C" in strand:
			strand="-"
		repName="satellite"
		#repType=cols[8].strip()
		col1keyTup=(start,end,strand)
		if contig not in tblDict:
			tblDict[contig]={}
		if col1keyTup in tblDict[contig]:
			print "ERROR: feature start, end and orientation have exact overlap - check for duplication", col1keyTup, contig, repName
		else:
			tblDict[contig][col1keyTup]=[{}] # each repeat element has just one feature key: "repeat_region"
			tblDict[contig][col1keyTup][0]['repeat_region']=[{},{}] # qualifier keys: "rpt_type" and "satellite"
			tblDict[contig][col1keyTup][0]['repeat_region'][0]['rpt_type']="tandem"
			tblDict[contig][col1keyTup][0]['repeat_region'][1]['satellite']=repName
	elif repCat=="Low_complexity":
		contig=cols[0]
		start=cols[1]
		end=cols[2]
		strand=cols[3]
		if "C" in strand:
			strand="-"
		repName=cols[8].strip()
		#print repName
		col1keyTup=(start,end,strand)
		if contig not in tblDict:
			tblDict[contig]={}
		if col1keyTup in tblDict[contig]:
			print "ERROR: feature start, end and orientation have exact overlap - check for duplication", col1keyTup, contig, repName
		else:
			tblDict[contig][col1keyTup]=[{}] # each repeat element has just one feature key: "misc_feature"
			tblDict[contig][col1keyTup][0]['misc_feature']=[{}] # qualifier key: "note"
			tblDict[contig][col1keyTup][0]['misc_feature'][0]['note']=repName

print "repeat features, mobile elements (retrotransposons VIPER, SIRE, CZAR and L1Tc) and misc features (e.g. at-rich) have been added"
repsInfo.close()

#pseudogenes
pseudoInfo = open(pseudosGFF)
lines = pseudoInfo.readlines()
for line in lines:
	line = line.strip()
	cols = line.split("\t")
	contig=cols[0]
	start=cols[3]
	end=cols[4]
	# check if the pseudogenes overlap rRNA coordinates - if they do, don't include that pseudogene in the annotation
	pseudoflag=0
	for rcoord in rrnaCoordsList:
		if contig==rcoord[0]:
			pr=range(int(start),int(end))
			rr=range(int(rcoord[1]),int(rcoord[2]))
			prs=set(pr)
			pseudorrnaOverlap=prs.intersection(rr)
			if len(pseudorrnaOverlap)>0:
				print "pseudogene overlaps rRNA coordinate, skipping this pseudogene for final annotation", contig, start, end
				pseudoflag=1
	if pseudoflag==1:
		continue
	strand=cols[6]
	desc=cols[8].strip()
	descParts=desc.split("|")
	descName=descParts[4]
	nrBLASThit=re.sub('\[.*?\]','', descName)
	nrHit=nrBLASThit.strip("_")
	if "hypothetical protein" in nrHit:
		nrHit="hypothetical protein"
	nrHitName=re.sub('_',' ', nrHit)
	#print nrHitName
	pseudoStats=descParts[0].split("_")
	pseudoFeature=pseudoStats[4]
	if "F" in pseudoFeature and "S" in pseudoFeature:
		pseudoNote="frameshift and premature termination codon"
	elif "F" in pseudoFeature and "S" not in pseudoFeature:
		pseudoNote="frameshift"
	elif "S" in pseudoFeature and "F" not in pseudoFeature:
		pseudoNote="premature termination codon"
	col1keyTup=(start,end,strand)
	if contig not in tblDict:
		tblDict[contig]={}
	if col1keyTup in tblDict[contig]:
		print "ERROR: feature start, end and orientation have exact overlap - check for duplication", col1keyTup, contig, nrHitName	
	else:
		tblDict[contig][col1keyTup]=[{}] # just gene feature key for pseudogenes
		tblDict[contig][col1keyTup][0]['gene']=[{},{},{},{}] # gene desc, locus tag and note qualifier keys and pseudo tag
		tblDict[contig][col1keyTup][0]['gene'][0]['gene_desc']=nrHitName
		tblDict[contig][col1keyTup][0]['gene'][1]['locus_tag']="X"
		tblDict[contig][col1keyTup][0]['gene'][2]['pseudo']=""
		tblDict[contig][col1keyTup][0]['gene'][3]['note']=pseudoNote

print "pseudogene features have been added"
pseudoInfo.close()

# all the information has now been saved - the contigs and coordinates now have to be sorted and any locus tags required have to be assigned (and protein_id/transcript_ids modified to incl. locus tags)
# simple mapping file also generated at the end

ofh = open(outfile,"w")

locusTagGenenameMapping={}
locusTagCounter=0

for c in sorted(tblDict):
	#print c
	contigDict=tblDict[c]
	contigLine=">Feature "+c+"\n"
	ofh.write(contigLine)
	#sd=sorted(contigDict.items(), key=lambda x: (int(str(x[0][0]).strip("<")),x[1]))
	sd=sorted(contigDict.items(), key=lambda x: (int(re.sub('[<>]', '', str(x[0][0]))),x[1]))
	for k in  sd: # for each coordinate for that contig
		if k[0][2]=="-":
			coordStr=k[0][1]+"\t"+k[0][0]
		else:
			coordStr=k[0][0]+"\t"+k[0][1]
		for d in k[1]:
			for ky in d:
				featureKey=coordStr+"\t"+ky+"\n"
				ofh.write(featureKey)
				if "gene" in ky:
					locusTagCounter+=1
					ltcStr=str(locusTagCounter)
					ltcPadded=ltcStr.zfill(5)
					locusTag=prefixDict['locus_tag']+"_"+ltcPadded
					if k[0] in geneCoordDict:
						locusTagGenenameMapping[geneCoordDict[k[0]]]=locusTag #save genename:locus_tag info for mapping file
				for arrayElem in d[ky]:
					for aEkey in arrayElem:
						if "locus_tag" in aEkey:
							arrayElem[aEkey]=locusTag
						if "protein_id" in aEkey:
							arrayElem[aEkey]=prefixDict['protein_id']+locusTag
						if "transcript_id" in aEkey:
							arrayElem[aEkey]=prefixDict['transcript_id']+locusTag
						qualifierKeyVal="\t\t\t"+aEkey+"\t"+arrayElem[aEkey]+"\n"
						ofh.write(qualifierKeyVal)
		

ofh.close()

ofh2 = open(mappingFile,"w")

for genekey in sorted(locusTagGenenameMapping):
	outStr=genekey+"\t"+locusTagGenenameMapping[genekey]+"\n"
	ofh2.write(outStr)

ofh2.close()


