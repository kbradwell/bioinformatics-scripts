#! /usr/bin/python

'''
Script parses yn00 results.
Input: 
yn00 results file 
Output:
A single "detailed" line with all the NG86 and yn00 results: NG86 omega, dN, dS, and YN00 S, N, t, kappa, omega, dN, dN SE, dS, dS SE
A single "brief" line with just the NG86 dN, dS and omega and yn00 dN, dS and omega.

Functions adapted from Module Bio.Phylo.PAML._parse_yn00 Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com) 
'''

from sys import argv
import string
import re
import os

script,results_file=argv

def parse_ng86(lines, results): #Nei & Gojobori (1986) results
    startFlag=0
    sequences = []
    for line in lines:
	if "runmode" in line:
		startFlag=1
		continue
        line_floats_res = re.findall("-*\d+\.\d+(?!_)", line)
	#print "line_float_res",line_floats_res
        line_floats = [float(val) for val in line_floats_res]
	#print "line floats",line_floats
	if startFlag==1:
        	matrix_row_res = re.match("([^\s]+)", line)
	else:
		continue
	#print "line",line
        if matrix_row_res is not None:
	    #print "not none",matrix_row_res.group(1)
            seq_name = matrix_row_res.group(1).strip()
	    #print "seq name",seq_name
            sequences.append(seq_name)
            results[seq_name] = {}
            for i in range(0, len(line_floats), 3):
                NG86 = {}
                NG86["omega"] = line_floats[i]
                NG86["dN"] = line_floats[i+1]
                NG86["dS"] = line_floats[i+2]
                results[seq_name][sequences[i//3]] = {"NG86": NG86}
                results[sequences[i//3]][seq_name] = {"NG86": NG86}
    return (results, sequences)


def parse_yn00(lines, results, sequences): # Yang & Nielsen (2000) results
    for line in lines:
        line_floats_res = re.findall("-*\d+\.\d+", line)
        line_floats = [float(val) for val in line_floats_res]
        row_res = re.match("\s+(\d+)\s+(\d+)", line)
        if row_res is not None:
            seq1 = int(row_res.group(1))
	    #print "seq1,",seq1
            seq2 = int(row_res.group(2))
            #print "seq2",seq2
            seq_name1 = sequences[seq1-1]
	    #print "seqname1",seq_name1
            seq_name2 = sequences[seq2-1]
            YN00 = {}
            YN00["S"] = line_floats[0]
            YN00["N"] = line_floats[1]
            YN00["t"] = line_floats[2]
            YN00["kappa"] = line_floats[3]
            YN00["omega"] = line_floats[4]
            YN00["dN"] = line_floats[5]
            YN00["dN SE"] = line_floats[6]
            YN00["dS"] = line_floats[7]
            YN00["dS SE"] = line_floats[8]
            results[seq_name1][seq_name2]["YN00"] = YN00
            results[seq_name2][seq_name1]["YN00"] = YN00
            seq_name1 = None
            seq_name2 = None
    return results


"""Parse yn00 results file.""" 
results = {} 
with open(results_file) as handle: 
	lines = handle.readlines() 
for line_num in range(len(lines)): 
	line = lines[line_num] 
	if "(A) Nei-Gojobori (1986) method" in line: 
		#print "NG"
		ng86_start = line_num + 1 
	elif "(B) Yang & Nielsen (2000) method" in line:
		#print "YN" 
		(results, sequences) = parse_ng86(lines[ng86_start:line_num], results) 
	        #print sequences
		yn00_start = line_num + 1 
	elif "(C) LWL85, LPB93 & LWLm methods" in line: 
		results = parse_yn00(lines[yn00_start:line_num], results, sequences) 
if len(results) == 0: 
	raise ValueError("Invalid results file.") 

seen=[]
for g in results:
	for h in results[g]:
		tup = (g,h)
		tup2 = (h,g)
		if tup not in seen and tup2 not in seen:
			#print str(results_file),g,h,results[g][h]["YN00"],results[g][h]["NG86"]
			outList= ["detailed",str(results_file),g,h,results[g][h]["NG86"]["omega"],results[g][h]["NG86"]["dN"],\
			results[g][h]["NG86"]["dS"],results[g][h]["YN00"]["S"],results[g][h]["YN00"]["N"],\
			results[g][h]["YN00"]["t"],results[g][h]["YN00"]["kappa"],results[g][h]["YN00"]["omega"],\
			results[g][h]["YN00"]["dN"],results[g][h]["YN00"]["dN SE"],results[g][h]["YN00"]["dS"],\
			results[g][h]["YN00"]["dS SE"]]
			outStr=[str(i) for i in outList]
			outRes=" ".join(outStr)
			print outRes
			# input file name must have the correct format to produce this final formatted results file
			outList2=["brief",str(results_file),g,h,results[g][h]["NG86"]["dN"],results[g][h]["NG86"]["dS"],\
                        results[g][h]["NG86"]["omega"],results[g][h]["YN00"]["dN"],results[g][h]["YN00"]["dS"],\
                        results[g][h]["YN00"]["omega"]]
			outStr2=[str(i) for i in outList2]
			outRes2=" ".join(outStr2)
			print outRes2
		seen.append(tup)
		seen.append(tup2)
