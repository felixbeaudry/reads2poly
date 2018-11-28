# -*- coding: utf-8 -*-
#python
##script for extracting Hudson data (from ms), edited from perl script by S. Wright
#Felix Beaudry May 28 2018

import fileinput, argparse
import sys

def arguments():
	parser  = argparse.ArgumentParser(description="extracting Hudson data (from ms)")
	parser.add_argument("-i","--input",help="input",required=True)
	parser.add_argument("-p","--pop",help="number of individuals in first population",required=True)
	parser.add_argument("-e","--end",help="total number of inds",required=True)
	args = parser.parse_args()
	return(args)

args = arguments()
inpath = args.input
infile = open(inpath,'r')
popSplit = int(args.pop)
popEnd = int(args.end)

def calculate_pi(seqs):
	diffs=0
	comps=0
	pi=0
	if len(seqs) > 0:
		for i in range(0,len(seqs),1):
			for t in range((i+1),len(seqs),1):
				comps+=1
				sqi = list(seqs[i])
				sqt = list(seqs[t])
				for j in range(0,len(seqs[0]),1):
					if sqi[j] != sqt[j]:
						diffs+=1
		if comps>0:
			pi= float(diffs)/float(comps)
		else :
			pi = 0
	else:
		pi = 0
	return pi

def calculate_dxy(seqOne,seqTwo):
	diffs=0
	comps=0
	dxy=0
	if (len(seqOne) + len(seqTwo)) > 0:
		for i in range(0,len(seqOne),1):
			sqO = list(seqOne[i])
			for t in range(0,len(seqTwo),1):
				sqT = list(seqTwo[t])
				comps+=1
				for j in range(0,len(seqOne[0]),1):
					if sqO[j] != sqT[j]:
						diffs+=1
		if comps>0:
			dxy= float(diffs)/float(comps)
		else :
			dxy = 0
	else:
		dxy = 0
	return dxy



sequence = [0]*popEnd
sequence_one = [0]*popSplit
sequence_two = [0]*(popEnd-popSplit)
i =0
pi_tot = 0
pi_one = 0
pi_two = 0
fst = 0
rep = 0
on = "F"
print "rep\tpi_one\tpi_two\tpi_tot\tdxy"; 

for line in infile:
	if not line.strip():
		on = "F"

	if "segsites" in line:
		segsite = line.split(': ')
		#print segsite[1]

	if "//"  in line:
		if rep == 0:
			rep+=1
		elif int(segsite[1]) == 0:
			rep+=1
			sys.stdout.write(str(rep-1)+"\tNA\tNA\tNA\n")

		else:
			rep+=1
 			pi_tot = calculate_pi(sequence)
 			pi_one = calculate_pi(sequence_one)
 			pi_two = calculate_pi(sequence_two)

 			dxyRep = calculate_dxy(sequence_one,sequence_two)
 			sys.stdout.write(str(rep-1)+"\t"+str(pi_one)+"\t"+str(pi_two)+"\t"+str(pi_tot)+"\t"+str(dxyRep)+"\n")

			sequence = [0]*popEnd
			sequence_one = [0]*popSplit
			sequence_two = [0]*(popEnd-popSplit)
 			pi_tot = 0
 			pi_one = 0
 			pi_two = 0
 			fst = 0
 			i = 0

	if on == "T":
		sequence[i]=line
		#print i
		#print sequence[i]
		if i < popSplit:
			sequence_one[i]=line
		else :
			sequence_two[i - popSplit]=line
		if i == popEnd:
			on == "F"
		i+=1

	if "positions" in line:
		on = "T"

