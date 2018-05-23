# -*- coding: utf-8 -*-
#python
#blast_cleaner: script to take blast output and make into a best-hit fasta file
#Felix Beaudry 22 May 2018

import fileinput, argparse
from numpy import array, zeros, concatenate
import sys

def arguments():
        parser  = argparse.ArgumentParser(description="take blast output and make into a best-hit fasta file")
        parser.add_argument("-i","--input",help="input",required=True)
        parser.add_argument("-o","--outdir",help="output directory",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
inpath = args.input
outpath = args.outdir
infile = open(inpath,'r')

def blastParser(incoming_line):
	incoming_line.rstrip()
	early = incoming_line.split('>')
	earlyString = str(early[1])
	late = earlyString.split('<')
	clean = str(late[0])
	return clean

sequence = ["N"]
start = 0
qlength = 0
hitStart= 0
hlength = 0
extraSpots = [1]
locus_name = ""

for line in infile:
	if "<Iteration_query-def>" in line :
		#recent variables at next locus
		sequence = ["N"]
		start = 0
		qlength = 0
		hlength = 0
		extraSpots = [1]
		#print locus name (will later turn this into file name)
		locus_name = blastParser(line)

	if "<Iteration_query-len>" in line :
		#reference length
		qlength = int(blastParser(line))

	if "<Hsp_hit-from>" in line:
		hitStart = int(blastParser(line))

	if "<Hsp_hit-to>" in line:
		hlength = int(blastParser(line)) - hitStart


	if "<Hsp_query-from>" in line:
		#where to start filling sequence array with hit sequence
		start = int(blastParser(line)) - 1  

	#make array of spots where the query was extended
	if "<Hsp_qseq>" in line:
		#split query where it has been extended
		lineSplit = blastParser(line).split("-")
		#make array the right length
		extraSpots = extraSpots * (len(lineSplit) -1)

		for j in range(0,len(lineSplit)-1,1):
			extraSpots[j] = len(lineSplit[j]) + extraSpots[j-1]
			if extraSpots[j] == extraSpots[j-1]:
				extraSpots[j] = extraSpots[j] +1

				#make reference sequence full of Ns
		if len(sequence) == 1:
			if hlength-qlength >= 0:
				sequence = sequence * (hlength + len(extraSpots))
			else:
				sequence = sequence * (qlength + len(extraSpots))
				
				
	if "<Hsp_hseq>" in line :
		#split hit into array
		lineSplit = list(blastParser(line))
		
		position = 0 


		if len(sequence) >= len(lineSplit):
			for i in range(0, len(lineSplit) -1 , 1):
				sequence[start+position] = lineSplit[i]
				position += 1

		if len(sequence) < len(lineSplit):
		#remove spots in the hit where the query was extended
			for i in range(0, len(sequence) -1 , 1):
				sequence[start+position] = lineSplit[i]
				position += 1
			for j in range(0,len(extraSpots),1):
				del sequence[extraSpots[j]]



	if "</Hit>" in line :	
		file_name = outpath+'/'+locus_name+'.fasta'
		file = open(file_name,'w') 
		file.write(">roth\n")
		for i in range(0, len(sequence), 1):
			file.write(sequence[i])
			if (i % 60) == 0 & i != 0:
				file.write('\n') 
		file.write('\n') 
		file.close() 

#print as cells, with frames of sixty




		


