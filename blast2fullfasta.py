# -*- coding: utf-8 -*-
#python
#blast_cleaner: script to take blast output and make into a best-hit fasta file
#python blast2fullfasta.py -i <outgroup>.blast -o . -s <outgroup>

#Felix Beaudry 22 May 2018

import fileinput, argparse
from numpy import array, zeros, concatenate
import sys

def arguments():
        parser  = argparse.ArgumentParser(description="take blast output and make into a best-hit fasta file")
        parser.add_argument("-i","--input",help="input",required=True)
        parser.add_argument("-o","--outdir",help="output directory",required=True)
        parser.add_argument("-s","--outgroup",help="name for outgroup sequence",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
inpath = args.input
outpath = args.outdir
outgroup = args.outgroup
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
hitCount = 0

for line in infile:
	if "<Iteration_query-def>" in line :
		#recent variables at next locus
		sequence = ["N"]
		start = 0
		qlength = 0
		hlength = 0
		extraSpots = [1]
		locus_name = blastParser(line)
		#print locus_name

	if "<Iteration_query-len>" in line :
		#reference length
		qlength = int(blastParser(line))
		#print "qlength"
		#print qlength

	if "<Hsp_num>" in line :
		hitCount = int(blastParser(line))

	if "<Hsp_query-from>" in line:
		#where to start filling sequence array with hit sequence
		if hitCount == 1 :
			start = int(blastParser(line)) - 1 
			#print "start"
			#print start
		else :
			pass 
		

	if "<Hsp_hit-from>" in line:
		if hitCount == 1 :
			hitStart = int(blastParser(line))
		else :
			pass 

	if "<Hsp_hit-to>" in line:
		if hitCount == 1 :
			hlength = int(blastParser(line)) - hitStart
			#print "hlength"
			#print hlength
		else :
			pass 

	#make array of spots where the query was extended
	if "<Hsp_qseq>" in line:
		if hitCount == 1 :
			#split query where it has been extended
			lineSplit = blastParser(line).split("-")
			#make array the right length
			extraSpots = extraSpots * (len(lineSplit) -1)

			for j in range(0,len(extraSpots),1):
				if j == 0:
					extraSpots[j] = len(lineSplit[j])
				elif extraSpots[j] == extraSpots[j-1]:
					extraSpots[j] = extraSpots[j-1] +1
				else:
					extraSpots[j] = len(lineSplit[j]) + extraSpots[j-1]


					#make reference sequence full of Ns
			if len(sequence) == 1:
				if hlength-qlength >= 0:
					sequence = sequence * (hlength + len(extraSpots))
				else:
					sequence = sequence * (qlength + len(extraSpots))
		else :
			pass 
				
	if "<Hsp_hseq>" in line :
		if hitCount == 1 :
			#split hit into array
			lineSplit = list(blastParser(line))
			
			position = 0 
			for i in range(0, len(lineSplit), 1):
				sequence[start+position] = lineSplit[i]
				position += 1

			for j in range(0,len(extraSpots),1):
				del sequence[int(extraSpots[j])]

			if hlength-qlength >= 0:
				while qlength < len(sequence):
					del sequence[len(sequence)-1]

		else :
			pass
					
	if "</Hit>" in line :
		if len(sequence) > 1:
			file_name = outpath+'/'+locus_name+'.fasta'
			file = open(file_name,'w')
			outgroup_name = ">"+str(outgroup)+"\n"
			file.write(outgroup_name)
			#sys.stdout.write('\n')
			for i in range(0, len(sequence), 1):
				file.write(sequence[i])
				#sys.stdout.write(sequence[i])
				if ((i+1)%60 == 0 and i != 0):
					#sys.stdout.write('\n')
					file.write('\n') 
			#sys.stdout.write('\n\n')
			file.write('\n\n') 
			file.close() 






		


