# -*- coding: utf-8 -*-
#python
#blast_cleaner: script to take blast output and make into a best-hit fasta file
#Felix Beaudry 22 May 2018

import fileinput, argparse
from numpy import array, zeros

def arguments():
        parser  = argparse.ArgumentParser(description="take blast output and make into a best-hit fasta file")
        parser.add_argument("-i","--input",help="input",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
inpath = args.input
infile = open(inpath,'r')

lociArray = (['empty'],['empty'])
rowCount = 0

for line in infile:
	if "<Iteration_query-def>" in line :
		line.rstrip()
		early = line.split('>')
		earlyString = str(early[1])
		late = earlyString.split('<')
		
		lociArray.append(str(late[0]))
		rowCount += 1
		#cut at ">" and "<"
		#store
	if "<Hsp_hseq>" in line :
		early = line.split('>')
		earlyString = str(early[1])
		late = earlyString.split('<')

		#make two celled array of name and sequence (including adding the sequence together)
		#append the two-celled array to the overall array

		
		name = str(lociArray[rowCount])
		file_name =  "test_files/roth/"+name+".fasta"

#print as cells, with frames of sixty

		file = open(file_name,'w') 
		file.write(">roth\n")
		file.write(late[0]) 
		file.close() 


		


