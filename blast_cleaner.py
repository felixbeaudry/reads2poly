# -*- coding: utf-8 -*-
#python
#blast_cleaner: script to take blast output and make into a best-hit fasta file
#Felix Beaudry 22 May 2018

import fileinput, argparse
from numpy import array, zeros, concatenate

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



for line in infile:
	if "<Iteration_query-def>" in line :
		
		#cut at ">" and "<"
		line.rstrip()
		early = line.split('>')
		earlyString = str(early[1])
		late = earlyString.split('<')
		#store
		
		locus_name = str(late[0])
		line_all = ''

	if "<Hsp_hseq>" in line :
		early = line.split('>')
		earlyString = str(early[1])
		late = earlyString.split('<')

		#make two celled array of name and sequence (including adding the sequence together)
		#append the two-celled array to the overall array
		line_clean = str(late[0])
		line_all = line_all + line_clean

	if "</Hit>" in line :
		
		file_name = outpath+locus_name+'.fasta'
		file = open(file_name,'w') 
		file.write(">roth\n")
		for i in range(0, len(line_all), 60):
			file.write(line_all[i:i+60])
			file.write('\n') 
		file.close() 

#print as cells, with frames of sixty




		


