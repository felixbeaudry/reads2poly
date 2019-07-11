
#python
#for converting .maf into .vcf-type vertical variant file
#Felix Beaudry April 16 2019

import re, math, fileinput, argparse, sys
import numpy as np

def arguments():
	parser  = argparse.ArgumentParser(description="converting .maf into .vcf-type vertical variant file")
	parser.add_argument("-i","--input",help="input .maf file",required=True)
	parser.add_argument("-q","--query",help="query assembly name",required=True)
	parser.add_argument("-p","--hit",help="hit assembly name",required=True)
	args = parser.parse_args()
	return(args)

args = arguments()
inpath = args.input
infile = open(inpath,'r')
queryName = str(args.query)
hitName = str(args.hit)

#set variables
result_array = np.array(["Query\tq_pos\tq_ncl\tq_alt\tHit\th_pos\th_ncl\th_alt\n"])
result_array = result_array.astype('U256')
content = infile.readlines()
qAss = 'empty'
pAss = 'empty'
dummyVar = "T"
qLength = 1

for line in content:  								#for each line in .maf
	if re.match('^#(.*)', line):					#skip header
	    continue 
  	if re.match('a', line): 						#reset variables
  		if re.match(str(dummyVar), str("T")): 		#print header
  			sys.stdout.write("\t".join(result_array))
  			dummyVar = "F"
  		else: 										#print matrix from last round
  			for i in range(0,qLength):
  				sys.stdout.write("\t".join(result_array[i]))
  				sys.stdout.write("\n")
   		qAss = 'empty'
   		pAss = 'empty'
	if re.match('s(.*)', line):
		maf_line = line.rstrip().split()
		#AssScaf = maf_line[1].split('.') 
		Ass = maf_line[1] #Ass = AssScaf[0]
		Scaf = maf_line[1] #Scaf = AssScaf[1]
		if re.match(str(qAss),'empty'):			#query
			qAss = Ass
			qString = list(maf_line[6])
			result_array = np.zeros([int(maf_line[3]),8])
			result_array = result_array.astype('U256')
			qLength = int(maf_line[3])   	
			for i in range(0,qLength):
				result_array[i,0] = Scaf
				result_array[i,1] = str(int(maf_line[2]) + i)
				result_array[i,2] = qString[i]
		elif re.match(str(qAss), str(Ass)):		#extra hits in query
			for i in range(0,qLength):
				result_array[i,3] = str(float(result_array[i,3]) + 1)
		elif re.match(str(pAss),'empty'):			#first hit
			pAss = Ass
			pString = list(maf_line[6])
			for i in range(0,qLength):
				result_array[i,4] = Scaf
				result_array[i,5] = str(int(maf_line[2]) + i)
				result_array[i,6] = pString[i]
		else:										#extra hits
			for i in range(0,qLength):
				result_array[i,7] = str(float(result_array[i,7]) + 1)

			


