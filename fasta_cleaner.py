#python
#script to output only useful fasta sequences
#Felix Beaudry 8 May 2018

import fileinput, argparse
import sys

def arguments():
        parser  = argparse.ArgumentParser(description="script to output only useful fasta sequences")
        parser.add_argument("-i","--input",help="input",required=True)
        parser.add_argument("-c","--cutoff",help="percent missing data cutoff",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
inpath = args.input
minC = args.cutoff
perc = (1- ( float(minC) / 100))

names = [None] * 2
alignment = [None] * 2
Ns = [None] * 2
length = [None] * 2

x=-1

def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

with open(inpath) as infile:
	for line in nonblank_lines(infile):
		if ">" in line :
			x += 1
			names[x] = str(line)
		else:
		 	if alignment[x] == None:
		 		alignment[x] =  str(line) 
		 		Ns[x] = float(line.count('N'))
		 		length[x] = float(len(line)) 
		 	else:
		 		alignment[x] =  str(alignment[x]) + str(line) 
		 		Ns[x] = float(line.count('N')) + float(Ns[x])
		 		length[x] = float(len(line)) + float(length[x])


for z in 0,1:
	if (length[z] > 0) :
		if( (Ns[z] / length[z]) < perc):
			print names[z]
			#sys.stdout.write('\n')
			curr = alignment[z]
			for i in range(0, len(curr), 1):
				sys.stdout.write(curr[i])
				if ((i+1)%60 == 0 and i != 0):
					sys.stdout.write('\n')
	sys.stdout.write('\n')
				
			



