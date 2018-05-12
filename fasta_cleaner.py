#python
#script to output only useful fasta sequences
#Felix Beaudry 8 May 2018

import fileinput, argparse
from numpy import array, zeros

def arguments():
        parser  = argparse.ArgumentParser(description="script to output only useful fasta sequences")
        parser.add_argument("-i","--input",help="input",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
inpath = args.input
try:
	infile = open(inpath,'r')
except IOError:
	print('There was an error opening the file!')
    return

sequence = array(
	[
	zeros(5),
	zeros(5)
	])
alignment = ['','']
x=-1
for line in infile:
	if ">" in line :
		x += 1

		values = line.split('_')
		ind = str(values[0])
		for w in range(1,len(values)-1):	
			ind = ind + "_" + str(values[w])
		sequence[x][0] = values[(len(values)-1)]
	else:
	 	sequence[x][1] = line.count('N') + sequence[x][1] 
	 	sequence[x][2] = line.count('-') + sequence[x][2] 
	 	sequence[x][3] = len(line) + sequence[x][3] 
	 	alignment[x] = alignment[x] + line

x=int(0)


for x in range(0, 2):
	if (sequence[x][3] >0) & ((sequence[x][1]/sequence[x][3]) < 0.6) & ((sequence[x][2]/sequence[x][3]) < 0.6):
		print str(ind)+"_"+str(int(sequence[x][0]))
		print alignment[x]
		x += 1
	else:
	 	x += 1


