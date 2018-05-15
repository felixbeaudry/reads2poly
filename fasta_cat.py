#python
#script to cat fasta sequences
#Felix Beaudry 15 May 2018

import fileinput, argparse, re
from numpy import array, zeros

def arguments():
        parser  = argparse.ArgumentParser(description="script to cat fasta sequences")
        parser.add_argument("-i","--input",help="input",required=True)
        parser.add_argument("-n","--ind",help="individual name",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
inpath = args.input
infile = open(inpath,'r')
indName = args.ind

print('>',indName)
#print ind name

for line in infile:
    if re.match(">",line):
        pass
    else:
		line = line.rstrip("\n")
		string = string + line
print (string)