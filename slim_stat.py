# -*- coding: utf-8 -*-
#python
#Felix Beaudry Aug 10 2018

import fileinput, argparse
import sys
import numpy
import math

def arguments():
        parser  = argparse.ArgumentParser(description="script for script slim scripts")
        parser.add_argument("-i","--input",help="input file",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
infile = args.input

with open(infile) as f:
    content = f.readlines()
content = [line.strip() for line in content] 
floaties = [float(i) for i in content]
sumOf = sum(floaties)
mean = sumOf / len(content)
	
slope = (floaties[len(content) -1 ] - floaties[0]) / len(content)

#se
totalDev = 0
for x in range(0,(len(content) -1)):
	dev = (floaties[x] - mean)**2
	totalDev = totalDev + dev
sd = math.sqrt(totalDev / (len(content) -1))
se = sd / math.sqrt((len(content) -1))
print mean, " ", se ," ", slope



