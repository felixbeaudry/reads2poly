# -*- coding: utf-8 -*-
#python
##script for script slim scripts
#Felix Beaudry June 22 2018

import fileinput, argparse
import sys

def arguments():
        parser  = argparse.ArgumentParser(description="script for script slim scripts")
        parser.add_argument("-c","--chrom",help="type of chromosome",required=True)
        parser.add_argument("-l","--len",help="number of basepairs",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
le = args.len
ch = args.chrom

print('initialize() {')
print('initializeMutationRate(1e-6);')
print('initializeMutationType("m1", 0.5, "f", 0.05);')
print('initializeGenomicElementType("g1", m1, 1.0);')
print('initializeGenomicElement(g1, 0, '+le+');')
print('initializeRecombinationRate(1e-8);')
print('initializeSex("'+ch+'");')
print('}')
print('1 {')
print('sim.addSubpop("p1", 1000);')
print('}')
print('6000 late() {')
print('s = 0.05; N = 1000; p_fix = (1 - exp(-2 * s)) / (1 - exp(-4 * N * s)); n_gens = 1000;   mu = 1e-6;  locus_size = 100000;')
print('expected = mu * locus_size * n_gens * 2 * N * p_fix;')
print('subs = sim.substitutions;')
print('actual = sum(subs.fixationGeneration >= 5000);')
print('cat("ratio "+(actual/expected) + "\\n");')
print('}')