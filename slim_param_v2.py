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
        parser.add_argument("-r","--recomb",help="recombination rate",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
le = args.len
ch = args.chrom
re = args.recomb

print('initialize() {')
print('initializeMutationRate(1e-7);')
print('initializeMutationType("m1", 0.5, "g", -0.0368,3.022); //deleterious')
print('initializeMutationType("m2", 0.5, "f", 0.05); //nearly neutral')
print('initializeGenomicElementType("g1", c(m1,m2), c(0.142,0.5));')
print('initializeGenomicElement(g1, 0, '+le+');')
print('initializeRecombinationRate('+re+');')
print('}')
print('1 {')
print('sim.addSubpop("p1", 500);')
print('}')
print('6000 late() {')
print('s = 0.05;')
print('N = 500;')
print('p_fix = (1 - exp(-2 * s)) / (1 - exp(-4 * N * s));')
print('n_gens = 1000;')
print('mu = 1e-7;')
print('locus_size = '+(le+1)+';')
print('expected = mu * locus_size * n_gens * 2 * N * p_fix;')
print('subs = sim.substitutions;')
print('actual = sum(subs.fixationGeneration >= 5000);')
print('cat("ratio: "+(actual/expected) + "\\n");')
print('}')