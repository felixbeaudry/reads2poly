# -*- coding: utf-8 -*-
#python
##script for script slim scripts
#Felix Beaudry June 22 2018

import fileinput, argparse
import sys

def arguments():
        parser  = argparse.ArgumentParser(description="script for script slim scripts")
        parser.add_argument("-s","--sel",help="selection coefficient",required=True)
        parser.add_argument("-p","--pro",help="proportion of sites under selection",required=True)
        parser.add_argument("-t","--typ",help="direction of selection",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()

di = args.typ


if di == "n":
	se = str(-1 * ( 10  ** (float(args.sel) / 100) ) )
else :
	se = str(1 * ( 10  ** (float(args.sel) / 100) ) )

pr = str( (10  ** (float(args.pro) / 10) ) /10 )

print('initialize()')
print('{')
print('initializeMutationRate(1.5e-8);')
print('initializeMutationType("m1", 0.5, "f", 0.0);     // PAR')
print('initializeMutationType("m2", 0.5, "f", 0.0);     // non-PAR')
print('initializeMutationType("m3", 1.0, "f", 0.0);     // Y marker')
print('initializeMutationType("m4", 0.5, "f", '+se+');     // selection')
print('initializeGenomicElementType("g1", m1, 1.0);     // PAR: m1 only')
print('initializeGenomicElementType("g2", c(m2,m4), c('+pr+','+str((float(1)-float(pr)))+'));     // non-PAR: m2 only')
print('initializeGenomicElement(g1, 0, 2699999);        // PAR')
print('initializeGenomicElement(g2, 2700000, 5999999);  // non-PAR')
print('initializeSex("A");')
print('initializeRecombinationRate(c(1e-8, 0), c(2699999, 5999999), sex="M");')
print('initializeRecombinationRate(1e-8, sex="F");')
print('}')
print('1 late() {')
print('sim.addSubpop("p1", 1000);')
print('i = p1.individuals;')
print('males = (i.sex == "M");')
print('maleGenomes = i[males].genomes;')
print('yChromosomes = maleGenomes[rep(c(F,T), sum(males))];')
print('yChromosomes.addNewMutation(m3, 0.0, 5999999);')
print('}')
print('modifyChild() {')
print('if (child.sex == "F")')
print('{')
print('if (childGenome1.containsMarkerMutation(m3, 5999999))')
print('return F;')
print('if (childGenome2.containsMarkerMutation(m3, 5999999))')
print('return F;')
print('return T;')
print('}')
print('else')
print('{')
print('if (childGenome1.containsMarkerMutation(m3, 5999999))')
print('return T;')
print('if (childGenome2.containsMarkerMutation(m3, 5999999))')
print('return T;')
print('return F;')
print('}')
print('}')
print('1:10000 late() {')
print('if (sim.generation % 1000 == 0) {')
print('numY = sum(p1.individuals.sex == "M");')
print('numX = 2 * size(p1.individuals) - numY;')
print('firstMale = p1.individuals[p1.individuals.sex == "M"][0];')
print('fMG = firstMale.genomes;')
print('if (fMG[0].containsMarkerMutation(m3, 5999999)) {')
print('firstY = fMG[0];')
print('firstX = fMG[1];')
print('} else if (fMG[1].containsMarkerMutation(m3, 5999999)) {')
print('firstY = fMG[1];')
print('firstX = fMG[0];')
print('} else')
print('stop("### ERROR: no Ys in first male");')
print('ymuts = firstY.mutationsOfType(m2);')
print('ycounts = sim.mutationCounts(NULL, ymuts);')
print('removeY = ymuts[ycounts == numY];')
print('xmuts = firstX.mutationsOfType(m2);')
print('xcounts = sim.mutationCounts(NULL, xmuts);')
print('removeX = xmuts[xcounts == numX];')
print('removes = c(removeY, removeX);')
print('sim.subpopulations.genomes.removeMutations(removes, T);')
print('}')
print('}')
print('6000:10000 late() {')
print('if (sim.generation % 100 == 0) {')
print('numMales =  sum(p1.individuals.sex == "M");')
print('maleInds = p1.individuals[p1.individuals.sex == "M"][(numMales-1)];')
print('Malecounts = sum(maleInds.genomes.countOfMutationsOfType(m4)); ')
print('maleFitness = 1+((Malecounts * '+se+')/numMales); ')
print('numFems =  sum(p1.individuals.sex == "F");')
print('femInds = p1.individuals[p1.individuals.sex == "F"][(numFems-1)];')
print('Femcounts = sum(femInds.genomes.countOfMutationsOfType(m4)); ')
print('femaleFitness = 1+((Femcounts * '+se+')/numFems); ')
print('cat("ratio: " + (femaleFitness/maleFitness) + "\\n");') 
print('}')
print('}')