#https://www.rdocumentation.org/packages/evobiR/versions/1.1/topics/CalcD

library(evobiR)

args = commandArgs(trailingOnly=TRUE) #<- for loop over several loci 
#CalcPopD(alignment = args, ambig="D", align.format='fasta') <- CalcPopD will not run in this format

CalcPopD(alignment = args)
