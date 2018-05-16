#https://www.rdocumentation.org/packages/evobiR/versions/1.1/topics/CalcD

library(evobiR)

args = commandArgs(trailingOnly=TRUE) #<- for loop over several loci 
#CalcPopD(alignment = args, ambig="D", align.format='fasta') <- CalcPopD will not run in this format

CalcD(alignment = args, sig.test = "J", block.size = 1000, replicate = 1000)

#CalcPopD(alignment = args)

#" we decided that using sites that weren't fixed didn't really add much " H. Blackmon