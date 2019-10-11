# reads2poly
In this repository are a bunch of scripts to go between file formats commonly used in genetics for specific analyses. Most of the software referenced in bash scripts is not mine, and I will add links to the makers website.

alignTemplate.sh is a bash script template for aligning reads to an assembly and getting SNPs output.

fasta_maker.sh takes a vcf and makes a fasta files (using a script from Wei Wang)
this script also runs polymorphurama_interpop.pl which is an edited version of polymorphurama by Bachtrog and Andolfatto to do between population statistics. polymorphurama_axy.R can calculate whether there are significant differences in the statistics calculated above between the A X and Y chromosomes. This script uses BeginPerlBioinfoB_1.pm, which is a module for polymorphurama written by S. Wright.

dnds_maker.sh takes the output from fasta_maker to align the within population sequences to an outgroup sequence, in which the outgroup sequence may be too divergent to map properly, but that can still be found with BLAST. It uses codoner.pl to make sure there are enough bp to be in codons, and also uses PRANK's codon model for alignment.
abba_maker.sh takes the output from dnds_maker to cat loci into sequences for analyses like abba-baba tests. It uses the python script fasta_cat.py

seer.pl and summer.pl calculate statistics for the output of ld.sh (forthcoming)
