import fileinput
import argparse
import sys
import re

def arguments():
        parser  = argparse.ArgumentParser(description="import .csv into SQL")
        parser.add_argument("-i","--input",help="file containing names of individuals to merge",required=True)
        parser.add_argument("-l","--lociNames",help="name of SQL database containing complete list of loci",required=True)
        parser.add_argument("-t","--tableName",help="name of SQL database into which to merge",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
inpath = args.input
infile = open(inpath,'r')
loc = args.lociNames
table = args.tableName

sys.stdout.write('echo "CREATE TABLE '+table+' AS (SELECT '+loc+'.*')
                
for line in infile:
        line = line.rstrip()
        line = re.sub('\.', '', line)
        sys.stdout.write(', '+line+'.mapped AS '+line)

sys.stdout.write(' FROM '+loc)

infile = open(inpath,'r')
for line in infile:
        line = line.rstrip()
        line = re.sub('\.', '', line)
        sys.stdout.write(' LEFT JOIN '+line+' ON '+loc+'.Loci = '+line+'.Locus') 

sys.stdout.write(');" | mysql -u felix.beaudry -D felix_rha')



