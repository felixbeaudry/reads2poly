import fileinput
import argparse

def arguments():
        parser  = argparse.ArgumentParser(description="import .csv into SQL")
        parser.add_argument("-i","--input",help="path for input",required=True)
        parser.add_argument("-c","--column",help="column name",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()
inpath = args.input
infile = open(inpath,'r')
col = args.column

for line in infile:
        line_without_newline = line.splitlines()[0]
        values = line_without_newline.split('\t')
        [Locus,
        length,
        mapped,
        unmapped
        ] = values


        output = \
        """
        INSERT INTO %s VALUES ('%s', '%s');
        """ %(col,
                Locus,
                mapped
                )
        print output