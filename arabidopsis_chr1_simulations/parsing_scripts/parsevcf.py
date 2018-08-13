#! /usr/bin/env python

import sys

def main():
    infilename = sys.argv[1]
    header = True
    outfilename = infilename.strip('.vcf') + '_summary.txt'
    infile = open(infilename, 'r')
    lines = infile.readlines()
    outfile = open(outfilename, 'w')
    outfile.write('position,anc,der\n')

    for line in lines:
        if not line.startswith('#'):
            line = line.strip('\n')
            position = float(line.split('\t')[1])
            gts = line.split('\t')[9:]
            gts = "|".join(gts)
            anc = gts.count("0")
            der = gts.count("1")
            outfile.write('{0:.0f},{1:.0f},{2:.0f}\n'.format(position,anc,der))
    
    outfile.close()
    infile.close()

if __name__ == "__main__":
    main()
