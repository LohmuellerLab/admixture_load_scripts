#! /usr/bin/python

import sys

infilename = sys.argv[1]

genomesize = 29100000
nummarkers = genomesize/1000 + 1

sum_not_int = 0
sum_int_derived = 0

with open(infilename, 'r') as lines:
    for line in lines:
        if not line.startswith('#'):
            chrom,pos,id,ref,alt,qual,filter,info,format = line.strip('\n').split('\t')[0:9]
            if not info.split(';')[5] == 'MT=10':
                pass
            else:
                gts = line.strip('\n').split('\t')[9:]
                gts = [x.replace('|','') for x in gts]
                N = len(gts)
                gts = ''.join(gts)
                #need to count only 1s because sites with all 0s are lost
                if 'p1' in infilename:
                    not_int = gts.count('1')
                    int_derived = (nummarkers * N * 2) - not_int
                elif 'p2' in infilename:
                    int_derived = gts.count('1')
                    not_int = (nummarkers * N * 2) - int_derived
                sum_not_int = sum_not_int + not_int
                sum_int_derived = sum_int_derived + int_derived

print sum_int_derived/(nummarkers * N * 2.)
