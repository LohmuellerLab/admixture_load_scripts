#! /usr/bin/env python

from decimal import *
getcontext().prec = 6
import sys
import numpy as np

def main():
    infilename = sys.argv[1]
    infile = open(infilename, 'r').readlines()
    outfilename = infilename.strip('.ancestry') + '.ancestryblocked'
    
    infile = [x.strip('\n').split(',') for x in infile[1:]]
    starts, ends, ancestry = list(zip(*infile))
    starts = [int(float(x)) for x in starts]
    starts[0] += 1
    ends = [int(float(x))-1 for x in ends]
    ancestry = [Decimal(x) for x in ancestry]
    simFolder = "admixture_simulations_fullsims"
    
    recomblist = open('/u/home/b/bkim331/project-klohmueldata/bernard_data/{0}/recomb.txt'.format(simFolder),'r').readlines()
    recomblist = [[int(float(x.strip('\n').split(' ')[1]))-9999,int(float(x.strip('\n').split(' ')[1])),Decimal(x.strip('\n').split(' ')[2])] for x in recomblist]
    recombcoords =  list(zip(*recomblist))[0:2]
    recombcoords = list(zip(recombcoords[0],recombcoords[1]))
    
    infilecoords = list(zip(*infile))[0]
    
    ancestries = []
    ii=0
    
    for l, r in recombcoords:
        print ii
        ii+=1        
        l_check = [(x >= l) and (x <= r) for x in starts]
        r_check = [(x >= l) and (x <= r) for x in ends]
        check = [i for i,(x,y) in enumerate(zip(l_check,r_check)) if (x or y)]
        print i
        if len(check) == 0:
            idx = [x >= l for x in starts].index(True) - 1
            starts_temp = [starts[idx]]
            ends_temp = [ends[idx]]
            ancestry_temp=ancestry[idx]
        else:
            starts_temp = list(np.array(starts)[check])
            ends_temp = list(np.array(ends)[check])
            ancestry_temp = list(np.array(ancestry)[check])
        if (starts_temp[0] < l):
            starts_temp[0] = l
        if (ends_temp[-1] > r):
            ends_temp[-1] = r
        weights = (np.subtract(ends_temp,starts_temp)+Decimal(1))/(r-l+1)
        
        a = sum(np.multiply(weights, ancestry_temp))
        ancestries.append(a)
        
        newstart_idx = [(x < r) for x in ends].index(False)
        starts = starts[newstart_idx:]
        ends = ends[newstart_idx:]
        ancestry = ancestry[newstart_idx:]
        
    outfile = open(outfilename, 'w')
    outline = 'start,end,ancestry\n'
    outfile.write(outline)
    
    for (s, e), a in zip(recombcoords, ancestries):
        outline = '{0},{1},{2}\n'.format(s,e,a)
        outfile.write(outline)
    
    outfile.close()

if __name__ == "__main__":
    main()
