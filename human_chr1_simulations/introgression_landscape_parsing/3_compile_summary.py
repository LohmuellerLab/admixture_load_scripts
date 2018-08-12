#! /usr/bin/env python

from decimal import *
getcontext().prec = 6
import sys

def variance(data):
    n = Decimal(len(data))
    mean = sum(data)/n
    return (1/(n-1))*sum([(x-mean)**2 for x in data])

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)

def main():
    h=sys.argv[1]
    repNum = 100
    simFolder = "admixture_simulations_fullsims"
        
    #open output file
    outfilename = 'summary.csv'
    outfile = open(outfilename ,'w')
    headerline = 'h,start,end,recRate,percentExon,' + ','.join(map(str,range(1,repNum+1))) + ',mean,var\n'
    outfile.write(headerline)
    
    if True:
        #open recombination file
        recomblist = open('/u/home/b/bkim331/project-klohmueldata/bernard_data/{0}/recomb.txt'.format(simFolder),'r').readlines()
        recomblist = [[int(float(x.strip('\n').split(' ')[1]))-9999,int(float(x.strip('\n').split(' ')[1])),Decimal(x.strip('\n').split(' ')[2])] for x in recomblist]

        #open exon files
        exonlist = open('/u/home/b/bkim331/project-klohmueldata/bernard_data/{0}/exons.csv'.format(simFolder),'r').readlines()
        exonlist = [[int(y) for y in x.strip('\n').split(',')] for x in exonlist[1:]]
        
        #open ancestry % files as array
        infilenames = ['h{0}_full_{1}_p1.ancestryblocked'.format(h,i) for i in range(1,repNum+1)]
        lineslist = [open(x, 'r').readlines() for x in infilenames]
        lineslist = [x[1:] for x in lineslist]
        lineslist = [[y.strip('\n').split(',') for y in x] for x in lineslist]
        perIntrogressedlist = [map(Decimal, list(zip(*x))[2]) for x in lineslist]
        perIntrogressedlist = zip(*perIntrogressedlist)
        
        #iterate through intervals defined by recomb file
        for chunk, perIntrogressed in zip(recomblist, perIntrogressedlist):
            start, stop, recRate = chunk
            exonoverlap = sum([getOverlap([start,stop],x) for x in exonlist])/Decimal(stop-start+1) #could optimize better
            mean = sum(perIntrogressed)/repNum
            var = variance(perIntrogressed)
            outline = '{0},{1},{2},{3},{4},{5},{6},{7}\n'.format(h,start,stop,recRate,exonoverlap,','.join(map(str,perIntrogressed)),mean,var)
            outfile.write(outline)
        
        outfile.close()

if __name__ == "__main__":
    main()
