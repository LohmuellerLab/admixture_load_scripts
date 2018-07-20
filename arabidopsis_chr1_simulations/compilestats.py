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

def getDerivedCounts(start, stop, data, markerIntervalSize):
    numSites = (stop-start+1)/markerIntervalSize
    if (start == 1):
        numSites += 1
    totalVars = numSites * 2 * 1000 #1000 is simulated population size
    totalDers = sum([der*((pos >= start) & (pos <= stop)) for pos,anc,der in data]) #could be optimized
    return Decimal(totalDers)/Decimal(totalVars)

def main():
    h = sys.argv[1]
    repNum = int(sys.argv[2])
    markerIntervalSize = 1000
    simFolder = "admixture_simulations_fullsims"
        
    #open output file
    outfilename = 'summary.csv'
    outfile = open(outfilename ,'w')
    headerline = 'h,start,end,recRate,percentExon,' + ','.join(map(str,range(1,repNum+1))) + ',mean,var\n'
    outfile.write(headerline)
    
    if True:
        #open ancestry % files as array
        infilenames = ['h{0}_full_{1}_p2_markers_summary.txt'.format(h,i) for i in range(1,repNum+1)]
        lineslist = [open(x, 'r').readlines() for x in infilenames]
        lineslist = [x[1:] for x in lineslist]
        lineslist = [[map(int,x.strip('\n').split(',')) for x in y] for y in lineslist]
        
        #open recombination file
        recomblist = open('/u/home/b/bkim331/project-klohmueldata/bernard_data/{0}/recomb.txt'.format(simFolder),'r').readlines()
        recomblist = [[int(float(x.strip('\n').split(' ')[1]))-9999,int(float(x.strip('\n').split(' ')[1])),Decimal(x.strip('\n').split(' ')[2])] for x in recomblist]

        #open exon files
        exonlist = open('/u/home/b/bkim331/project-klohmueldata/bernard_data/{0}/exons.csv'.format(simFolder),'r').readlines()
        exonlist = [map(int,x.strip('\n').split(',')) for x in exonlist[1:]]
        
        #iterate through intervals defined by recomb file
        for chunk in recomblist:
            start, stop, recRate = chunk
            exonoverlap = sum([getOverlap([start,stop],x) for x in exonlist])/Decimal(stop-start+1) #could optimize better
            perIntrogressed = [getDerivedCounts(start, stop, x, markerIntervalSize) for x in lineslist]
            mean = sum(perIntrogressed)/repNum
            var = variance(perIntrogressed)
            outline = '{0},{1},{2},{3},{4},{5},{6},{7}\n'.format(h,start,stop,recRate,exonoverlap,','.join(map(str,perIntrogressed)),mean,var)
            outfile.write(outline)
        
        outfile.close()

if __name__ == "__main__":
    main()
