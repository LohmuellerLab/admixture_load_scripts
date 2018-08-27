"""
usage: parseancestrytrees.py treefilename.trees
"""

import sys, msprime, pyslim
import pandas, numpy, random
#import matplotlib.pyplot

def ancestry2(ts):
    """
    This function takes a tree sequence object and parses out local ancestry exactly.
    It can do so because the populations of all nodes were recorded at a specific 
    generation after the population split but before admixture. Thus, we just need to 
    take the current-day nodes and see which ones track their ancestry back to those 
    recorded nodes. Thanks to Peter Ralph and Ben Haller for help making this 
    function run faster.
    """
    #sample size
    #note this should be 1500. for sex chromosome
    n=2000.
    #time at which population info was recorded
    #note that times can differ depending on when the tree sequences were recorded,
    #e.g. during an early() event
    mixtime=2001.
    #get nodes from p1 at mixtime-recorded nodes 
    p1 = [x.id for x in ts.nodes() if ((x.population == 1) and (x.time == mixtime))] 
    # get nodes alive today
    today = [x.id for x in ts.nodes() if x.time == 0.0]
    # how many currently-alive nodes descend from pop 1 in each tree
    tree_p1 = [sum([t.num_tracked_samples(u) for u in p1])/n
               for t in ts.trees(tracked_samples=today, sample_counts=True)]
    return tree_p1

#
infilename = sys.argv[1]

# Load the tree sequence file 
ts = pyslim.load(infilename, slim_format=True)

#get vector that describes the % of ancestry derived from p1
p1ancestry = ancestry2(ts)

#record positions
starts=[]
ends=[]

for x in ts.trees():
    starts.append(x.interval[0])
    ends.append(x.interval[1])

#write ancestry information to a file
outfilename = infilename.strip('.trees') + '.ancestry'
outfile = open(outfilename, 'w')
outfile.write('start,end,ancestry\n')

for start, end, anc in zip(starts, ends, p1ancestry):
    outfile.write('{0},{1},{2}\n'.format(start, end, anc))

outfile.close()