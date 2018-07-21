"""
usage: parseancestrytrees.py treefilename.trees
"""

import sys, msprime, pyslim
import pandas, numpy, random
#import matplotlib.pyplot

def ancestry2(ts):
    n=2000.
    #time at which population info was recorded
    mixtime=2001.
    #get nodes from p1 at mixtime-recorded nodes 
    p1 = [x.id for x in ts.nodes() if ((x.population == 1) and (x.time == mixtime))] 
    # get nodes alive today
    today = [x.id for x in ts.nodes() if x.time == 0.0]
    # how many currently-alive nodes descend from pop 1 in each tree
    tree_p1 = [sum([t.num_tracked_samples(u) for u in p1])/n
               for t in ts.trees(tracked_samples=today, sample_counts=True)]
    return tree_p1

infilename = sys.argv[1]
#infilename = 'h0.0_full_1_p1.trees'

# Load the .trees file and assess true local ancestry
ts = pyslim.load(infilename, slim_format=True)
p1ancestry = ancestry2(ts)
starts=[]
ends=[]

for x in ts.trees():
    starts.append(x.interval[0])
    ends.append(x.interval[1])

#tree = ts.first()
#print(tree.draw(format="unicode"))

#x = [x for pair in zip(starts, ends) for x in pair]
#y = [x for x in p1ancestry for _ in (0, 1)]
#matplotlib.pyplot.plot(x, y)
#matplotlib.pyplot.show()

outfilename = infilename.strip('.trees') + '.ancestry'

outfile = open(outfilename, 'w')
outfile.write('start,end,ancestry\n')

for start, end, anc in zip(starts, ends, p1ancestry):
    outfile.write('{0},{1},{2}\n'.format(start, end, anc))

outfile.close()

df = pandas.DataFrame(data={'start':starts, 'end':ends, 'ancestry':p1ancestry})
df['block'] = (df.ancestry.shift(1) != df.ancestry).astype(int).cumsum()

desert_lengths = []

for ii in list(df.block.unique()):
    df_subset = df.loc[df['block'] == ii]
    starts_subset = numpy.array(df_subset.start)
    ends_subset = numpy.array(df_subset.end)
    block_length = numpy.sum(ends_subset - starts_subset)
    desert_lengths.append(block_length)
    desert_lengths = sorted(desert_lengths, reverse=True)

outfilename = infilename.strip('.trees') + '.deserts'
outfile = open(outfilename, 'w')
outfile.write('length\n')

for length in desert_lengths:
    outfile.write('{0}\n'.format(length))

outfile.close()
