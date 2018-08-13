"""
usage: python blocklengthdistribution.py <infilename>
"""

import sys, pandas, numpy

#use pandas data frame to record the deserts of introgressed ancestry
#df = pandas.DataFrame(data={'start':starts, 'end':ends, 'ancestry':p1ancestry})

infilename = sys.argv[1]
h = infilename.strip('.ancestry').split('_')[0]
rep = infilename.strip('.ancestry').split('_')[2]
df = pandas.read_csv(infilename)

#assign numbers to contiguous same ancestry blocks
df['block'] = (df.ancestry.shift(1) != df.ancestry).astype(int).cumsum()

block_pis=[]
block_lengths = []

#go through each unique block and get length
for ii in list(df.block.unique()):
    df_subset = df.loc[df['block'] == ii]
    starts_subset = numpy.array(df_subset.start)
    ends_subset = numpy.array(df_subset.end)
    block_length = numpy.sum(ends_subset - starts_subset)
    block_lengths.append(block_length)
    block_pis.append(df_subset.ancestry.unique()[0])
    #desert_lengths = sorted(desert_lengths, reverse=True)

outfilename = infilename.strip('.ancestry') + '_blockdist.csv'
outfile = open(outfilename, 'w')
outfile.write('h,rep,pi,length\n')

for pi, length in zip(block_pis, block_lengths):
    outfile.write('{0},{1},{2},{3}\n'.format(h,rep,pi,length))

outfile.close()