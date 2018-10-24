#mu=0.7e-8
#ns:s = 2.31

#set to
setwd(dirname(parent.frame(2)$ofile))
options(scipen=999)

data = read.csv('chr1_rmap_cM_cleaned.csv')

samples=c("P2","P3","P6","P7","P8","P9","P10","P12","P15","P17","P19","P20","P35","P66","P129","P145","P169")

starts = as.numeric(sapply(samples, function(x) data[!is.na(data$start),][data[!is.na(data$start),]$sample == x,]$start[1]))
stops = as.numeric(sapply(samples, function(x) tail(data[!is.na(data$start),][data[!is.na(data$start),]$sample == x,]$stop,n=1)))

start = max(starts)
stop = min(stops)

data = data[!is.na(data$start),]
#data = data[data$start >= start,]
#data = data[data$stop <= stop,]

f <- function(d){(1/2) * (1-exp(-2*d/100))} #haldane's map function

options(warn=2)

rate.compute <- function(data,sample,coords){
    data_subset=data[data$population==sample,]
    #n = sample_counts[[popn]]
    #data_subset$rate = as.numeric(data_subset$XOs)/(n*(data_subset$stop - (data_subset$start + 1)))
    data_subset$rate = sapply(data_subset$H, f)/(data_subset$stop-data_subset$start)
    
    recRates = c()
    
    options(warn=2)
    for (i in 1:(length(coords)-1)){
        a = coords[i]; b=coords[i+1]
        chunk1 = data_subset[a >= data_subset$start & a < data_subset$stop,]
        chunk2 = data_subset[b >= data_subset$start & b < data_subset$stop,]
        if (chunk1$start == chunk2$start){
            rate = chunk1$rate
        } else if (chunk1$start != chunk2$start){
            overlap1=(chunk1$stop-a)/by
            overlap2=(b-chunk2$start)/by
            rate = (overlap1 * chunk1$rate) + (overlap2 * chunk2$rate)
        }
        
        recRates = append(recRates, rate)
    }
    options(warn=1)
    
    return(recRates)
}

by=100000;start=488426;stop=29595658
stop = stop - ((stop-start) %% by)
coords = seq(start,stop,by=by)

data.rates = data.frame(start=coords[1:length(coords)-1],stop=coords[2:length(coords)])
for (i in samples){
    data.rates[[i]] = rate.compute(data,i,coords)
}

data.rates$means = rowMeans(data.rates[samples])
data.rates$stopslim = data.rates$stop - 488426 # -1 sould be here
data.rates$lab = c("recRate")

chr.exon = read.csv('araport_exons_merged.bed',sep="\t",header=FALSE)
names(chr.exon) = c("chr","start","stop")
chr.neutral = read.csv('araport_exons_merged_complement.bed',sep="\t",header=FALSE)
names(chr.neutral) = c("chr","start","stop")
chr.exon$stop=chr.exon$stop-1
chr.neutral$stop=chr.neutral$stop-1

chr.exon$lab = "exon"
chr.neutral$lab = "neutral"

chr.exon = chr.exon[chr.exon$stop >= start,]
chr.neutral = chr.neutral[chr.neutral$stop >= start,]
chr.exon = chr.exon[chr.exon$start <= stop,]
chr.neutral = chr.neutral[chr.neutral$start <= stop,]

try(chr.exon[chr.exon$start < start,]$start <- start)
try(chr.neutral[chr.neutral$start < start,]$start <- start)

try(chr.exon[chr.exon$stop > stop,]$stop <- stop)
try(chr.neutral[chr.neutral$stop > stop,]$stop <- stop)

chr.exon$start = chr.exon$start - start; chr.neutral$start = chr.neutral$start - start
chr.exon$stop = chr.exon$stop - start; chr.neutral$stop = chr.neutral$stop - start

write.table(chr.exon[c("lab","start","stop")], file="sim_seq_info.txt", append=FALSE, quote=FALSE, sep=" ",row.names=FALSE,col.names=FALSE)
write.table(chr.neutral[c("lab","start","stop")], file="sim_seq_info.txt", append=TRUE, quote=FALSE, sep=" ",row.names=FALSE,col.names=FALSE)
write.table(data.rates[c("lab","stopslim","means")], file="sim_seq_info.txt", append=TRUE, quote=FALSE, sep=" ",row.names=FALSE,col.names=FALSE)