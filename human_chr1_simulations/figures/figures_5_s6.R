setwd('~/Dropbox/project_admixture_load/fixed_chr')
library(ggplot2)
library(reshape2)

#multi plotting function from R cookbook: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }


 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

condense_data <- function(data, windowsize=100000,end=100000000){
	startposition = 1
	new.data <- data[0,]
	while ((startposition + windowsize - 1) <= end){
		data.s <- data[(data$start >= startposition) & (data$start <= startposition + windowsize - 1),]
		endposition = as.numeric(data.s[10,][3])
		data.s <- cbind(cbind(data.s[1,][1:2],data.s[10,][3]),data.frame(t(colMeans(data.s)[4:length(colMeans(data.s))])))
		new.data <- rbind(new.data,data.s)
		startposition = startposition + windowsize
	}
	
	return(new.data)
}

#genome size
maxcoord = 100000000
#plotting window size
intervalsize = 100000
#number of simulation replicates
numreps=100

###Parse all simulation data into single data frame
#make empty data frame
header = c("h","start","end","recRate","percentExon", sprintf("X%s",seq(1,numreps,by=1)), "mean", "var")

#model 0
model="0"

data_additive = read.csv(sprintf('summary_h0.5_model%s.csv',model))[header]
data_additive = condense_data(data_additive, windowsize=intervalsize, end=maxcoord)
data_additive = data_additive[data_additive$start < maxcoord,]

data_recessive = read.csv(sprintf('summary_h0.0_model%s.csv',model))[header]
data_recessive = condense_data(data_recessive, windowsize=intervalsize, end=maxcoord)
data_recessive = data_recessive[data_recessive$start < maxcoord,]

data_hs = read.csv(sprintf('summary_hs_model%s.csv',model))[header]
data_hs$h <- c(0)
data_hs = condense_data(data_hs, windowsize=intervalsize, end=maxcoord)
data_hs = data_hs[data_hs$start < maxcoord,]
data_hs$h <- c("h(s)")

data_model0 = rbind(data_additive,data_recessive,data_hs)
data_model0$model = c("Model 0")

#model 2
model="2"

data_additive = read.csv(sprintf('summary_h0.5_model%s.csv',model))[header]
data_additive = condense_data(data_additive, windowsize=intervalsize, end=maxcoord)
data_additive = data_additive[data_additive$start < maxcoord,]

data_recessive = read.csv(sprintf('summary_h0.0_model%s.csv',model))[header]
data_recessive = condense_data(data_recessive, windowsize=intervalsize, end=maxcoord)
data_recessive = data_recessive[data_recessive$start < maxcoord,]

data_hs = read.csv(sprintf('summary_hs_model%s.csv',model))[header]
data_hs$h <- c(0)
data_hs = condense_data(data_hs, windowsize=intervalsize, end=maxcoord)
data_hs = data_hs[data_hs$start < maxcoord,]
data_hs$h <- c("h(s)")

data_model2 = rbind(data_additive,data_recessive,data_hs)
data_model2$model = c("Model 2")

#model 4

model="4"

data_additive = read.csv(sprintf('summary_h0.5_model%s.csv',model))[header]
data_additive = condense_data(data_additive, windowsize=intervalsize, end=maxcoord)
data_additive = data_additive[data_additive$start < maxcoord,]

data_recessive = read.csv(sprintf('summary_h0.0_model%s.csv',model))[header]
data_recessive = condense_data(data_recessive, windowsize=intervalsize, end=maxcoord)
data_recessive = data_recessive[data_recessive$start < maxcoord,]

data_hs = read.csv(sprintf('summary_hs_model%s.csv',model))[header]
data_hs$h <- c(0)
data_hs = condense_data(data_hs, windowsize=intervalsize, end=maxcoord)
data_hs = data_hs[data_hs$start < maxcoord,]
data_hs$h <- c("h(s)")

data_model4 = rbind(data_additive,data_recessive,data_hs)
data_model4$model = c("Model 4")

###figure 5
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)

exons = read.csv('genes_merged_shifted.csv')
n_ex = length(exons[,1])
exons = exons[exons$start < maxcoord,]

#ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.1) + geom_point(aes(x=data_model0$start+50000,y=data_model0$mean,colour=factor(data_model0$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=2) + geom_smooth(aes(x=data_model0$start+50000,y=data_model0$mean,colour=factor(data_model0$h)),method="loess",span=0.01, se=FALSE) + scale_y_continuous(" ") + scale_x_continuous("") + theme_bw() + theme(legend.position="None")

cols = c("red","#FF9933", "blue") #00000000
scaleFUN <- function(x) sprintf("%.3f", x)

ptop <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.1) + geom_point(aes(x=data_model0$start+50000,y=data_model0$mean,colour=factor(data_model0$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=2) + geom_smooth(aes(x=data_model0$start+50000,y=data_model0$mean,colour=factor(data_model0$h)),method="loess",span=0.01, se=FALSE) + scale_y_continuous(" ", labels=scaleFUN) + scale_x_continuous("") + theme_bw() + theme(legend.position=c(0.975,0.975),legend.justification=c(1,1),legend.background = element_rect(fill=alpha('grey', 0.4))) + scale_colour_manual(expression(paste("Dominance coefficient, (",italic("h"),")")),values=cols) 

pmid <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.1) + geom_point(aes(x=data_model2$start+50000,y=data_model2$mean,colour=factor(data_model2$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=2) + geom_smooth(aes(x=data_model2$start+50000,y=data_model2$mean,colour=factor(data_model2$h)),method="loess",span=0.01, se=FALSE) + scale_y_continuous("", labels=scaleFUN) + scale_x_continuous("") + theme_bw() + scale_colour_manual(values=cols) + guides(colour=FALSE)

#expression(paste(italic("p"["I"])," 10,000 generations after admixture"))

pbot <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.1) + geom_point(aes(x=data_model4$start+50000,y=data_model4$mean,colour=factor(data_model4$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=2) + geom_smooth(aes(x=data_model4$start+50000,y=data_model4$mean,colour=factor(data_model4$h)),method="loess",span=0.01, se=FALSE) + scale_y_continuous(" ", labels=scaleFUN) + scale_x_continuous("Genomic position") + theme_bw() + scale_colour_manual(values=cols) + guides(colour=FALSE)

grid.newpage()
pdf("figure_5.pdf",width=10,height=10)
lg <- tableGrob(c("Model 0","Model 2","Model 4"))
rg <- arrangeGrob(ptop,pmid,pbot,ncol=1,nrow=3)
grid.draw(cbind(lg,rg,size="last"))
dev.off()

q <- ggplot() + geom_point(aes(x=data_model0$start+50000,y=data_model2$recRate),size=0.2,colour="grey") + geom_smooth(aes(x=data_model2$start+50000,y=data_model2$recRate),method="loess",span=0.01, se=FALSE, colour="black") + scale_y_continuous("") + scale_x_continuous("") + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text.x=element_blank(),panel.border=element_blank())+ scale_colour_manual(values=c("grey")) + guides(colour=FALSE)
ggsave("figure_5_recomb.pdf",width=10,height=1.5,units="in")

#figure s6
which="X1"

ptop <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.1) + geom_point(aes(x=data_model0$start+50000,y=data_model0[[which]],colour=factor(data_model0$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=1, alpha=0.75) + geom_smooth(aes(x=data_model0$start+50000,y=data_model0[[which]],colour=factor(data_model0$h)),method="loess",span=0.01, se=FALSE) + scale_y_continuous(" ") + scale_x_continuous("") + theme_bw() + theme(legend.position=c(0.975,0.975),legend.justification=c(1,1), legend.background = element_rect(fill=alpha('grey', 0.4))) + scale_colour_manual(values=cols) + scale_colour_manual(expression(paste("Dominance coefficient, (",italic("h"),")")),values=cols) + geom_hline(yintercept=mean(data_model0[data_model0$h == "0",]$mean), colour=cols[1], linetype=2) + geom_hline(yintercept=mean(data_model0[data_model0$h == "0.5",]$mean), colour=cols[2], linetype=2) + geom_hline(yintercept=mean(data_model0[data_model0$h == "hs",]$mean), colour=cols[3], linetype=2)

pmid <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.1) + geom_point(aes(x=data_model2$start+50000,y=data_model2$X1,colour=factor(data_model2$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=1, alpha=0.75) + geom_smooth(aes(x=data_model2$start+50000,y=data_model2$X1,colour=factor(data_model2$h)),method="loess",span=0.01, se=FALSE) + scale_y_continuous(" ") + scale_x_continuous("") + theme_bw() + scale_colour_manual(values=cols) + guides(colour=FALSE) + geom_hline(yintercept=mean(data_model2[data_model2$h == "0",]$mean), colour=cols[1], linetype=2) + geom_hline(yintercept=mean(data_model2[data_model2$h == "0.5",]$mean), colour=cols[2], linetype=2) + geom_hline(yintercept=mean(data_model2[data_model2$h == "hs",]$mean), colour=cols[3], linetype=2)

pbot <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.1) + geom_point(aes(x=data_model4$start+50000,y=data_model4$X1,colour=factor(data_model4$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=1, alpha=0.75) + geom_smooth(aes(x=data_model4$start+50000,y=data_model4$X1,colour=factor(data_model4$h)),method="loess",span=0.01, se=FALSE) + scale_y_continuous(" ") + scale_x_continuous("Genomic position") + theme_bw() + scale_colour_manual(values=cols) + guides(colour=FALSE) + geom_hline(yintercept=mean(data_model4[data_model4$h == "0",]$mean), colour=cols[1], linetype=2) + geom_hline(yintercept=mean(data_model4[data_model4$h == "0.5",]$mean), colour=cols[2], linetype=2) + geom_hline(yintercept=mean(data_model4[data_model4$h == "hs",]$mean), colour=cols[3], linetype=2)

grid.newpage()
pdf(sprintf("figure_s6_%s.pdf",which),width=10,height=10)
lg <- tableGrob(c("Model 0","Model 2","Model 4"))
rg <- arrangeGrob(ptop,pmid,pbot,ncol=1,nrow=3)
grid.draw(cbind(lg,rg,size="last"))
dev.off()