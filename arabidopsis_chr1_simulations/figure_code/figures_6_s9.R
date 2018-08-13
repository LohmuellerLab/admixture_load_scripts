setwd('~/Dropbox/project_admixture_load/arabidopsis_sims')
#setwd('/mnt/e/Dropbox/project_admixture_load/fixed_chr')

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


maxcoord = 29100000
intervalsize = 100000
repnum = 100

#model 0
model = "0"

data_hs = read.csv(sprintf('summary_hs_model%s.csv',model))
data_hs$h <- c(0)
data_hs = condense_data(data_hs, windowsize=intervalsize, end=maxcoord)
data_hs = data_hs[data_hs$start < maxcoord,]
data_hs$h <- c("hs")
data_hs$model <- c(sprintf("model%s",model))

data_model0 = data_hs

#model 2

model = "2"

data_hs = read.csv(sprintf('summary_hs_model%s.csv',model))
data_hs$h <- c(0)
data_hs = condense_data(data_hs, windowsize=intervalsize, end=maxcoord)
data_hs = data_hs[data_hs$start < maxcoord,]
data_hs$h <- c("hs")
data_hs$model <- c(sprintf("model%s",model))

data_model2 = data_hs


#model 4

model = "4"

data_hs = read.csv(sprintf('summary_hs_model%s.csv',model))
data_hs$h <- c(0)
data_hs = condense_data(data_hs, windowsize=intervalsize, end=maxcoord)
data_hs = data_hs[data_hs$start < maxcoord,]
data_hs$h <- c("hs")
data_hs$model <- c(sprintf("model%s",model))

data_model4 = data_hs

###figure 6

library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)

exons = read.csv('exons.csv')
n_ex = length(exons[,1])
exons = exons[exons$start < maxcoord,]

cols = c("blue") #00000000
repid = "mean"

ptop <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.01) + geom_point(aes(x=data_model0$start+50000,y=data_model0$mean,colour=factor(data_model0$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=2) + geom_smooth(aes(x=data_model0$start+50000,y=data_model0[[repid]],colour=factor(data_model0$h)),method="loess",span=0.05, se=FALSE) + scale_y_continuous(" ") + scale_x_continuous("") + theme_bw() + theme(legend.position=c(0.975,0.5),legend.justification=c(1,1),legend.background = element_rect(fill=alpha('grey', 0.4))) + scale_colour_manual(values=cols) + scale_colour_manual(expression(paste("Dominance coefficient, (",italic("h"),")")),values=cols)

pmid <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.01) + geom_point(aes(x=data_model2$start+50000,y=data_model2[[repid]],colour=factor(data_model2$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=2) + scale_y_continuous(" ") + scale_x_continuous("") + theme_bw() + scale_colour_manual(values=cols) + guides(colour=FALSE)

pbot <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.01) + geom_point(aes(x=data_model4$start+50000,y=data_model4[[repid]],colour=factor(data_model4$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=2) + scale_y_continuous("") + scale_x_continuous("Genomic position") + theme_bw() + scale_colour_manual(values=cols) + guides(colour=FALSE)

grid.newpage()
pdf("figure_6.pdf",width=10,height=10)
lg <- tableGrob(c("Model 0","Model 2","Model 4"))
rg <- arrangeGrob(ptop,pmid,pbot,ncol=1,nrow=3)
grid.draw(cbind(lg,rg,size="last"))
dev.off()

q <- ggplot() + geom_point(aes(x=data_model0$start+50000,y=data_model2$recRate),size=0.2,colour="grey") + geom_smooth(aes(x=data_model2$start+50000,y=data_model2$recRate),method="loess",span=0.01, se=FALSE, colour="black") + scale_y_continuous("") + scale_x_continuous("") + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text.x=element_blank(),panel.border=element_blank())+ scale_colour_manual(values=c("grey")) + guides(colour=FALSE)
ggsave("figure_6_recomb.pdf",width=10,height=1.5,units="in")

#figure s9

library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)

exons = read.csv('exons.csv')
n_ex = length(exons[,1])
exons = exons[exons$start < maxcoord,]

cols = c("blue") #00000000
repid="X4"

ptop <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.01) + geom_point(aes(x=data_model0$start+50000,y=data_model0[[repid]],colour=factor(data_model0$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=2) + geom_smooth(aes(x=data_model0$start+50000,y=data_model0[[repid]],colour=factor(data_model0$h)),method="loess",span=0.05, se=FALSE) + scale_y_continuous(" ") + scale_x_continuous("") + theme_bw() + theme(legend.position=c(0.975,0.5),legend.justification=c(1,1),legend.background = element_rect(fill=alpha('grey', 0.4))) + scale_colour_manual(values=cols) + scale_colour_manual(expression(paste("Dominance coefficient, (",italic("h"),")")),values=cols)

pmid <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.01) + geom_point(aes(x=data_model2$start+50000,y=data_model2[[repid]],colour=factor(data_model2$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=2) + scale_y_continuous(" ") + scale_x_continuous("") + theme_bw() + scale_colour_manual(values=cols) + guides(colour=FALSE)

pbot <- ggplot() + geom_rect(aes(xmin=exons$start,xmax=exons$stop,ymin=0,ymax=Inf),fill="blue",alpha=0.01) + geom_point(aes(x=data_model4$start+50000,y=data_model4[[repid]],colour=factor(data_model4$h)),size=0.2) + geom_hline(yintercept=0.05, colour="black", linetype=2) + scale_y_continuous("") + scale_x_continuous("Genomic position") + theme_bw() + scale_colour_manual(values=cols) + guides(colour=FALSE)

grid.newpage()
pdf("figure_s9.pdf",width=10,height=10)
lg <- tableGrob(c("Model 0","Model 2","Model 4"))
rg <- arrangeGrob(ptop,pmid,pbot,ncol=1,nrow=3)
grid.draw(cbind(lg,rg,size="last"))
dev.off()