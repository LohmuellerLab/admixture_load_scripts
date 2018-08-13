setwd('~/Dropbox/project_admixture_load/sims_final/vary_splittime_new')

library(reshape2)
library(ggplot2)

###
### code for Figure 5
###

# set scaling factor to what was used in simulations
scalingfactor=5

#prepare and build data frames
data = read.csv('model0_split.csv')
splittimes = unique(data$splittime)
data = data[0,]
data$model = character()
data_fst = data[0,]
data_fst$model = character()

for (model in c("0","4")){
	data_m = read.csv(sprintf('./model%s_split.csv',model))
    for (splittime in splittimes){
        data_i = data_m[data_m$splittime == as.integer(splittime),]
        if (model == "0"){
            data_i$model = rep("0: No size change", length(data_i[[1]]))
        } else if (model == "1"){
            data_i$model = rep("1: Short bottleneck", length(data_i[[1]]))
        } else if (model == "4"){
            data_i$model = rep("4: Long bottleneck", length(data_i[[1]]))
        }
        data = rbind(data, data_i[data_i$gen == (100000/scalingfactor) + as.integer(splittime)/scalingfactor + (10000/scalingfactor),])
        data_fst = rbind(data_fst, data_i[data_i$gen == ((100000 + splittime)/scalingfactor - 1),])
    }
}

#make figure
library(grid)
p<-ggplot(data, aes(x=factor(splittime), fill=factor(model), colour=factor(model), p2Fraction)) + geom_violin() + stat_summary(fun.y=mean, fun.ymin=function(x) mean(x) - sd(x), fun.ymax = function(x) mean(x) + sd(x), geom="pointrange", position=position_dodge(width=.9), colour="black") + geom_hline(yintercept=0.05, linetype="dashed") + scale_y_continuous("Proportion of ancestry from introgression") + scale_x_discrete("Time between split and admixture (generations)") + theme_bw(base_size=12) + scale_fill_manual("Model",values=c("blue","red")) + scale_colour_manual(values=c("blue","red")) + guides(colour=FALSE) + theme(plot.margin=margin(t=30,l=5,b=5,r=5))


#get average Fst for the simulation generation immediately before the admixture event
library(plyr)
dt <- data_fst[c("generation","FST_bhatia","model")]
gen = ddply(dt[dt$model == "0: No size change",],~generation,summarise,mean=mean(FST_bhatia))$generation
mean0 = ddply(dt[dt$model == "0: No size change",],~generation,summarise,mean=mean(FST_bhatia))$mean
mean4 = ddply(dt[dt$model == "4: Long bottleneck",],~generation,summarise,mean=mean(FST_bhatia))$mean
dt = data.frame(splittimes, mean0, mean4)
dt.m <- melt(dt, id.vars=c("splittimes"))

xpos = (1:13)-0.3
mean0 = round(dt$mean0, digits=3)
mean4 = round(dt$mean4, digits=3)

#add annotations and numbers to figure
p <- p + annotation_custom(
    grob = textGrob(label=expression(paste(italic("F"["ST"]))), hjust=0, gp=gpar(cex=1.5, fontsize=10, col="black")),
    ymin = 1.09,
    ymax = 1.09,
    xmin = -0.3,
    xmax = -0.3
)

for (i in 1:13){
p <- p + annotation_custom(
    grob = textGrob(label=mean0[i], hjust=0, gp=gpar(cex=1.5, fontsize=6, col="blue")),
    ymin = 1.125,
    ymax = 1.125,
    xmin = xpos[i],
    xmax = xpos[i]
)
}

for (i in 1:13){
p <- p + annotation_custom(
    grob = textGrob(label=mean4[i], hjust=0, gp=gpar(cex=1.5, fontsize=6, col="red")),
    ymin = 1.07,
    ymax = 1.07,
    xmin = xpos[i],
    xmax = xpos[i]
)
}

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

#save figure as pdf
ggsave('Figure_5.pdf',plot=gt,height=10, width=20, units="cm")
