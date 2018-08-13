setwd('~/Dropbox/project_admixture_load/sims_final/vary_splittime_new')

library(reshape2)
library(ggplot2)

# line plot

data = read.csv('model0_split.csv')
data = data[0,]
data$model = character()

for (model in c("0","4")){
  for (splittime in c("4000")){
        data_m = read.csv(sprintf('./model%s_split.csv',model))
        data_i = data_m[data_m$splittime == as.integer(splittime),]
        if (model == "0"){
            data_i$model = rep("0: No size change", length(data_i[[1]]))
        } else if (model == "1"){
            data_i$model = rep("1: Short bottleneck", length(data_i[[1]]))
        } else if (model == "4"){
            data_i$model = rep("4: Long bottleneck", length(data_i[[1]]))
        }
        data = rbind(data, data_i)
  }
}

data$generation = (data$generation-10000)*10

data$privateToP1 = (data$totunique - data$p2unique)
data$privateToP2 = (data$totunique - data$p1unique)
data$shared = (data$p1unique + data$p2unique - data$totunique)
data.melted = data[c("generation","model","privateToP1","privateToP2","shared")]
data.melted = melt(data.melted, id.vars=c("generation","model"))
data.melted$variable <- factor(data.melted$variable)
levels(data.melted$variable) <- c("private to source", "private to recipient", "shared")

data.melted = data.melted[data.melted$generation != 40000,]

q2 <- ggplot(data.melted, aes(x=generation, y=value, linetype=factor(variable), colour=factor(variable))) + facet_grid(model~., scales="free_x") + stat_summary(fun.y=mean, geom="line", alpha=0.75) + scale_x_continuous("Time after split (generations)",limits=c(0,40000)) + geom_vline(xintercept=c(100,250,500,1000,2500,5000,10000,15000,20000,25000,30000,35000,40000), linetype=2, colour="grey", alpha=0.75)+ theme_bw() + scale_colour_manual("", values=c("red","blue","black")) + scale_linetype_manual(values=c(1,2,1)) + guides(linetype=FALSE) + scale_y_continuous("Number of sites") 
ggsave('Figure_5b.pdf',plot=q2,height=10, width=20, units="cm")

#FST

data = read.csv('model0_split.csv')
data = data[0,]
data$model = character()
scalingfactor=5

for (model in c("0","4")){
        data_m = read.csv(sprintf('./model%s_split.csv',model))
        data_i = data_m[data_m$splittime == 40000,]
        data_i$generation = data_i$generation*scalingfactor - 100000
        data_i = data_i[data_i$generation < 40000,]
        if (model == "0"){
            data_i$model = rep("0: No size change", length(data_i[[1]]))
        } else if (model == "4"){
            data_i$model = rep("4: Long bottleneck", length(data_i[[1]]))
        }
        data = rbind(data, data_i)
}

library(reshape2)
library(ggplot2)
data.melted = data[c("generation","model","FST_bhatia")]
names(data.melted) = c("generation","model","FST")
data.melted = melt(data.melted, id.vars=c("generation","model"))
data.melted$variable <- factor(data.melted$variable)
levels(data.melted$variable) <- c("FST")

q2 <- ggplot(data.melted, aes(x=generation, y=value, linetype=factor(model), colour=factor(model))) + stat_summary(fun.y=mean, geom="line") + scale_x_continuous("Time after split (generations)",limits=c(0,40000)) + geom_vline(xintercept=c(100,250,500,1000,2500,5000,10000,15000,20000,25000,30000,35000,40000), linetype=2, colour="grey", alpha=0.75)+ theme_bw() + scale_colour_manual("", values=c("black","black","black")) + scale_linetype_manual("Model",values=c(1,2,1)) + scale_y_continuous(expression(paste(italic("F"["ST"]))))  + guides(colour=FALSE)
ggsave('Figure_5c.pdf',plot=q2,height=10, width=20, units="cm")