setwd('~/Dropbox/project_admixture_load/sims_final/')
library(ggplot2)
library(reshape2)

###figure s2

scalingfactor=5

stattitle = "Mean fitness"
xtitle = "Time after split (thousands of generations)"

df <- data.frame(r=character(), gen=character(), stat_lower=character(), stat_upper=character(), stat_mean=character(), stat_name=character(), sim=character(), h=character())

for (sim in c("model0","model1","model2","model3","model4")){
    
    for (h in c("0.0","0.5")){
      for (stat in c("meanFitnessP1","meanFitnessP2")){
        data = read.csv(sprintf('./h_%s_%s.csv',h,sim))
        data$fitnessRatio = data$meanFitnessP1/data$meanFitnessP2
        data_fit = data[c("rep","r","generation",stat)]
        
        generations = 2000:2600*5
        generations = generations*10/scalingfactor
        rs = c(1e-6, 1e-7, 1e-8, 1e-9)

        data_ribboned = data.frame(r=character(), gen=character(), stat_lower=character(), stat_median=character(), stat_upper=character(), stat_mean=character(), stat_name=character())
        data_names = names(data_ribboned)
        
        for (r in rs){
            for (i in generations){
                qqs = quantile(data_fit[data_fit$generation == i & data_fit$r == r,][c(stat)][[1]], probs=c(0.25,0.5,0.75),na.rm=TRUE)
                stat_lower = qqs[1]
                stat_median = qqs[2]
                stat_upper = qqs[3]
                stat_mean = data_fit[data_fit$generation == i & data_fit$r == r,][c(stat)][[1]]
                stat_mean[stat_mean <= 0] <- NA
                stat_mean = mean(stat_mean, na.rm=TRUE)	
                data_ribboned = rbind(data_ribboned,data.frame(r=r, gen=i, stat_lower, stat_median, stat_upper, stat_mean, stat))
            }
        }

        data_ribboned$gen <- (data_ribboned$gen - (100000/scalingfactor)) * scalingfactor
        rownames(data_ribboned) <- c()
        
        data_ribboned$h = rep(h, length(data_ribboned[,1]))
        data_ribboned$sim = rep(sim, length(data_ribboned[,1]))
        df = rbind(df, data_ribboned)
      }
    }
}

df$h <- factor(df$h)
levels(df$h) <- c("paste(italic(h),'=0.0')","paste(italic(h),'=0.5')")

df$sim <- factor(df$sim)
levels(df$sim) <- c("'Model 0'", "'Model 1'", "'Model 2'", "'Model 3'", "'Model 4'")

df$stat <- factor(df$stat)
levels(df$stat) <- c("italic(w[D])","italic(w[R])")


df.melted <- df
q <- ggplot(df.melted, aes(colour=factor(r))) + facet_grid(sim~h+stat, scales="free_y", labeller=label_parsed) + geom_line(aes(x=gen/1000,y=stat_mean),linetype=1) + geom_vline(xintercept=20000/1000,linetype=2,colour="grey",alpha=0.75) + scale_colour_discrete("recombination rate",labels=c(expression(10^"-9"),expression(10^"-8"),expression(10^"-7"),expression(10^"-6"))) + scale_fill_discrete("recombination rate",labels=c(expression(10^"-9"),expression(10^"-8"),expression(10^"-7"),expression(10^"-6"))) + scale_y_continuous(stattitle) + scale_x_continuous(xtitle) + theme_bw(base_size=14) + theme(axis.title.y=element_text(angle=90,vjust=0.5))
setwd('~/Dropbox/project_admixture_load/sims_final/')
ggsave("figure_S2.pdf",q,width=12, height=6, units="in")


######Fig S3

stat = "p2MeanDel"
stattitle = "Mean derived deleterious variants per haplotype"
xtitle = "Time after split (generations)"

df <- data.frame(r=character(),gen=character(), stat_mean=character(),sim=character(),h=character())

for (sim in c("model0","model1","model2","model3","model4")){
    for (h in c("0.0","0.5")){
        data = read.csv(sprintf('./h_%s_%s.csv',h,sim))
        data_fit = data[c("rep","r","generation",stat)]
        
        generations = 2000:2600*5
        generations = generations*10/scalingfactor
        rs = c(1e-6, 1e-7, 1e-8, 1e-9)
        
        data_ribboned = data.frame(r=character(), gen=character(), stat_mean=character())
        data_names = names(data_ribboned)
        
        for (r in rs){
            for (i in generations){
                stat_mean = data_fit[data_fit$generation == i & data_fit$r == r,][c(stat)][[1]]
            	stat_mean[stat_mean==0] <- NA
                stat_mean = mean(stat_mean, na.rm=TRUE)
                #stat_lower = stat_mean - sd(data_fit[data_fit$generation == i & data_fit$r == r,][c(stat)][[1]], na.rm=TRUE)
                #stat_upper = stat_mean + sd(data_fit[data_fit$generation == i & data_fit$r == r,][c(stat)][[1]], na.rm=TRUE)
                data_ribboned = rbind(data_ribboned,data.frame(r=r, gen=i, stat_mean))
            }
        }

        data_ribboned$gen <- (data_ribboned$gen - (100000/scalingfactor)) * scalingfactor
        rownames(data_ribboned) <- c()
                
        data_ribboned$h = rep(h, length(data_ribboned[,1]))
        data_ribboned$sim = rep(sim, length(data_ribboned[,1]))
        df = rbind(df, data_ribboned)
    }
}

df$h <- factor(df$h)
levels(df$h) <- c("paste(italic(h),'=0.0')","paste(italic(h),'=0.5')")


df$sim <- factor(df$sim)
levels(df$sim) <- c("'Model 0'", "'Model 1'", "'Model 2'", "'Model 3'", "'Model 4'")


q <- ggplot(df, aes(colour=factor(r))) + facet_grid(sim~h, scales="free_y", labeller=label_parsed) + geom_line(aes(x=gen,y=stat_mean),linetype=1,alpha=0.75) + geom_vline(xintercept=20000,linetype=2,colour="grey",alpha=0.6) + scale_colour_discrete("recombination rate", labels=c(expression(10^"-9"),expression(10^"-8"),expression(10^"-7"),expression(10^"-6"))) + scale_fill_discrete("recombination rate") + scale_y_continuous(stattitle) + scale_x_continuous(xtitle, limits=c(0,30000)) + theme_bw(base_size=16) + theme(axis.title.y=element_text(angle=90,vjust=0.5))
setwd('~/Dropbox/project_admixture_load/sims_final/')
ggsave("figure_s3.pdf",q,width=10, height=7, units="in")


###Figure S4
setwd('~/Dropbox/project_admixture_load/sims_final/')

library(reshape2)
library(ggplot2)

stat = "p2MeanHomDel"
stattitle = "Mean homozygous derived deleterious sites per individual"
xtitle = "Time after split (generations)"

df <- data.frame(r=character(),gen=character(), stat_mean=character(),sim=character(),h=character())

for (sim in c("model0","model1","model2","model3","model4")){
    for (h in c("0.0","0.5")){
        data = read.csv(sprintf('./h_%s_%s.csv',h,sim))
        data_fit = data[c("rep","r","generation",stat)]
        
        generations = 2000:2600*5
        generations = generations*10/scalingfactor
        rs = c(1e-6, 1e-7, 1e-8, 1e-9)

        data_ribboned = data.frame(r=character(), gen=character(), stat_mean=character())
        data_names = names(data_ribboned)
        
        for (r in rs){
            for (i in generations){
                stat_mean = data_fit[data_fit$generation == i & data_fit$r == r,][c(stat)][[1]]
            	stat_mean[stat_mean==0] <- NA
                stat_mean = mean(stat_mean, na.rm=TRUE)
                #stat_lower = stat_mean - sd(data_fit[data_fit$generation == i & data_fit$r == r,][c(stat)][[1]], na.rm=TRUE)
                #stat_upper = stat_mean + sd(data_fit[data_fit$generation == i & data_fit$r == r,][c(stat)][[1]], na.rm=TRUE)
                data_ribboned = rbind(data_ribboned,data.frame(r=r, gen=i, stat_mean))
            }
        }

        #data_ribboned$r <- data_ribboned$r/10
        data_ribboned$gen <- (data_ribboned$gen - (100000/scalingfactor)) * scalingfactor
        rownames(data_ribboned) <- c()
        
        data_ribboned$h = rep(h, length(data_ribboned[,1]))
        data_ribboned$sim = rep(sim, length(data_ribboned[,1]))
        df = rbind(df, data_ribboned)
    }
}

df$h <- factor(df$h)
levels(df$h) <- c("paste(italic(h),'=0.0')","paste(italic(h),'=0.5')")

df$sim <- factor(df$sim)
levels(df$sim) <- c("'Model 0'", "'Model 1'", "'Model 2'", "'Model 3'", "'Model 4'")

df <- df[df$gen > 19000,]
df$gen = df$gen - (200000/scalingfactor)

q <- ggplot(df, aes(colour=factor(r))) + facet_grid(sim~h, scales="free_y") + geom_line(aes(x=gen,y=stat_mean),linetype=1,alpha=0.75) + geom_vline(xintercept=20000,linetype=2,colour="grey",alpha=0.6) + scale_colour_discrete("recombination rate", labels=c(expression(10^"-9"),expression(10^"-8"),expression(10^"-7"),expression(10^"-6"))) + scale_fill_discrete("recombination rate") + scale_y_continuous(stattitle) + scale_x_continuous(xtitle, limits=c(0,30000)) + theme_bw(base_size=16) + theme(axis.title.y=element_text(angle=90,vjust=0.5))
setwd('~/Dropbox/project_admixture_load/sims_final/')
ggsave("figure_S4.pdf",q,width=10, height=7, units="in")
