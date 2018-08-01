setwd('~/Dropbox/project_admixture_load/sims_final/')

library(reshape2)
library(ggplot2)

scalingfactor = 5

###figure 2
stattitle = expression(paste(frac(italic("w"["R"]),italic("w"["S"]))))
stat = "fitnessRatio"
xtitle = "Time after split (generations)"

df <- data.frame(r=character(), gen=character(), stat_lower=character(), stat_upper=character(), stat_mean=character(), stat_name=character(), sim=character(), h=character())

scalingfactor = 5

for (sim in c("model0","model1","model2","model3","model4")){
    for (h in c("0.0","0.5")){
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
                stat_mean = mean(data_fit[data_fit$generation == i & data_fit$r == r,][c(stat)][[1]], na.rm=TRUE)
                data_ribboned = rbind(data_ribboned,data.frame(r=r, gen=i, stat_lower, stat_median, stat_upper, stat_mean, stat))
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

df.melted <- df

q <- ggplot(df.melted, aes(colour=factor(r))) + facet_grid(sim~h, scales="free_y", labeller=label_parsed) + geom_ribbon(aes(x=gen,ymin=1/stat_lower,ymax=1/stat_upper,fill=factor(r)),alpha=0.2) + geom_line(aes(x=gen,y=1/stat_median),linetype=2) + geom_hline(yintercept=1.0,linetype=2) + geom_vline(xintercept=20000,linetype=2,colour="grey",alpha=0.8) + scale_colour_discrete("recombination rate",labels=c(expression(10^"-9"),expression(10^"-8"),expression(10^"-7"),expression(10^"-6"))) + scale_fill_discrete("recombination rate",labels=c(expression(10^"-9"),expression(10^"-8"),expression(10^"-7"),expression(10^"-6"))) + scale_y_continuous(stattitle) + scale_x_continuous(xtitle) + theme_bw(base_size=14) + theme(axis.title.y=element_text(angle=0,vjust=0.5)) 
ggsave("figure_2.pdf",q,width=10, height=6, units="in")

###figure 3
stat = "p2Fraction"
stattitle = expression(paste(italic("p"["I"])))
xtitle = "Time after split (generations)"

df <- data.frame(r=character(), gen=character(), stat_lower=character(), stat_upper=character(), stat_mean=character(), stat_name=character(), sim=character(), h=character())

for (sim in c("model0","model1","model2","model3","model4")){
    for (h in c("0.0","0.5")){
        data = read.csv(sprintf('./h_%s_%s.csv',h,sim))
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
                stat_mean = mean(data_fit[data_fit$generation == i & data_fit$r == r,][c(stat)][[1]], na.rm=TRUE)
                data_ribboned = rbind(data_ribboned,data.frame(r=r, gen=i, stat_lower, stat_median, stat_upper, stat_mean, stat))
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

q <- ggplot(df, aes(colour=factor(r))) + facet_grid(sim~h, scales="free_y") + geom_ribbon(aes(x=gen,ymin=stat_lower,ymax=stat_upper,fill=factor(r)),alpha=0.2) + geom_line(aes(x=gen,y=stat_mean),linetype=2) + geom_hline(yintercept=0.05,linetype=2) + geom_vline(xintercept=20000,linetype=2,colour="grey",alpha=0.8) + scale_colour_discrete("recombination rate", labels=c(expression(10^"-9"),expression(10^"-8"),expression(10^"-7"),expression(10^"-6"))) + scale_fill_discrete("recombination rate", labels=c(expression(10^"-9"),expression(10^"-8"),expression(10^"-7"),expression(10^"-6"))) + scale_y_continuous(stattitle) + scale_x_continuous(xtitle) + theme_bw(base_size=16) + theme(axis.title.y=element_text(angle=0,vjust=0.5)) 
setwd('~/Dropbox/project_admixture_load/sims_final/')
ggsave("figure_3.pdf",q,width=10, height=7, units="in")


