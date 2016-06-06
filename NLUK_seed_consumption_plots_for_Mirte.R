rm(list=ls())

library(plyr)
library(nlme)
library(ggplot2)
library(Rmisc)

dd <- read.csv('NLUK_seed_consumption_genotype.csv')


# Get total seed consumption per season for each bird ---------------------

totals <- ddply(dd,
                .(ring,season,geno,genon,sex),
                summarise,
                totalseeds = sum(seeds),
                month = mean(o_month),
                ndays= length(seeds))
totals <- subset(totals,totalseeds<3000)
totals$seedsperday <- totals$totalseeds/totals$ndays


# Plot seed consumption and number of visits against genotype -------------

FigA <- ggplot(totals,aes(x = season, y = seedsperday,fill = geno))+
  geom_jitter(aes(col = geno),position = position_jitterdodge(jitter.width = 0.2),shape = 1)+
  geom_boxplot(alpha = 0,outlier.colour = NA,coef = 0)+
  theme_bw()

FigB <- ggplot(totals,aes(x = season, y = ndays,fill = geno))+
  geom_jitter(aes(col = geno),position = position_jitterdodge(jitter.width = 0.2),shape = 1)+
  geom_boxplot(alpha = 0,outlier.colour = NA,coef = 0)+
  theme_bw()

FigC <- ggplot(totals,aes(x = season, y = totalseeds,fill = geno))+
  geom_jitter(aes(col = geno),position = position_jitterdodge(jitter.width = 0.2),shape = 1)+
  geom_boxplot(alpha = 0,outlier.colour = NA,coef = 0)+
  theme_bw()

multiplot(FigA,FigB,FigC)

