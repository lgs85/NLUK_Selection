rm(list=ls())

library(plyr)
library(nlme)
library(lme4)
library(ggplot2)
library(Rmisc)
library(arm)

dd <- read.csv('NLUK_seed_consumption_genotype.csv')
reprod <- read.csv('NLUK_Reproductive_success_and_genotype.csv')

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

ggplot(totals,aes(x = as.numeric(geno), y = log10(totalseeds)))+
  stat_smooth(method = 'lm',col = 'darkgrey')+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'),
                legend.title = element_blank())+
  ylab('log10(seeds consumed)')+
  xlab('Number of copies of long-billed allele')+
  scale_x_continuous(breaks = c(1:3),labels = c(0:2))



#Models
m1 <- lme(log(totalseeds)~genon+sex+month,random = ~1|season,data=totals)
summary(m1)

m1 <- lme((seedsperday)~genon+sex+month,random = ~1|season,data=totals)
summary(m1)

m1 <- lme(log(ndays)~genon+sex+month,random = ~1|season,data=totals)
summary(m1)


