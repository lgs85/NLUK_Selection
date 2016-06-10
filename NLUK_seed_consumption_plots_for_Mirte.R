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






#NEED TO LOOK AT REPRODUCTIVE SUCCESS, FEEDING AND GENOTYPE
reprod$ryear <- with(reprod,paste(Ring,Year))
totals$fyear <- with(totals,paste(ring,substring(season,6)))

inboth <- intersect(reprod$ryear,totals$fyear)

reprod2 <- subset(subset(reprod,ryear %in% inboth),!duplicated(ryear))
totals2 <- subset(totals,fyear %in% inboth)

reprod2 <- reprod2[order(reprod2$ryear),]
totals2 <- totals2[order(totals2$fyear),]

totals2$N_Fledglings <- reprod2$N_Fledglings


ggplot(totals2,aes(y = N_Fledglings, x = totalseeds,col = geno))+
         geom_point()+
         stat_smooth(method = 'lm')



