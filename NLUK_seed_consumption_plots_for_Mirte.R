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
                .(ring,season,geno,genon,sex,haplosum),
                summarise,
                totalrecs = sum(recs),
                month = mean(o_month),
                ndays= length(recs))
#totals <- subset(totals,totalrecs<4000)
totals$recsperday <- totals$totalrecs/totals$ndays


# Plot seed consumption and number of visits against genotype -------------

ggplot(totals,aes(x = as.numeric(geno), y = log10(totalrecs)))+
  stat_smooth(method = 'lm',col = 'darkgrey')+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'),
                legend.title = element_blank())+
  ylab('log10(number feeder visits)')+
  xlab('Number of copies of long-billed allele')+
  scale_x_continuous(breaks = c(1:3),labels = c(0:2))

ggplot(totals,aes(x = haplosum, y = (totalrecs)))+
  geom_point()+
  stat_smooth(col = 'darkgrey')


reprod <- subset(reprod,!is.na(haplosum))


myhaps <- reprod[,c(9:25)]

temp <- data.frame(nhaps = unlist(myhaps),
                   hapID = rep(colnames(myhaps),each = nrow(myhaps)),
                   N_Fledglings = rep(reprod$N_Fledglings,ncol(myhaps))
                   )
ggplot(temp,aes(x = nhaps, y = (N_Fledglings),col = hapID))+
  stat_smooth(method = 'lm',se = F)



uhaps <- as.character(unique(temp$hapID))
estimates <- rep(NA,length(uhaps))
pvals <- rep(NA,length(uhaps))

for(i in 1:length(uhaps))
{
  temp2 <- subset(subset(temp,hapID == uhaps[i]))
  m1 <- summary(glm(N_Fledglings~nhaps,data = temp2,family = 'poisson'))
  pvals[i] <- m1$coefficients[2,4]
  estimates[i] <- m1$coefficients[2,1]
}

t1 <- cbind(uhaps,estimates,pvals)

df = 2*length(pvals)
pchisq( -2*sum(log(pvals)), df, lower.tail=FALSE)


sigs <- data.frame(subset(t1,as.numeric(pvals)<0.05))
sigs$L95 <- NA
sigs$U95 <- NA
sigs$LP <- NA
sigs$UP <- NA

nreps = 1000

for(i in 1:nrow(sigs))
{
  t_out <- rep(NA,nreps)
  p_out <- rep(NA,nreps)
  cd <- subset(temp,hapID == paste(sigs$uhaps[i]))
  for(j in 1:nreps)
  {
    cd$N_Fledglings <- cd$N_Fledglings[sample(c(1:length(cd$N_Fledglings)),length(cd$N_Fledglings),replace = F)]
    m1 <- summary(glm(N_Fledglings~nhaps,data = cd,family = 'poisson'))
    t_out[j] <- m1$coefficients[2,1]
    p_out[j] <- m1$coefficients[2,4]
  }
  sigs$U95[i] <- quantile(t_out,0.95)
  sigs$L95[i] <- quantile(t_out,0.05)
  sigs$UP[i] <- quantile(p_out,0.95)
  sigs$LP[i] <- quantile(p_out,0.05)
}

sigs




ggplot(temp2,aes(x = factor(nhaps), y = (N_Fledglings)))+
  geom_boxplot(notch=T)

totals$month2 <- poly(totals$month,2)[,2]
#Models
m1 <- lme(log(recs)~genon+sex+month+season,random = ~1|ring,data=dd)
summary(m1)

m1 <- lme(log(totalrecs)~genon+sex,random = ~month2|season,data=totals)
summary(m1)

m1 <- lme(log(totalrecs)~haplosum+sex,random = ~month2|season,data=totals)
summary(m1)

m1 <- lme((recsperday)~genon+sex+month,random = ~1|season,data=totals)
summary(m1)

m1 <- lme(log(ndays)~genon+sex,random = ~month2|season,data=totals)
summary(m1)

hist(log(totals$ndays ))
