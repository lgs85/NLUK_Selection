rm(list=ls())

library(ggplot2)
library(lme4)
library(glmmADMB)
library(plyr)

dd <- read.csv('NLUK_genotypes_and_bill_length.csv')
dd1 <- read.csv('NLUK_Reproductive_success_and_genotype.csv')

dd1$Year <- factor(dd1$Year)
dd1 <- subset(dd1,(!(is.na(N_Fledglings))))


# Bill length and genotype ------------------------------------------------

#Plot Bill length vs genotype for both pops
ggplot(dd,aes(x = as.numeric(geno),y = BillLength,fill = Pop,col = Pop))+
  stat_smooth(method = 'lm')+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'),
                legend.title = element_blank())+
  ylab('Bill length (mm)')+
  xlab('Number of copies of long-billed allele')+
  scale_colour_manual(values = c('gold','darkgrey'))+
  scale_fill_manual(values = c('gold','darkgrey'))+
  scale_x_continuous(breaks = c(1:3),labels = c(0:2))

  




# Genotype and Number of fledglings ---------------------------------------



temp <- subset(dd1,N_Fledglings>0)
ggplot(temp,aes(x = as.numeric(geno), y = N_Fledglings,col = Pop,fill = Pop))+
  stat_smooth(method = 'glm',method.args = list(family = 'poisson'))+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'),
                legend.title = element_blank())+
  ylab('Number of fledglings')+
  xlab('Number of copies of long-billed allele')+
  scale_colour_manual(values = c('gold','darkgrey'))+
  scale_fill_manual(values = c('gold','darkgrey'))+
  scale_x_continuous(breaks = c(1:3),labels = c(0:2))




#MODELS - warning, these are fairly slow

#Zero inflated model for Number of Fledglings

m1 <- glmmadmb(N_Fledglings~Sex+as.numeric(geno)*Pop+(1|Year)+(1|BirdID),
               family = 'gaussian', #Also tried running with poisson - still significant
               data=dd1,
               zeroInflation = TRUE)
summary(m1) #Interaction significant


#Separate binomial and poisson models
dd1$NF2 <- ifelse(dd1$N_Fledglings>0,1,0)
m2 <- glmer(NF2~Sex + as.numeric(geno)*Pop + (1|Year) + (1|BirdID),
               family = 'binomial',
               data=dd1)
summary(m2) #Interaction not significant


m3 <- lmer(N_Fledglings~Sex + as.numeric(geno)*Pop + (1|Year) + (1|BirdID),
               data=subset(dd1,NF2 == 1),
           family = 'poisson')
summary(m3) #Interaction significant - so it's the number of fledglings, not a presence/absence effect








