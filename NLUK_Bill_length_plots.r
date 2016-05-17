rm(list=ls())

library(ggplot2)
library(lme4)
library(glmmADMB)

dd <- read.csv('NLUK_genotypes_and_bill_length.csv')
dd1 <- read.csv('NLUK_Reproductive_success_and_genotype.csv')

dd1$Year <- factor(dd1$Year)




# Bill length and genotype ------------------------------------------------

#Plot Bill length vs genotype for both pops
ggplot(dd,aes(x = geno,y = BillLength,fill = Pop))+
  geom_boxplot()+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'),
                legend.title = element_blank())+
  ylab('Bill length (mm)')+
  xlab('Genotype')+
  scale_colour_manual(values = c('gold','darkgrey'))+
  scale_fill_manual(values = c('gold','darkgrey'))

  
  
#lm of bill-length vs genotype
summary(lm(BillLength~geno_n+Pop+Sex,data=dd))




# Genotype and Number of fledglings ---------------------------------------

#Plot Fledglings vs genotype for both pops
ggplot(dd1,aes(x = geno,y = N_Fledglings,fill = Pop))+
  geom_boxplot(notch=T)+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'))+
  ylab('Number of fledglings')+
  xlab('Genotype')+
  scale_colour_manual(values = c('gold','darkgrey','black'))+
  scale_fill_manual(values = c('gold','darkgrey','black'))


#MODELS - warning, these are fairly slow

#Look at the N_Fledglings data

hist(dd1$N_Fledglings) #Looks like a normal distribution with zero inflation.

#Zero inflated model for Number of Fledglings
dd1 <- subset(dd1,(!(is.na(N_Fledglings))))
m1 <- glmmadmb(N_Fledglings~Sex+as.numeric(geno)*Pop+(1|Year)+(1|BirdID),
               family = 'gaussian', #Also tried running with poisson - still significant
               data=dd1,
               zeroInflation = TRUE)
summary(m1) #Interaction significant


#Separate binomial and poisson models
dd1$NF2 <- ifelse(dd1$N_Fledglings>0,1,0)
m2 <- glmmadmb(NF2~Sex + as.numeric(geno)*Pop + (1|Year) + (1|BirdID),
               family = 'binomial',
               data=dd1,
               zeroInflation = FALSE)
summary(m2) #Interaction not significant

m3 <- glmmadmb(N_Fledglings~Sex + as.numeric(geno)*Pop + (1|Year) + (1|BirdID),
               family = 'gaussian',
               data=subset(dd1,NF2 == 1),
               zeroInflation = FALSE)
summary(m3) #Interaction significant - so it's the number of fledglings, not a presence/abssence effect



