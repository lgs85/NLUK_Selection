rm(list=ls())

library(ggplot2)
library(lme4)
library(glmmADMB)

dd <- read.csv('NLUK_genotypes_and_bill_length.csv')
dd1 <- read.csv('NLUK_Reproductive_success_and_genotype.csv')
wsurv <- read.csv('Wytham_survival.csv',colClasses = 'character')

dd1$Year <- factor(dd1$Year)




# Bill length and genotype ------------------------------------------------

#Plot Bill length vs genotype for both pops
ggplot(dd,aes(x = geno,y = BillLength,fill = Pop))+
  geom_boxplot(alpha = 0.2)+
  geom_jitter(aes(col = Pop),position = position_jitterdodge(jitter.width = 0.2))+
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





#Error bar Plot Bill length vs genotype for both pops


temp <- ddply(dd,
              .(geno,Pop),
              summarise,
              Bill = mean(BillLength),
              BillSE = se(BillLength))

ggplot(temp,aes(x = geno,y = Bill,col = Pop))+
  geom_point(position = position_dodge(width = 0.6))+
  geom_errorbar(position = position_dodge(width = 0.6),aes(ymax = Bill+BillSE,ymin = Bill-BillSE))+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'),
                legend.title = element_blank())+
  ylab('Bill length (mm)')+
  xlab('Genotype')+
  scale_colour_manual(values = c('gold','darkgrey'))+
  scale_fill_manual(values = c('gold','darkgrey'))





# Genotype and Number of fledglings ---------------------------------------

#Plot Fledglings vs genotype for both pops
ggplot(subset(dd1,N_Fledglings>0),aes(x = geno,y = N_Fledglings,fill = Pop))+
  geom_boxplot()+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'))+
  ylab('Number of fledglings')+
  xlab('Genotype')+
  scale_colour_manual(values = c('gold','darkgrey','black'))+
  scale_fill_manual(values = c('gold','darkgrey','black'))


dd1 <- subset(dd1,(!(is.na(N_Fledglings))))
m1 <- glmmadmb(N_Fledglings~Sex+(1|Year),
               family = 'gaussian', #Also tried running with poisson - still significant
               data=dd1,
               zeroInflation = TRUE)

dd1$CF <- m1$residuals

ggplot(dd1,aes(x = geno,y = CF,fill = Pop))+
  geom_boxplot(notch=T)+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'))+
  ylab('Number of fledglings')+
  xlab('Genotype')+
  scale_colour_manual(values = c('gold','darkgrey','black'))+
  scale_fill_manual(values = c('gold','darkgrey','black'))





temp <- ddply(dd1,
              .(geno,Pop),
              summarise,
              NF = mean(NF2),
              NFSE = se(NF2))

ggplot(temp,aes(x = geno,y = NF,col = Pop))+
  geom_point(position = position_dodge(width = 0.6))+
  geom_errorbar(position = position_dodge(width = 0.6),aes(ymax = NF+NFSE,ymin = NF-NFSE))+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'),
                legend.title = element_blank())+
  ylab('Number of fledglings')+
  xlab('Genotype')+
  scale_colour_manual(values = c('gold','darkgrey'))+
  scale_fill_manual(values = c('gold','darkgrey'))






#MODELS - warning, these are fairly slow

#Look at the N_Fledglings data

hist(dd1$N_Fledglings) #Looks like a normal distribution with zero inflation.

#Zero inflated model for Number of Fledglings

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
summary(m3) #Interaction significant - so it's the number of fledglings, not a presence/absence effect


