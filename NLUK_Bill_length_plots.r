rm(list=ls())

library(ggplot2)
library(lme4)
library(glmmADMB)
library(plyr)
library(Rmisc)

dd <- read.csv('NLUK_genotypes_and_bill_length.csv')
dd1 <- read.csv('NLUK_Reproductive_success_and_genotype.csv')

dd1$Year <- factor(dd1$Year)
dd1 <- subset(dd1,(!(is.na(N_Fledglings))))


# Bill length and genotype ------------------------------------------------

dd$CenBill <- dd$BillLength-tapply(dd$BillLength,dd$Ringer,mean)[dd$Ringer] 

boxplot(CenBill~geno,data=subset(dd,Pop == 'UK'))

#Plot Bill length vs genotype for both pops
ggplot(dd,aes(x = as.numeric(geno),y = CenBill,fill = Pop,col = Pop))+
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



temp <- subset(dd1,N_Fledglings>-1)
ggplot(temp,aes(x = (geno), y = N_Fledglings,fill = Pop))+
  geom_boxplot(notch=T)+
  theme_classic()+
  theme(        axis.line.x = element_line(colour = 'black'),
                axis.line.y = element_line(colour = 'black'),
                legend.title = element_blank())+
  ylab('Number of fledglings')+
  xlab('Number of copies of long-billed allele')+
  scale_colour_manual(values = c('gold','darkgrey'))+
  scale_fill_manual(values = c('gold','darkgrey'))




#MODELS - warning, these are fairly slow

#Zero inflated model for Number of Fledglings

m1 <- glmmadmb(N_Fledglings~Sex+as.numeric(geno)*Pop+(1|Year)+(1|BirdID),
               family = 'gaussian', #Also tried running with poisson - still significant
               data=dd1,
               zeroInflation = TRUE)
summary(m1) #Interaction significant

temp <- subset(dd1,!is.na(haplosum))
m1 <- glmmadmb(N_Fledglings~Sex+haplosum+(1|Year)+(1|BirdID),
               family = 'gaussian', #Also tried running with poisson - still significant
               data=temp,
               zeroInflation = TRUE)
summary(m1) #Interaction significant



#Separate binomial and poisson models
dd1$NF2 <- ifelse(dd1$N_Fledglings>0,1,0)
m2 <- glmer(NF2~Sex + as.numeric(geno)*Pop + (1|Year) + (1|BirdID),
               family = 'binomial',
               data=dd1)
summary(m2) #Interaction not significant


m3 <- glmer(N_Fledglings~Sex + as.numeric(geno)*Pop + (1|Year) + (1|BirdID),
               data=subset(dd1,NF2 == 1),
           family = 'poisson')
summary(m3) #Interaction significant - so it's the number of fledglings, not a presence/absence effect


#
# uk <- subset(dd1,Pop == 'UK')
# m4 <- glmmadmb(N_Fledglings~Sex+as.numeric(geno)+(1|Year)+(1|BirdID),
#                family = 'gaussian', #Also tried running with poisson - still significant
#                data=uk,
#                zeroInflation = TRUE)
# summary(m4)
# 
# 
nl <- subset(dd1,Pop == 'NL')
m5 <- glmmadmb(N_Fledglings~as.numeric(geno)+(1|Year)+(1|BirdID),
               family = 'gaussian', #Also tried running with poisson - still significant
               data=nl,
               zeroInflation = TRUE)
summary(m5)


bills <- read.csv('NLUK_bill_length_temporal.csv')
bills <- subset(bills,bl<15.3)
bills <- subset(bills,social_immigrant_status == 'local')

bills$cenbills <- bills$bl-tapply(bills$bl,bills$sex,mean)[bills$sex]

summary(lm(bl~yrbirth,data=bills))



ggplot(bills,aes(x=(yrbirth),y=cenbills))+
  geom_jitter(col = 'grey',width = 0.3)+
  theme_bw() +
  stat_smooth(col = grey(0.2),fill = grey(0.4),se=T,method = 'lm')

bills$bd <- ifelse(is.na(bills$male_bd),bills$female_bd,bills$male_bd)


library(plyr)
temp <- ddply(bills,
              .(yrbirth,sex),
              summarise,
              bill_length = mean(bl))

ggplot(temp,aes(x=yrbirth,y=bill_length,col=sex)) + geom_point() +
  theme_bw() +
  stat_smooth()
  #geom_dotplot(stackdir="center",binaxis="y",binwidth = 0.02,colour="gold",fill="gold") +
  xlab("Year") +
  ylab ("Bill Length") +
  ylim (12.5,15) +
  theme(axis.title.x = element_text(size =24))  +
  theme(axis.title.y = element_text(size =24))  +
  theme(axis.text.x = element_text(size =18))  +
  theme(axis.text.y = element_text(size =18)) +
  scale_x_discrete(breaks=c(1982,1987,1992,1997,2002,2007)) +
  geom_smooth(aes(group=1), method="lm", se=T)

  
  
mbills1 <- read.csv('Tringbills.csv',stringsAsFactors = F)
mbills2 <- read.csv('OUNHMbills.csv',stringsAsFactors = F)

b3 <- read.csv('livebirds.csv')
b3 <- subset(b3,!is.na(Year))

mbills <- rbind(mbills1,mbills2)

temp <- rbind(mbills1,mbills2,b3)

plot(Bill~Year,data = temp)
temp$Live <- ifelse(temp$Year == 2016,'old','new')


temp <- subset(temp,Subspecies == 'newtoni')
anova(lm(Bill~Live,data =temp))

mbills$Month <- sapply(mbills$Date, function(x) strsplit(x,split = '-')[[1]][2])

mbills$NMonth <- NA
mbills$NMonth[mbills$Month == 'Jan'] <- 1
mbills$NMonth[mbills$Month == 'Feb'] <- 2
mbills$NMonth[mbills$Month == 'Mar'] <- 3
mbills$NMonth[mbills$Month == 'Apr'] <- 4
mbills$NMonth[mbills$Month == 'May'] <- 5
mbills$NMonth[mbills$Month == 'Jun'] <- 6
mbills$NMonth[mbills$Month == 'Jul'] <- 7
mbills$NMonth[mbills$Month == 'Aug'] <- 8
mbills$NMonth[mbills$Month == 'Sep'] <- 9
mbills$NMonth[mbills$Month == 'Oct'] <- 10
mbills$NMonth[mbills$Month == 'Nov'] <- 11
mbills$NMonth[mbills$Month == 'Dec'] <- 12


mbills$Season <- ifelse(mbills$NMonth %in% c(3:7),'Summer','Winter')

#mbills <- subset(mbills,Year %in% c(1800:2016))
mbills$UKL <- ifelse(mbills$Subspecies == 'newtoni','UK','Mainland Europe')
mbills$UKLS <- with(mbills,paste(Subspecies,Sex))


mbills <- subset(mbills,!is.na(Bill))
mbills <- subset(mbills,!is.na(Billw))
mbills <- subset(mbills,!is.na(Billd))
pc1 <- princomp(mbills[,c('Bill','Billw','Billd')])

str(pc1)

mbills$pc1 <- pc1$scores[,1]
mbills$pc2 <- pc1$scores[,2]
mbills$pc3 <- pc1$scores[,3]

mbills$cenbills <- mbills$Bill-mean(mbills$Bill)

temp <- subset(mbills,Subspecies!='corsus')
temp <- subset(temp,!is.na(Year))


head(temp)

Fig4A <- ggplot(temp,
       aes(x = UKL,y = cenbills))+
  geom_boxplot(outlier.color = NA,
               notch = T,
               col = 'gold',
               fill = grey(0.4))+
  geom_jitter(col = 'gold')+
  theme_bw()+
  ylab('Mean centred bill length (mm)')+
  xlab('')+
  theme(axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13))+
  ylim(-1.5,1.5)


Fig4B <- ggplot(bills,aes(x=(yrbirth),y=cenbills))+
  geom_jitter(col = 'gold',width = 0.3)+
  theme_bw() +
  stat_smooth(col = grey(0.2),fill = grey(0.4),se=T,method = 'lm')+
  ylab('Mean centred bill length (mm)')+
  xlab('Year of birth')+
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))+
  ylim(-1.5,1.5)


multiplot(Fig4A,Fig4B,layout = matrix(1:2,nrow = 1))








m1 <- lm(Bill~(Year)+Subspecies+Sex+Age,data=temp)
summary(m1)

uk <- subset(mbills,Subspecies == 'newtoni')
summary(lm(Bill~Year+Age+Sex,data = uk))
summary(lm(Bill~Year+Age+Sex+UKL,data = mbills))


ggplot(uk,aes(x = Year,y = Bill,col = Sex))+
  geom_point()+
  stat_smooth(method = 'lm')


summary(lm(bl~yrbirth,data=bills))
summary(lm(Bill~Year,data=mbills))

library(pwr)

myf2 = 0.007474/(1-0.007474)

pwr.f2.test(f2 = myf2,u = 1,power = 0.75)
pwr.f2.test(f2 = myf2,u = 1,v = 149)

