rm(list=ls())
library(data.table)
library(ggplot2)
library(plyr)
library(RMark)
library(lme4)
library(glmmADMB)
library(GenABEL)


wsurv <- read.csv('Wytham_survival.csv',colClasses = 'character')
W <- read.csv('Wytham_morphometrics.csv',stringsAsFactors = F)
NL <- read.csv('NL_morphometrics.csv',stringsAsFactors = F)
geno1 <- fread('AX-100866146NLUK.ped',stringsAsFactors = F)
genoNL <- fread('AX-100866146NLrecent.ped',stringsAsFactors = F)
haps <- read.table('Counts_all_Haps_all_Indis.txt',header = T,stringsAsFactors = F)
wrs <- read.csv('Wytham_reprod_success.csv',stringsAsFactors = F)
nrs <- read.csv('NL_reprod_success.csv',stringsAsFactors = F)
NLid <- read.csv('NL_codes.csv',stringsAsFactors = F)
feed1115 <- read.csv('Lewis 2011-2015 seeds.csv',stringsAsFactors = F)
feed710<- read.csv('Lewis 2007-2010 seeds.csv',stringsAsFactors = F)

NLid <- subset(NLid,!(duplicated(RingNumber)))

#Sort out BirdID
NL$BirdID <- NA
for(i in 1:nrow(NL))
{
  if(NL$indid[i] %in% NLid$IndividualID)
  {
    NL$BirdID[i] <- paste0('NL',subset(NLid,IndividualID == NL$indid[i])$GwasID)
  }
}
NL <- subset(NL,!is.na(BirdID))


genoNLBirdID <- rep(NA,nrow(genoNL))
for(i in 1:nrow(genoNL))
{
  if(genoNL$V2[i] %in% NLid$IndividualID)
  {
    genoNLBirdID[i] <- paste0('NL',subset(NLid,IndividualID == genoNL$V2[i])$GwasID)
  }
}
genoNL$V2 <- genoNLBirdID

#Functions
se <- function(x) sd(x,na.rm=T)/(sqrt(length(x)))


#Combine genotype and bill length data
temp <- data.frame(BirdID = c(W$bto_ring,NL$BirdID),
                   Pop = c(rep('UK',nrow(W)),
                           rep('NL',nrow(NL))),
                   BillLength = c(W$bill_length,NL$BillLength2_gos),
                   Sex = c(W$sex,rep('F',nrow(NL))),
                   Ringer = c(W$initials,rep('NL1',nrow(NL))))

temp <- subset(temp,!(is.na(temp$BillLength)))

geno2 <- rbind(geno1,genoNL)
geno2$geno <- paste0(geno2$V7,geno2$V8)

inboth <- intersect(geno2$V2,temp$BirdID)
temp <- subset(temp,BirdID %in% inboth)
geno3 <- subset(geno2,V2 %in% inboth)
temp <- subset(temp,!duplicated(BirdID))
geno3 <- subset(geno3,!duplicated(V2))
temp <- temp[order(temp$BirdID),]
geno3 <- geno3[order(geno3$V2),]

genos <- temp
genos$geno <- geno3$geno
genos$o_pop <- geno3$V1

genos <- subset(genos,geno != '00')
genos <- subset(genos,BillLength < 15)

genos$geno <- factor(genos$geno,levels = c('TT','CT','CC'))
genos$PopSex <- ifelse(genos$Pop == 'UK',
                       ifelse(genos$Sex == 'M',
                              'UK (males)',
                              'UK (females)'),
                       'NL (females)')

genos$geno_n <- as.numeric(genos$geno)-1


#Write to file
write.csv(genos,'NLUK_genotypes_and_bill_length.csv',row.names = F)
write.table(genos[,c('o_pop','BirdID')],'NLUK_bill_length_birds.txt',row.names = F,col.names = F,quote = F)







#Compare reproductive success vs genotype (both pops)
nrs <- subset(subset(nrs,Include == 1),BroodType == 0)
nrs$BirdID <- paste0('NL',nrs$GWASId)

rs <- data.frame(BirdID = c(wrs$BirdID,nrs$BirdID),
                 Ring = c(wrs$BirdID,nrs$RingNumber),
                 Year = factor(c(wrs$Year,nrs$YearOfBreeding)),
                 N_Fledglings = c(wrs$N_Fledglings,nrs$NoFledged),
                 N_Recruits = c(wrs$N_Recruits,nrs$NoRecruits),
                 Pop = c(rep('UK',nrow(wrs)),rep('NL',nrow(nrs))),
                 Sex = c(wrs$Sex,rep('F',nrow(nrs))))



geno2$BirdID <- geno2$V2
haps$haplosum <- apply(haps[,2:ncol(haps)],1,sum)
head(haps)

rs$geno <- NA
newhap <- mat.or.vec(nrow(rs),ncol(haps)-1)
colnames(newhap) <- colnames(haps)[2:ncol(haps)]
newhap[newhap == 0] <- NA
rs <- cbind(rs,newhap)

for(i in 1:nrow(rs))
{
  if(rs$BirdID[i] %in% geno2$BirdID)
  {
    rs$geno[i] <- subset(geno2,BirdID == rs$BirdID[i])$geno[1]
  }
  if(rs$Ring[i] %in% haps$ID)
  {
    rs[i,9:ncol(rs)] <- subset(haps,ID == rs$Ring[i])[2:ncol(haps)]
  }
}


rs <- subset(rs,geno != '00')
rs$geno <- factor(rs$geno,levels = c('TT','CT','CC'))



rs$PopSex <- ifelse(rs$Pop == 'UK',
                    ifelse(rs$Sex == 'M',
                           'UK (males)',
                           'UK (females)'),
                    'NL (females)')


write.csv(rs, 'NLUK_Reproductive_success_and_genotype.csv',row.names =  F)




lrs <- ddply(rs,
             .(BirdID,geno,X33_1,Pop,Sex,PopSex),
             summarise,
             LRS = sum(N_Recruits,na.rm=T),
             Year = min(as.numeric(paste(Year)),na.rm=T))
lrs$Year <- factor(lrs$Year)

write.csv(lrs, 'NLUK_Lifetime_reproductive_success_and_genotype.csv',row.names =  F)



# Genotype and feeding behaviour ------------------------------------------
feedx <- rbind(feed710[,2:5],feed1115[,2:5])



feedx$month <- as.numeric(sapply(feedx$date,function(x) strsplit(x,split = '-')[[1]][2]))
feedx$yr <- as.numeric(sapply(feedx$date,function(x) strsplit(x,split = '-')[[1]][1]))

feedx$ring <- toupper(feedx$ring)


inboth <- intersect(geno1$V2,feedx$ring)
feed2 <- subset(feedx,ring %in% inboth)
geno1$geno <- with(geno1,paste0(V7,V8))
feed2$season <- NA

for(i in 1:nrow(feed2))
{
  if(feed2$month[i] > 6)
  {
    feed2$season[i] <- paste0(feed2$yr[i],'-',feed2$yr[i]+1)
  } else
  {
    feed2$season[i] <- paste0(feed2$yr[i]-1,'-',feed2$yr[i])
  }
  feed2$geno[i] <- subset(geno1,as.character(paste(V2)) == as.character(paste(feed2$ring[i])))$geno
}

feed2$cenrecs <- log(feed2$recs)-tapply(feed2$recs,feed2$season,function(x) mean(log(x)))[feed2$season] 
feed2$censeeds <- log(feed2$seeds)-tapply(feed2$seeds,feed2$season,function(x) mean(log(x)))[feed2$season] 



feed2$rate <- with(feed2,seeds/recs)

head(feed2) 




table(feed2$season)
feed2 <- subset(feed2,geno !='00')


feed2$genon <- as.numeric(factor(feed2$geno))

summary(lm(censeeds~genon+season,data=feed2))
summary(lm(cenrecs~genon+season,data=feed2))
summary(lm(rate~genon,data=feed2))


m1 <- lmer(cenrecs~month+genon*season+(1|ring)+(1),data=feed2)
m2 <- lmer(cenrecs~month+genon+season+(1|ring),data=feed2)

anova(m1,m2)
confint.merMod(m1)


temp <- ddply(feed2,
              .(ring,season,geno,genon),
              summarise,
              recs = mean(cenrecs),
              seeds = mean(censeeds),
              rate = mean(rate),
              n= length(cenrecs))

hist(temp$n,breaks = 100)

temp <- subset(temp,recs > -1.5)

anova(lm(recs~genon+season,data=temp))



ggplot(temp,aes(x = season, y = recs,fill = geno))+
  geom_jitter(aes(col = geno),position = position_jitterdodge(jitter.width = 0.2))+
  geom_boxplot(notch=F,alpha = 0,col = 'black',outlier.colour = NA)+
  theme_bw()

ggplot(temp,aes(x = season, y = seeds,fill = geno))+
  geom_jitter(aes(col = geno),position = position_jitterdodge(jitter.width = 0.2))+
  geom_boxplot(notch=F,alpha = 0,col = 'black',outlier.colour = NA)+
  theme_bw()

ggplot(temp,aes(x = recs, y = seeds,col = geno))+
  geom_point()+
  stat_smooth(method = 'lm')+
  theme_bw()