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
                   Sex = c(W$sex,rep('F',nrow(NL))))

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

genos <- subset(genos,geno != '00')
genos <- subset(genos,BillLength < 15)

genos$geno <- factor(genos$geno,levels = c('TT','CT','CC'))
genos$PopSex <- ifelse(genos$Pop == 'UK',
                       ifelse(genos$Sex == 'M',
                              'UK (males)',
                              'UK (females)'),
                       'NL (females)')

genos$geno_n <- as.numeric(genos$geno)


#Write to file
write.csv(genos,'NLUK_genotypes_and_bill_length.csv',row.names = F)








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


lrs <- ddply(rs,
             .(BirdID,geno,X33_1,Pop,Sex),
             summarise,
             LRS = sum(N_Recruits,na.rm=T),
             Year = min(as.numeric(paste(Year)),na.rm=T))
lrs$Year <- factor(lrs$Year)


rs$PopSex <- ifelse(rs$Pop == 'UK',
                       ifelse(rs$Sex == 'M',
                              'UK (males)',
                              'UK (females)'),
                       'NL (females)')





