#
# Quick Sleuth of New Genomic Data
# Plates 92-97
# 2021-11-02 SEJ


library(GenABEL)
library(magrittr)
library(plyr)
library(tidyverse)
library(reshape2)
library(here)

#setwd("../20211109_Plates92-97_for_sleuthing/")

load("20211109_Plates92-97_for_sleuthing/20200323_Plates_1to91.GenABEL.RData")
load("20211109_Plates92-97_for_sleuthing/20211109_Plates_92-97.GenABEL.RData")

idsuccess <- read.table("20211109_Plates92-97_for_sleuthing/20211109_Plates_92-97_ID_Success_Rates.txt", header= T,
                        stringsAsFactors = F, sep = "\t")[,1:5]

plateinfo <- read.table("20211109_Plates92-97_for_sleuthing/SoayPlates92-97ForGenotypingLab.txt", header = T, stringsAsFactors = F, sep = "\t")

pedsnps <- readLines("20211109_Plates92-97_for_sleuthing/ParentageSNPs_431.txt")

PARsnps <- readLines("20211109_Plates92-97_for_sleuthing/PseudoAutosomalSNPs.txt")


clonalsheep <- read.table("20211109_Plates92-97_for_sleuthing/Plates_1to91_Clonal_Sheep.txt", header = T,stringsAsFactors = F)

idsuccess

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Check for cross-plate duplicate samples #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

repeatids <- idnames(soay92_97)[which(idnames(soay92_97) %in% idnames(soay91))]

#~~ Create a genotype table to compare

x91 <- soay91[repeatids,]
x97 <- soay92_97[repeatids,]

x91 <- as.character.gwaa.data(x91) %>% t %>% data.frame
x97 <- as.character.gwaa.data(x97) %>% t %>% data.frame
names(x97) <- paste0(names(x97), "_97")

x91$SNP.Name <- row.names(x91)
x97$SNP.Name <- row.names(x97)

x <- left_join(x91, x97)
x <- x[,sort(names(x))]

for(i in 2:ncol(x)){
  
  x[which(x[,i] == "C/A"),i] <- "A/C"
  x[which(x[,i] == "G/A"),i] <- "A/G"
  x[which(x[,i] == "G/C"),i] <- "C/G"
  x[which(x[,i] == "T/A"),i] <- "A/T"
  
}

for(i in repeatids){
  print(i)
  x1 <- x[,c(paste0("X",i), paste0("X",i,"_97"))]
  x1 <- na.omit(x1)
  print(table(x1[,1] == x1[,2]))
}

# Only a few mismatches are present. Could remove the samples from the newer
# plate.

soay92_97 <- soay92_97[which(!idnames(soay92_97) %in% repeatids),]

rm(x, x1,x91, x97, i, repeatids)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Which IDs are genetically identical?    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ merge the two GenABEL objects together

soaymerge <- merge.gwaa.data(soay91, soay92_97)

#~~ remove poorly typed IDs

x <- perid.summary(soaymerge)
x <- subset(x, CallPP > 0.9)

soaymerge <- soaymerge[row.names(x),]

#~~ Make a rudimentary kinship matrix with pedigree SNPs, flatten & tidy it

pedsnps <-pedsnps[which(pedsnps %in% snpnames(soaymerge))]

gkin <- ibs(soaymerge[,pedsnps])
gkin[upper.tri(gkin)] <- NA
gkin <- data.frame(gkin)
gkin <- cbind(ID1 = row.names(gkin), gkin)
gkin <- melt(gkin)
gkin <- na.omit(gkin)

gkin$variable <- gsub(".", "-", gkin$variable, fixed = T)
gkin$variable <- gsub("X", "", gkin$variable, fixed = T)
gkin$ID1 <- as.character(gkin$ID1)
names(gkin)[2] <- "ID2"
attr(gkin,"na.action") <- NULL

gkin <- subset(gkin, ID1 !=ID2)

str(gkin)

gkin

#~~ Distribution of relatedness

hist(gkin$value, breaks = 100)

#~~ Look at highly related and remove known "clonal" sheep

x <- subset(gkin, value > 0.9)
x <- subset(x, !ID1 %in% clonalsheep$ID1)

x$ID1 %in% idsuccess$Sample.ID
x$ID2 %in% idsuccess$Sample.ID

#~~ remove anything that is 7685 matching itself

temp1 <- grep("7685", x$ID1)
temp2 <- grep("7685", x$ID2)

removelines <- temp1[which(temp1 %in% temp2)]

x <- x[-removelines,]

7658 %in% idsuccess$Sample.ID
7658 %in% idnames(soay91)

#~~ There is an errant 7658 in the data which I will have to look at.

x <- subset(x, ID2 != 7658)

7434 %in% idsuccess$Sample.ID
7434 %in% idnames(soay91)

#~~ The only other Golden Sheep issue is that 12569 is the golden sheep
x

temp1 <- grep("7685", x$ID1)

x <- x[-temp1,]

temp2 <- grep("7685", x$ID2)

x <- x[-temp2,]

x

#12061 11599  - Plate 97 A11, 12061 is matching mum 11599 exactly
#12323 12322  - Plate 93 A10 and B10 are identical
#12861 12167  - Plate 97 E11 and D11 are identical
#12604 12564  - plate 95 D5 and Plate 94 B12 are genetically identical (should be different sexes)
#12834 12684  - Plate 97 F10 and Plate 96 E2 are genetically identical (should be different sexes and different cohorts)

rm(x, temp1, temp2, removelines, x)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Which sexes are wrong?                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

sexsoay <- soaymerge[,which(chromosome(soaymerge) == 27)]
sexsoay <- sexsoay[,which(!snp.names(sexsoay) %in% PARsnps)]

sexF <- perid.summary(sexsoay)
hist(sexF$F)
sexF <- arrange(sexF, F)
sexF$Order <-1:nrow(sexF)

plot(F ~ Order, data =sexF)
