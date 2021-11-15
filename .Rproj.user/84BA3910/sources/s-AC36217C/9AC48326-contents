# This script is used for assessing the new SNP chip data
# and for sleuthing.

# 15/11/2021
# Authors: Caelinn / Martin

library(tidyverse)
library(janitor)
library(data.table)
library(WriteXLS)

# Step:
# 1) Filter plink datasets for parentage SNPs, 
# 2) Merge datasets (previous SNP chip and new plate) 
# 3) Calculate relatedness
# 4) Pick relatedness between IDs on current plate and their mothers 
# 5) If they don't match, mixup likely.

# Data:
# - ParentageSNPs_431.txt is file containing one parentage SNP name per line
# - binary PLINK files 20200323_Plates_1to91 (old), 20211109_Plates_92-97 (new)
# - SoayPlates92-97ForGenotypingLab.txt, file containing IDs on new plate and
# potential mother ids from database


# filter main Plink file (containing all previous SNP chip data)
system(paste0("/usr/local/bin/plink --bfile data/20200323_Plates_1to91 --sheep ",     
              "--make-bed --out data/1to91_parentage -extract data/ParentageSNPs_431.txt "))

# filter new SNP chip data
system(paste0("/usr/local/bin/plink --bfile data/20211109_Plates_92-97 --sheep ",    
              "--make-bed --out data/92to97_parentage -extract data/ParentageSNPs_431.txt "))

# merge both into new plink files called "merged.bim/bed/fam"
system(paste0("/usr/local/bin/plink --bfile data/1to91_parentage --sheep ",      
              "--bmerge data/92to97_parentage.bed data/92to97_parentage.bim data/92to97_parentage.fam ",
              "--make-bed --out data/merged "))

# calculate relatedness matrix using merged SNP chip data 
system(paste0("/usr/local/bin/plink --bfile data/merged --sheep ",
              "--make-rel --out data/relmat"))

# Caelinn's code to load lower triangular relatedness matrix (plink output)
# and annotate with IDs | might need to be run on Eddie

relID <- read.table("data/relmat.rel.id")

U <- matrix(0, nrow(relID), nrow(relID))
U[outer(1:nrow(relID), 1:nrow(relID), '<=')] <- scan("data/relmat.rel")

row.names(U) <- relID$V2
colnames(U) <- relID$V2

IDs <- read_delim("data/SoayPlates92-97ForGenotypingLab.txt")
IDs <- data.frame(IDs[,c(3,7)])
All <- unique(c(IDs$ID,IDs$Parent.1..Mum.))

U <- U[which(row.names(U) %in% All),]
U <- U[,which(colnames(U) %in% All)]

# Code to pick Mother-offspring pairs in 
mo <- read_delim("data/SoayPlates92-97ForGenotypingLab.txt") %>% 
        clean_names()

add_rel <- function(id, mum_id) {
        if (!is.na(mum_id)) {
                out <- U[rownames(U) == as.character(id), colnames(U) == as.character(mum_id)][1]
                out  
        } else {
                out <- NA
        }
}
# find relatedness for mother offspring pairs in matrix
rel <- map2(mo$id, mo$parent_1_mum, add_rel) %>% 
        unlist()

# add to mo data.frame
mo$rel <- rel

# check whether it makes sense
hist(mo$rel)

# write to table
WriteXLS(mo, "mo_rel.xls")
              