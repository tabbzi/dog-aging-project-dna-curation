#!/usr/bin/env Rscript

# Libraries:
library(caret)
library(randomForest)
library(scales)
library(tidyverse)
library(data.table)
library(stringr)

# Model stored in .RData must have the following:
# `model` - a randomForest that can be called by predict()
# `modelMap` - a data frame containing: CHR, POS, modSNP, modID, COUNTED, ALT
# modID is the ID used in the randomForest definition

# Arguments:
args <- commandArgs(TRUE)
rdata <- as.character(args[2])
load(rdata)

args <- commandArgs(TRUE)
predfile <- as.character(args[1])
topsnpfile <- as.character(args[3])
mod <- as.character(args[4])
print(args)

# read in .traw file
dataQuery <- read.table(predfile, stringsAsFactors=F, header=T)
dataQuery = dataQuery %>% filter(ALT != "*")

# read in top snps
topSNPs <- read.table(topsnpfile, sep = ",", stringsAsFactors=F, header = T)

# get the list of IIDs; assume number of dogs >= 1; in traw, FID_IID starts at 7th col
# extract FID_IID
dogs <- colnames(dataQuery)[-c(1:6)]
dogs <- gsub("X","",dogs)

# extract IID
for(i in 1:length(dogs)){
  dogs[i] <- strsplit(dogs[i],"_")[[1]][2]
  }

# CALLED & MATCH ALLELES: join to query .traw by chromosome, position, counted allele, alt allele
dataMatch = modelMap %>%
  merge(dataQuery,
        by = c("CHR","POS","COUNTED","ALT"),
        all.x = T)

missingsnp = dataMatch %>%
  filter(is.na(SNP) & !is.na(modSNP)) %>%
  pull(modSNP)

dataMatch = dataMatch %>%
  filter(!is.na(SNP) & !is.na(modSNP))

# CALLED & SWAPPED ALLELES: for missing SNPs that did not match by counted and alt allele, check if reverse matches, then reverse genotype calls
dataSwap = modelMap %>%
  filter(modSNP %in% missingsnp) %>%
  merge(dataQuery,
        by.x = c("CHR","POS","ALT","COUNTED"),
        by.y = c("CHR","POS","COUNTED","ALT")) %>%
  mutate_at(vars(-c("CHR","POS","ALT","COUNTED","modSNP","modID","SNP","X.C.M")),
            ~case_when(as.character(.) == "0" ~ 2,
                       as.character(.) == "1" ~ 1,
                       as.character(.) == "2" ~ 0))

# UNCALLED: SNPs that remain in model and are uncalled in query dogs
missingsnp = setdiff(missingsnp,dataSwap$modSNP)
dataMiss = modelMap %>%
  merge(dataQuery,
        by = c("CHR","POS","COUNTED","ALT"),
        all.x = T) %>%
  filter(modSNP %in% missingsnp)

# combine results:
data = rbindlist(list(dataMiss,dataMatch,dataSwap), use.names = T) %>% arrange(CHR,POS)

# extract the allele counts information from new sample
# change it to the format that can be taken in for predict()
# for the variants not called in this sample, replace with zero
d = data %>%
  arrange(modID) %>%
  select(-c("CHR","POS","modSNP","modID","SNP","X.C.M","COUNTED","ALT")) %>%
  replace(is.na(.), 0)
rownamesSNP = data$modID
rownamesSNP = data %>%
  arrange(modID) %>%
  pull(modID)
d <- t(d)
colnames(d) = rownamesSNP

# size prediction
y <- predict(model, d, type = "response")

# generate phenotype table
phenoTable <- data.frame(id = dogs, prediction = y)
write.table(phenoTable, file = paste0(predfile, ".", mod, ".phenotypes.csv"), sep = ",", quote = F, row.names = F)

# generate genotype table
genoTable = data %>%
  replace(is.na(.), 0) %>%
  merge(topSNPs, by = c("CHR","POS")) %>%
  mutate(model_effect = if_else(COUNTED==alt,
                 alt_allele_effect,
                 ref_allele_effect)) %>%
  pivot_longer(cols = -c(CHR,
                         POS,
                         COUNTED,
                         ALT,
                         modSNP,
                         modID,
                         SNP,
                         X.C.M,
                         alt,
                         ref,
                         alt_allele_effect,
                         ref_allele_effect,
                         model_effect,
                         name,
                         gene,
                         possibleAlleles),
               names_to = "id") %>%
  mutate(effect = if_else(value==1,
                          0,
                          if_else(value==0,
                                  if_else(model_effect=="increase",
                                          -1,
                                          1),
                                  if_else(model_effect=="increase",
                                          1,
                                          -1))),
         firstCopy = if_else(value==0,
                             ALT,
                             if_else(value==1,
                                     ALT,
                                     COUNTED)),
         secondCopy = if_else(value==0,
                              ALT,
                              if_else(value==1,
                                      COUNTED,
                                      COUNTED))) %>%
  mutate(id = gsub("X0_","",id)) %>%
  select(id,name,gene,possibleAlleles,firstCopy,secondCopy,effect)
write.table(genoTable, file = paste0(predfile, ".", mod, ".genotypes.csv"), sep = ",", row.names = F)
