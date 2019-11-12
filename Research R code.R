library(ggplot2)
library(dplyr)
library(reshape2)

# Read mutation data for Head and Neck Carcinoma
HNC.mutationdata <- read.csv("/Users/assefam/OneDrive - Eastern Connecticut State University/CSC 450/CSC 450 Reserch/Mutated Gene data/Head and Neck Squamous Cell Carcinoma (TCGA, PanCancer Atlas).csv", TRUE, sep=",")

# Read mutation data for Invasive breast carcinoma
BIC.mutationdata <- read.csv("/Users/assefam/OneDrive - Eastern Connecticut State University/CSC 450/CSC 450 Reserch/Mutated Gene data/Mutated Genes Breast Invasive Carcinoma TCGA PanCancer Atlas.csv", TRUE, sep=",")

#Read mutation data for Ovarian Serous Cystadenocarcinoma
OVSC.mutationdata <- read.csv("/Users/assefam/OneDrive - Eastern Connecticut State University/CSC 450/CSC 450 Reserch/Mutated Gene data/Muatated Genes Ovarian Serous Cystadenocarcinoma (TCGA, PanCancer Atlas).csv", TRUE, sep=",")


# Extract mutated Genes for cancer mutation data: 

# All mutated genes in Head and Neck Carcinoma
HNC.mutationdata$Gene

  # There are 15900 mutated genes in HNC
  length(HNC.mutationdata$Gene)

# All mutated genes in Invasive breast carcinoma
BIC.mutationdata$ï..Gene

  # There are 16538 mutated genes in IBC
  length(BIC.mutationdata$ï..Gene)


# All mutated genes in ovarian Serous Cystadenocarcinoma
OVSC.mutationdata$Gene

  # There are 12856 mutated genes in OSV
  length(OVSC.mutationdata$Gene)
  
# Extracting Cancer genes from the original dataset using filter the dplyr library.
  cancer.genes.HNC.data <- filter(HNC.mutationdata, Is.Cancer.Gene..source..OncoKB. == "Yes")
  cancer.genes.BIC.data <- filter(BIC.mutationdata, Is.Cancer.Gene..source..OncoKB. == "Yes")
  cancer.genes.OVSC.data <- filter(OVSC.mutationdata, Is.Cancer.Gene..source..OncoKB. == "Yes")

# Extracting cancer genes for each cancer types
  cancer.genes.HNC <- cancer.genes.HNC.data$Gene
  # There are 951 cancer genes for HNC
  length(cancer.genes.HNC)
  
  cancer.genes.BIC <- cancer.genes.BIC.data$ï..Gene
  # There are 946 cancer genes for BIC
  length(cancer.genes.BIC)
  
  cancer.genes.OVSC <- cancer.genes.OVSC.data$Gene
  # There are 821 cancer genes for OVSC
  length(cancer.genes.OVSC)
  
  # Find common genes between HNC and BIC;  This returns index number for common genes 
  common.HNC.and.BIC.index <- match(cancer.genes.HNC, cancer.genes.BIC)
  
  common.HNC.and.BIC <- cancer.genes.BIC[common.HNC.and.BIC.index]
  
  # Getting rid of NA values
  common.HNC.and.BIC.withoutNA<- common.HNC.and.BIC[!is.na(common.HNC.and.BIC)] 
  # There are 889 common cancer genes between BIC and HNC.
  
  
  
  # Find common between HNC, BIC and OVSC
  common.HNC.BIC.OVSC.index <- match( common.HNC.and.BIC.withoutNA, cancer.genes.OVSC)
  
  common.HNC.BIC.OVSC <- cancer.genes.OVSC[common.HNC.BIC.OVSC.index]
  
  # Getting rid of NA values
  common.HNC.BIC.OVSC.withoutNA<- common.HNC.BIC.OVSC[!is.na(common.HNC.BIC.OVSC)] 
  # There are 755 common cancer genes between BIC, HNC, OVSC.
  length(common.HNC.BIC.OVSC.withoutNA)
  
  
  # Exporting commonly muated gene to a .csv file (commonly mutated genes)
  write.csv(common.HNC.BIC.OVSC.withoutNA, file = "common.HNC.BIC.OVSC.withoutNA.csv")
  