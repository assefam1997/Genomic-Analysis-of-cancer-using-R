library(ggplot2)
library(dplyr)
library(reshape2)

# Read mutation data for Head and Neck Carcinoma
HNC.mutationdata <- read.csv("/Users/assefam/OneDrive - Eastern Connecticut State University/CSC 450/CSC 450 Reserch/Mutated Gene data/Head and Neck Squamous Cell Carcinoma (TCGA, PanCancer Atlas).csv", TRUE, sep=",")

# Read mutation data for Invasive breast carcinoma
BIC.mutationdata <- read.csv("/Users/assefam/OneDrive - Eastern Connecticut State University/CSC 450/CSC 450 Reserch/Mutated Gene data/Mutated Genes Breast Invasive Carcinoma TCGA PanCancer Atlas.csv", TRUE, sep=",")

#Read mutation data for Ovarian Serous Cystadenocarcinoma
OVSC.mutationdata <- read.csv("/Users/assefam/OneDrive - Eastern Connecticut State University/CSC 450/CSC 450 Reserch/Mutated Gene data/Muatated Genes Ovarian Serous Cystadenocarcinoma (TCGA, PanCancer Atlas).csv", TRUE, sep=",")



# Extracting Cancer genes from the original dataset using filter the dplyr library.
cancer.genes.HNC.data <- filter(HNC.mutationdata, Is.Cancer.Gene..source..OncoKB. == "Yes")
cancer.genes.BIC.data <- filter(BIC.mutationdata, Is.Cancer.Gene..source..OncoKB. == "Yes")
cancer.genes.OVSC.data <- filter(OVSC.mutationdata, Is.Cancer.Gene..source..OncoKB. == "Yes")


# Extract all mutated Genes for cancer mutation data: 

# All mutated genes in Head and Neck Carcinoma
  c <- HNC.mutationdata$Gene

  # There are 15900 mutated genes in HNC
  length(HNC.mutationdata$Gene)

# All mutated genes in Invasive breast carcinoma
  d<- BIC.mutationdata$ï..Gene

  # There are 16538 mutated genes in IBC
  length(BIC.mutationdata$ï..Gene)


  # Finding commonly mutated genes(both cancer and non cancer)
  ll <- match(c,d)
  t <- d[ll]
  t.withoutNa <- t[!is.na(t)]
  lll<- match(t.withoutNa, e)
  tt<-e[lll]
  tt.withoutNa <- tt[!is.na(tt)]
  length(tt.withoutNa)
  
  
# All mutated genes in ovarian Serous Cystadenocarcinoma
  e <- OVSC.mutationdata$Gene

  # There are 12856 mutated genes in OSV
  length(OVSC.mutationdata$Gene)


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
  
  
  # Extracting cancer genes with frequency mutation rate 3% and above. 
  
  # There are 78 cancer genes with frequency mutation of 3% and above in HNC.
  above3percent.HNC <- cancer.genes.HNC.data[1:78, ]
  nrow(above3percent.HNC)
  cancergenes.above3percent.HNC <-  above3percent.HNC$Gene
  
  # There are 18 genes with frequency mutation rate of 3% and above in BIC
  above3percent.BIC <- cancer.genes.BIC.data[1:18, ]
  nrow(above3percent.BIC)
  cancergenes.above3percent.BIC <-  above3percent.BIC$ï..Gene
  
  #There are 13 genes with frequency mutation rate of 3% and above. 
  above3percent.OVSC <- cancer.genes.BIC.data[1:13,]
  nrow(above3percent.OVSC)
  cancergenes.above3percent.OVSC <- above3percent.OVSC$ï..Gene
  
  
  # Find common genes between HNC and BIC with frequency above 3%;  This returns index number for common genes 
  common.HNC.and.BIC.index <- match(cancergenes.above3percent.HNC, cancergenes.above3percent.BIC)
  
  common.HNC.and.BIC <- cancergenes.above3percent.BIC[common.HNC.and.BIC.index]
  
  # Getting rid of NA values
  common.HNC.and.BIC.withoutNA<- common.HNC.and.BIC[!is.na(common.HNC.and.BIC)] 
  # There are 889 common cancer genes between BIC and HNC.
  
  
  
  # Find common between HNC, BIC and OVSC
  common.HNC.BIC.OVSC.index <- match( common.HNC.and.BIC.withoutNA, cancergenes.above3percent.OVSC)
  
  common.HNC.BIC.OVSC <-   cancergenes.above3percent.OVSC[common.HNC.BIC.OVSC.index]
  
  # Getting rid of NA values
  common.HNC.BIC.OVSC.withoutNA<- common.HNC.BIC.OVSC[!is.na(common.HNC.BIC.OVSC)] 
  # There are 755 common cancer genes between BIC, HNC, OVSC.
  length(common.HNC.BIC.OVSC.withoutNA)
  
  
  # Exporting commonly muated gene to a .csv file (commonly mutated genes)
  write.csv(common.HNC.BIC.OVSC.withoutNA, file = "updated common.HNC.BIC.OVSC.withoutNA.csv")
  
  # Generating stacked frequency bargraph for different cancer types 
  
  # generate sample data for three cancer types (named one, two, and three) 
  
  one <- data.frame(genes = c("TP53", "PIK3CA", "PCLO", "KMT2C", "ARID1A", "NCOR1"), Freq = c(0.693,0.175, 0.153, 0.072, 0.035,0.033), cancer = "Head and Neck Squamous carcinoma")
  
  two <- data.frame(genes = c("TP53", "PIK3CA", "PCLO", "KMT2C", "ARID1A", "NCOR1"), Freq = c(0.326 ,0.326,0.036,0.093, 0.038, 0.047), cancer = "Invasive Breast Carcinoma")
  
  three <- data.frame(genes = c("TP53", "PIK3CA", "PCLO", "KMT2C", "ARID1A", "NCOR1"), Freq = c(0.713,0.015,0.027,0.046, 0.008, 0.015), cancer = "Ovarian Serous Cystadenocarcinoma")
  
  
  r <- rbind(one, two, three)
  
  
  
  library(ggplot2)
  
  ggplot(data = r, aes(genes, Freq, fill = cancer)) + geom_col(position = 'fill') + ggtitle("Mutation frequency of each gene in various types of cancer")
  
  
  
  # Generating box plot for three types of cancer and total number of genes mutated
   
  four <- data.frame (number = length(cancer.genes.HNC) , cancer = "Head and Neck Squamous carcinoma")
  five <- data.frame (number = length(cancer.genes.BIC) , cancer = "Invasive Breast Carcinoma")
  six <- data.frame ( number= length(cancer.genes.OVSC) , cancer = "Ovarian Serous Cystadenocarcinoma")
  
  library(gtools)
  r1 <- smartbind(four, five, six, fill = 0)
  
  
  
  
  ggplot(data = r1, aes(cancer, number)) + geom_boxplot() + ggtitle("number of mutated cancer genes for each types of cancer" ) 