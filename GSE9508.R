#########################################################
## Victoria Bartolotta                                   #
# CSC 450 Senior Research Project                       #
# Analyzing GSE9508                                     #
# comparing metastatic and non-metastatic               #
# osteosarcoma samples to find differentially           #
# expressed genes using p value 0.001                   #
#########################################################

# load libraries 
library(GEOquery)
library(limma)
library(ggplot2)

# load data
GSE9508 <- getGEO("GSE9508")
GSE9508.p <- pData(GSE9508[[1]])
GSE9508.expr <- exprs(GSE9508[[1]])

# # of probes, # of samples
dim(GSE9508.expr)

# only look at metastatic vs. nonmetastic only, not control group
b <- grep("metastatic", GSE9508.p$title)
GSE9508.expr <- GSE9508.expr[,b]
GSE9508.p <- GSE9508.p[b,]

# # of probes, # of samples
dim(GSE9508.expr)

# show behaivor of data
hist(GSE9508.expr, col = "gray", main="GSE9508 - Histogram")


# Boxplot to show normalized data
boxplot(GSE9508.expr, col=c("darkgreen", "darkgreen", "darkgreen",
                            "darkred", "darkred", "darkred"),
        ylab = "normalized value", 
        main="GSE9508 - boxplots", las=2)

# create column "new" to set each patient status to 
# NM = nonmetastic, MET = metastatic
# the last 13 samples are nonmetastic, so this works
GSE9508.p$new <- c("MET")
GSE9508.p$new[22:34] <- c("NM")

# how many MET vs. NM patients are there? 21 / 13
data <- as.character(GSE9508.p$new)
table(data)

# construct design matrix
design <- model.matrix(~0+data)
colnames(design) <- c("metastic", "non_metastatic")
head(design)

# fit linear model
fit <- lmFit(GSE9508.expr, design)

# specify contrasts 
contrast.matrix <- makeContrasts(metastic - non_metastatic, levels=design)
contrast.matrix

## fit model based on contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
head(fit2$coefficients)
head(fit2$sigma)
fit2 <- eBayes(fit2)

## get top probes, sorted by p-value (gives top 10 genes by default)
tt <- topTable(fit2,sort.by = "p")
head(tt)

# total number of probes in dataset
# and p value
tt.05 <- topTable(fit2,sort.by = "p", p.value = 0.05, number = nrow(GSE9508.expr))
head(tt.05)
nrow(tt.05)

# create and output data frame contraining the top 3 probes
# w/ each logFC and FDR only
df <- data.frame(tt.05$logFC, tt.05$adj.P.Val)
df %>% top_n(3)

# boxplot showing probe with lowest adj.p value 
probe <- rownames(tt.05)[1]
probe
logFC <- tt.05[1,]$logFC
2**logFC
FC <- paste0("FC = ", round(2**logFC, 2))
main <- paste0("Expression of ", probe, ", ", FC)
m <- match(probe, rownames(GSE9508.expr))
df1 <- data.frame(expr = GSE9508.expr[m,], FC = FC)

ggplot(df1, aes(x = FC, y = expr, fill = FC)) + geom_boxplot() +
  ylab("FC and FDR for probe") + ggtitle(main) +
  scale_fill_manual(values = c("pink", "blue")) +
  theme_classic() + theme(legend.position = "none")

# GPL platform of the data being analyzed
platform <- annotation(GSE9508[[1]])   
platform

# using getGEO, download the platform (GPL) for this data
pl <- getGEO(platform)
pl <- Table(pl)

# find the gene names corresponding to all probes and create a table
# w/ corresponding gene names, probe names, logFC, and adjusted p-values 
# only 
probe <- rownames(tt.05)[1]
m <- match(probe, pl$ID)
pl$`GENE_SYMBOL`[m]
# wdr1

g <- grep("WDR1", pl$`GENE_SYMBOL`)
pl <- pl[g,]
g <- grep("\\bWDR1\\b", pl$`GENE_SYMBOL`)
pl$ID[g]
# "8038" "29344" "36049" "38571"

genes <- tt.05[g,]
# t1 <- data.frame(genes$logFC, genes$adj.P.Val, pl$`GENE_NAME`, pl$ID)
# View(t1)


# using DAVID, identify Gene Otology terms and KEGG pathways
# associated w/ the differentially expressed probes 
pl <- getGEO(platform)
pl <- Table(pl)
probes <- rownames(tt.05)
m <- match(probes, pl$ID)
genes <- pl$`GENE_SYMBOL`[m]
keep <- genes!=""
genes <- genes[keep]
genes <- strsplit(genes, " /// ")
genes <- unlist(genes)
genes <- unique(genes) 

# list top 10 differentially expressed genes
genes[1:10]

# export list of genes to .csv to compare lists of genes between data sets
write.csv(genes,'C:\\Users\\torib\\Documents\\school\\fall 19\\csc450\\final\\gse9508.csv', row.names = FALSE)
