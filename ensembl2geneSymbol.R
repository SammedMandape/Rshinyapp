## biomart conversion of prostate_fpkm file from ensembl to gene symbol
## and keeping only unique gene symbols
## Author: @Sammed Mandape, April 2019
## smandape@email.arizona.edu

library(biomaRt)
library(dplyr)
tab2rows <- read.table("Prostate_FPKM", header = TRUE, nrows = 2)
#browser()
classes <- sapply(tab3rows, class)
fpkm <- read.table("Prostate_FPKM", header = TRUE, row.names=NULL, colClasses = classes,comment.char = "")
colnames(fpkm)[colnames(fpkm) == 'row.names'] <- 'ENSEMBL'
#fpkm[1:2,1:2]

#biomart ensembl gene id to gene symbol 
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl2gene<-getBM(filters = 'ensembl_gene_id', attributes = c('ensembl_gene_id','hgnc_symbol'),values = fpkm$ENSEMBL, mart=mart)

#write.table(ensembl2gene,file="Prostate_FPKM_Ensembl2GeneSymbol" , sep = "\t", row.names = FALSE)

# to see which rows are empty (ensemble id not matched to any gene symbol)
head(ensembl2gene[which(ensembl2gene$hgnc_symbol == ""),])

# remove rows where hgnc_genesymbol is empty
ensembl2gene_nonempty <- ensembl2gene[!(ensembl2gene$hgnc_symbol == ""),]
nrow(ensembl2gene_nonempty) #36854
head(ensembl2gene_nonempty)

#write.table(ensembl2gene_nonempty,file="Prostate_FPKM_Ensembl2GeneSymbol_nonempty" , sep = "\t", row.names = FALSE)
#fpkm[1:2,1:2]

#merge biomart results with prostate_fpkm
fpkm_geneSymbol <- merge(ensembl2gene_nonempty, fpkm, by.x = 'ensembl_gene_id', by.y = 'ENSEMBL')
nrow(fpkm_geneSymbol) #36854 

# keep only unique gene symbols
fpkm_geneSymbol_unique  <- (fpkm_geneSymbol[!duplicated(fpkm_geneSymbol[,2]),])
#head(fpkm_geneSymbol[,-1])

# remove first column with ensembl gene ids and keep only gene symbols
fpkm_geneSymbol_unique <- fpkm_geneSymbol_unique[,-1]
#fpkm_geneSymbol_unique[1:2,1:3]
write.table(fpkm_geneSymbol_unique,file="Prostate_FPKM_Ensembl2GeneSymbol_unique" , sep = "\t", row.names = FALSE)
