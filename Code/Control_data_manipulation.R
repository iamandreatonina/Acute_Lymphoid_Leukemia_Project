
# Load the package
library(biomaRt)
library(edgeR)
library(tidyverse)

#### Data manipulation GSE84445, 20 healthy samples ####
setwd("~/Desktop/clone_ Acute_Lymphoid_Leukemia_Project /data/Datasets/Controls")

GSE84445 <- read.table("GSE84445_Raw_counts.txt", header=T)

# Selecting only the samples healthy we need, the other are discarded 
GSE84445<-GSE84445[-12:-21]
GSE84445<-GSE84445[-22:-31]

# create the ensembl object that points to the H.sapiens database
ensembl <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")

#retrive the infromation we want 
convert <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                 filters="hgnc_symbol", 
                 values=GSE84445$Transcript_ID,
                 mart = ensembl)

# we add this info to the initial file -> use merge
GSE84445 <- merge(GSE84445,convert,by.x="Transcript_ID",by.y="hgnc_symbol")

# changing the columns cause we do not need the Hugo-symbol
GSE84445[1]<-GSE84445[22]
GSE84445<-GSE84445[1:21]

colnames(GSE84445)[1]<-"ensembl_gene_id"

##### Data manipulation for GSE227832, 10 healthy ####

# We select only the healty samples we need
GSE227832 <- read.table("../Tumors/GSE227832_RNAseq_read_counts.txt", header=T)
GSE227832_2 <- GSE227832[332:341]
GSE227832_2[11]<-GSE227832[1]
GSE227832_2[1]<-GSE227832_2[11]
GSE227832_2[11]<-GSE227832[332]
colnames(GSE227832_2)[11]<-"X817_B"
colnames(GSE227832_2)[1]<-"ensembl_gene_id"

#### Data manipulation GSE139073,all 40 healthy controls ####

GSE139073<- read.table('GSE139073_raw_counts.txt',header = T)
GSE139073$ensembl_gene_id <- row.names(GSE139073)
rownames(GSE139073)<- NULL
GSE139073 <- relocate(GSE139073,'ensembl_gene_id',.before='S1207')


#### Data manipulation GSE190269, all 5 samples are healthy controls, mesenchimal cells, no gene in common with the other cells.####

# GSE190269<- read.table('GSE190269.txt/GSE190269_SJCD34_n5_hg19_counts.txt', header = T)
# GSE190269_2 <-GSE190269[7:12] 
# GSE190269_2[1] <- GSE190269[1]
# colnames(GSE190269_2)[1] <- 'ensembl_gene_id'
# -> not take into consideration because lead everything to zero !!

#### Data manipulation GSE115736, extraction of 18 samples which are healthy controls ####
GSE115736<- read.table('GSE115736_Haemopedia-Human-RNASeq_raw.txt/GSE115736_Haemopedia-Human-RNASeq_raw.txt',header = T)
GSE115736 <- select(GSE115736,CD4T.1,CD4T.2,CD4T.3,CD4T.4,CD4T.5,CD8T.1,CD8T.2,CD8T.3,CD8T.4,CD8T.5,MemB.1,MemB.2,MemB.3,NveB.1,NveB.2,NveB.3,NveB.4,NveB.5)
GSE115736$ensembl_gene_id <- row.names(GSE115736)
rownames(GSE115736)<- NULL
GSE115736 <- relocate(GSE115736,'ensembl_gene_id',.before='CD4T.1')

#### Mergin the controls -> new data begin in column 12 ####
momentary_list <- list(GSE84445,GSE227832_2,GSE139073,GSE115736) 
final <- momentary_list %>%  reduce(full_join,by = 'ensembl_gene_id') %>% drop_na()


# write new table
write.table(final,"../Post_manipulation/Controls_merged.csv", col.names = T, sep=",")
