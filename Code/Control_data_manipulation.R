
# Load the package
library(biomaRt)
library(edgeR)

#####
## Data manipulation GSE84445

GSE84445 <- read.table("GSE84445_Raw_counts/GSE84445_Raw_counts.txt", header=T)
GSE84445_2<-GSE84445[-12:-21]
GSE84445_final<-GSE84445_2[-22:-31]

ensembl <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl") # create the ensembl object that points to the H sapiens database


convert <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                 filters="hgnc_symbol", 
                 values=GSE84445_final$Transcript_ID,
                 mart = ensembl)

# we add this info to the initial file -> use merge
complete <- merge(GSE84445_final,convert,by.x="Transcript_ID",by.y="hgnc_symbol")

# changing the columns cause we do not need the hugo symbol
complete[1]<-complete[22]
final<-complete[1:21]

colnames(final)[1]<-"ensembl_gene_id"

# write new table for the control
write.table(final,"GSE84445_final.csv", col.names = T, sep=",")

#####
## Data manipulation for GSE227832

GSE227832 <- read.table("GSE227832_RNAseq_read_counts.txt", header=T)
GSE227832_2 <- GSE227832[332:341]
GSE227832_2[11]<-GSE227832[1]
GSE227832_2[1]<-GSE227832_2[11]
GSE227832_2[11]<-GSE227832[332]
colnames(GSE227832_2)[11]<-"X817_B"
colnames(GSE227832_2)[1]<-"ensembl_gene_id"



#####
## mergin the controls -> new data begin in column 12
Controls_merg<-merge(GSE227832_2,final,by.x="ensembl_gene_id",by.y="ensembl_gene_id")

# write new table
write.table(Controls_merg,"Complete_datasets/Controls_merged.csv", col.names = T, sep=",")
