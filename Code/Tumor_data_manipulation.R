library(dplyr)
# library(EDASeq)
library(biomaRt)
library(edgeR)
library(readxl)

##### Make a unique merged dataframe for GSE181157

# Specify the directory path 
getwd()
dire

# List all files that end with ".txt"
file_names <- list.files(directory_path, pattern = "\\.txt$", full.names = TRUE)

# Extract file names without extension and remove ".counts" part
file_names_cleaned <- sub("\\.counts$", "", tools::file_path_sans_ext(basename(file_names)))

# Initialize
final_data2 <- data.frame()

# Loop through each file
for (i in seq_along(file_names)) {
  # Read the file into a temporary data frame
  temp_data <- read.table(file_names[i], col.names = c("ENSEMBL_ID", file_names_cleaned[i]))
  
  # Merge the temporary data frame with the final data frame
  if (nrow(final_data2) == 0) {
    final_data2 <- temp_data
  } else {
    final_data2 <- merge(final_data2, temp_data, by = "ENSEMBL_ID", all = TRUE)
  }
}

# View the final data frame
View(final_data2)

# Write the final data frame to a text file
write.csv(final_data2, file = "merged_GSE181157.csv", row.names = FALSE)



##### Data manipulation GSE227832
GSE227832 <- read.csv('GSE227832_RNAseq_read_counts.txt',header = T,sep = '\t')

# Now we select just the useful data by eliminating the others 
GSE227832 <- GSE227832 %>% dplyr::select(-c(332:341)) %>% 
              dplyr::select(-c(ALL_validation_12,ALL_Validation_20,Constitutional_458rep3,Constitutional_559rep2,ALL_317_2,ALL_317_3,ALL_468_2,ALL_555_2,ALL_680_2))

# substitution of first column with ensembl_gene_id
colnames(GSE227832)[1] <- 'ensembl_gene_id'

#####  Data manipulation GSE133499

GSE133499 <- readxl::read_xlsx('GSE133499_count_matrix_samples_sciarrillo.xlsx')

# remove useless data 
GSE133499 <- GSE133499 %>% dplyr::select(-c(40:43))

# substitution of first column with ensembl_gene_id
colnames(GSE133499)[1] <- 'ensembl_gene_id'

# create the ensembl object that points to the Hsapiens database
ensembl<- biomaRt::useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")

# retrieve the corresponding ensembl_gene_id from the hugo symbol of our dataset
convert_GSE133499 <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol",'tras'),
                 filters="hgnc_symbol", 
                 values=GSE133499$ensembl_gene_id,
                 mart = ensembl)

# we add this info to the initial file -> use merge
GSE133499 <- merge(GSE133499,convert_GSE133499,by.x="ensembl_gene_id",by.y="hgnc_symbol")
# we substitute the hugo symbol with the ensembl genes and eliminate the column with hugo symbol
GSE133499[1]<- GSE133499[40]
GSE133499 <- GSE133499 %>% dplyr::select(-ensembl_gene_id.y)

####### Data manipulation T_all

T_all <- read.csv('T-ALL-RawCount-Cohort7_8.txt',sep= '\t',header = T)

# create the ensembl object that points to the Hsapiens database
ensembl<- biomaRt::useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")

# retrieve the corresponding ensembl_gene_id from the hugo symbol of our dataset
convert_T_all <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                           filters="hgnc_symbol", 
                           values=T_all$GeneSymbol,
                           mart = ensembl)

# we add this info to the initial file -> use merge
T_all <- merge(T_all,convert_T_all,by.x="GeneSymbol",by.y="hgnc_symbol")
T_all[1] <- T_all[110]
T_all <- T_all %>% dplyr::select(-110)
colnames(T_all)[1] <- 'ensembl_gene_id'

###### Data manipulation GSE181157

GSE181157 <- read.csv('GSE181157/f_merged_GSE181157.csv')
colnames(GSE181157)[1]<- 'ensembl_gene_id'

##### Create a final dataframe with all the data 

merg1 <- merge(GSE133499,GSE181157,by = 'ensembl_gene_id',all.x = T,all.y = T)
merg2 <- merge(GSE227832,T_all, by='ensembl_gene_id',all.x = T,all.y = T)
Tumor_Dataframe <- merge(merg1,merg2,by= 'ensembl_gene_id',all.x = T,all.y = T)
Tumor_Dataframe<- na.omit(Tumor_Dataframe)

write.csv(Tumor_Dataframe,file = 'Tumor_dataframe.csv',row.names = F)
 