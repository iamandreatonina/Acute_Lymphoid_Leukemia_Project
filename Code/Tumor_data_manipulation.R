library(dplyr)
# library(EDASeq)
library(biomaRt)
library(edgeR)
library(readxl)

##### Make a unique merged dataframe for GSE181157

# Specify the directory path, to change! 
directory_path <- 'path of the direcotry where you containe the GSE181157_RAW folder '

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
write.csv(final_data2, file = "f_merged_GSE181157.csv", row.names = FALSE)

###### Data manipulation GSE181157, 173 pediatric tumor samples 

# setwd("~/Desktop/clone_ Acute_Lymphoid_Leukemia_Project /data/Datasets/Tumors")

GSE181157 <- read.csv('f_merged_GSE181157.csv')
colnames(GSE181157)[1]<- 'ensembl_gene_id'


##### Data manipulation GSE227832, 321 tumor samples, both pediatric and adults, without the control samples and the 'double' samples for some patients 
GSE227832 <- read.csv('GSE227832_RNAseq_read_counts.txt',header = T,sep = '\t')

# Now we select just the useful data by eliminating the others which are controls 
GSE227832 <- GSE227832 %>% dplyr::select(-c(332:341)) %>% 
              dplyr::select(-c(ALL_validation_12,ALL_Validation_20,Constitutional_458rep3,Constitutional_559rep2,ALL_317_2,ALL_317_3,ALL_468_2,ALL_555_2,ALL_680_2))

# substitution of first column with ensembl_gene_id
colnames(GSE227832)[1] <- 'ensembl_gene_id'

#####  Data manipulation GSE133499, 38 tumor samples pediatric
GSE133499 <- readxl::read_xlsx('GSE133499_count_matrix_samples_sciarrillo.xlsx')

# remove useless data, which are other samples unknown 
GSE133499 <- GSE133499 %>% dplyr::select(-c(40:43))

# substitution of first column with ensembl_gene_id
colnames(GSE133499)[1] <- 'ensembl_gene_id'

# create the ensembl object that points to the Hsapiens database
ensembl<- biomaRt::useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")

# retrieve the corresponding ensembl_gene_id from the hugo symbol of our dataset
convert_GSE133499 <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                 filters="hgnc_symbol", 
                 values=GSE133499$ensembl_gene_id,
                 mart = ensembl)

# we add this info to the initial file -> use merge
GSE133499 <- merge(GSE133499,convert_GSE133499,by.x="ensembl_gene_id",by.y="hgnc_symbol")
# we substitute the hugo symbol with the ensembl genes and eliminate the column with hugo symbol
GSE133499[1]<- GSE133499[40]
GSE133499 <- GSE133499 %>% dplyr::select(-ensembl_gene_id.y)

####### Data manipulation T_all, 107 tumor samples, both pediatric and adults. 

T_all <- read.csv('T-ALL-RawCount-Cohort7_8.txt',sep= '\t',header = T)

# create the ensembl object that points to the Hsapiens database
ensembl<- biomaRt::useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")

# retrieve the corresponding ensembl_gene_id from the hugo symbol of our dataset
convert_T_all <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                           filters="hgnc_symbol", 
                           values=T_all$GeneSymbol,
                           mart = ensembl)

# we add this info to the initial file -> use merge and get ensembl gene id 
T_all <- merge(T_all,convert_T_all,by.x="GeneSymbol",by.y="hgnc_symbol")
T_all[1] <- T_all[110]
T_all <- T_all %>% dplyr::select(-110)
colnames(T_all)[1] <- 'ensembl_gene_id'


###### Data manipulation GSE26530, tumor samples are 28 age unknown. -> 
GSE26530 <- read.table('GSE26530/raw_counts_GSE26530.csv', sep = ',',header = T) %>% t()

colnames(GSE26530) <- GSE26530[1,] # setted the colnames as samples 
GSE26530 <- GSE26530[-1,] # eliminate the first row, which containe the samples name 
GSE26530 <- as.data.frame(GSE26530) #retrasformina into a dataframe 
GSE26530$ensembl_gene_id <- rownames(GSE26530) # create column ensembl_gene_id
rownames(GSE26530) <- NULL
GSE26530 <- relocate(GSE26530,'ensembl_gene_id',.before='SRR627491')

# Need to remove control samples inside the dataframe, which are SRR627491,SRR627492,SRR627493,SRR627494,SRR627495,SRR627496,SRR627497,SRR627498
`%nin%` <- Negate(`%in%`)
controls <- c('SRR627491','SRR627492','SRR627493','SRR627494','SRR627495','SRR627496','SRR627497','SRR627498')

GSE26530 <- GSE26530[which(colnames(GSE26530) %nin% controls)]


##### Data manipulation GSE228632, 65 samples all tumor samples, bone marrow, all pediatric 

GSE228632 <- read.table('GSE228632_RNAseq_read_counts.txt', sep = '\t',header = T)

#change colname of the genes
colnames(GSE228632)[1] <- 'ensembl_gene_id' 

##### Create a final dataframe with all the tumor data 
momentary_list <- list(GSE181157,GSE227832,GSE133499,T_all,GSE228632) 
final <- momentary_list %>%  reduce(full_join,by = 'ensembl_gene_id') %>% drop_na()

# create merged dataframe of all the data take into consideration 
write.csv(final,file = '../Post_manipulation/Tumor_dataframe.csv',row.names = F)
 