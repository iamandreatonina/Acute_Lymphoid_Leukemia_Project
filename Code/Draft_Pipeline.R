# BiocManager::install('sva')
library(sva)
# library(tidyverse)
library(dplyr)
library(readxl)
library(ggplot2)
library(stringr) 
library(magrittr)
library(umap)
library(edgeR)
#BiocManager::install("DESeq2")
library(DESeq2)
#install.packages('factoextra')
#install.packages('cluster')
library(factoextra)
library(cluster)
library(RColorBrewer)
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(tibble)
library(plotly)
# library(rstudioapi)


`%nin%` <- Negate(`%in%`)

# Set Session to source file location and then:
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd("../data/Datasets/Post_manipulation")

# 88 controls and 705 tumor samples (now 1082)

##### Upload human specific genes and Data ####

Human_genes <- readxl::read_xlsx('Human-specific.xlsx')
#Load datasets 
Tumor <- read.csv('Tumor_dataframe.csv',sep =',',header = T) # needs to be unziped 
Control <- read.csv('Controls_merged.csv',sep = ',',header = T)

##### Batch effect correction ####
# We found out a duplicated ensembl_gene_id due to the fact there isn't a 1 to 1 mapping from ensembl_gene_id and hugo_symbols
# so we are going to eliminate the less informative one 

duplicato <- Tumor$ensembl_gene_id[duplicated(Tumor$ensembl_gene_id)]
duplicato2 <- Control$ensembl_gene_id[duplicated(Control$ensembl_gene_id)] # <- empty! Perfect!
#sum <- Tumor %>% dplyr::filter(Tumor$ensembl_gene_id == duplicato) 
# rowSums(sum[2:641]) # the first one is the most informative so we use distinct()
Tumor <- distinct(Tumor,ensembl_gene_id,.keep_all =T )

Tumor_2 <- as.matrix(sapply(Tumor[2:1083], as.numeric))
Control_2 <- as.matrix(sapply(Control[2:89], as.numeric))

PreCombat_control_df <- tidyr::gather(as.data.frame(Control_2),key = 'sample',value = 'read_number')
PreCombat_tumor_df_subset <- tidyr::gather(as.data.frame(Tumor_2)[1:20],key = 'sample',value = 'read_number')

jpeg(filename = '../images/control_Pre_Combat_boxplot.jpeg')
ggplot() +
  geom_boxplot(colour = 'darkgreen',fill = 'olivedrab',alpha = 0.7, mapping = aes(sample,read_number+1),data = PreCombat_control_df)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_y_log10()
dev.off()

# in this case there are too many samples so we are gonging to plot 20 samples instead of 640 
jpeg(filename = '../images/tumor_Pre_Combat_boxplot_subset.jpeg')
ggplot() +
  geom_boxplot(colour = 'darkred',fill = 'indianred',alpha = 0.5, mapping = aes(sample,read_number+1),data = PreCombat_tumor_df_subset, width = 0.5)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_y_log10()
dev.off()

#creation of batch for tumor and control, so creation of the vectors for batch separation, so check the data!!
batch_tumor <- c(rep(1,173),rep(2,321),rep(3,38),rep(4,108),rep(5,65),rep(6,377))
batch_control <- c(rep(1,20),rep(2,10),rep(3,40), rep(4,18))

# application of Combat-Seq and creation of adjusted dataframes 
tumor_adjusted <- as.data.frame(ComBat_seq(Tumor_2,batch = batch_tumor,group = NULL))
control_adjusted <- as.data.frame(ComBat_seq(Control_2, batch = batch_control, group = NULL))

# adding the ensembl_gene_id column 
colnames(control_adjusted)
control_adjusted <- add_column(control_adjusted,'ensembl_gene_id' =Control$ensembl_gene_id, .before = 'TU0049_CD4_HC')
colnames(tumor_adjusted)
tumor_adjusted <- add_column(tumor_adjusted, 'ensembl_gene_id' = Tumor$ensembl_gene_id, .before = 'GSM5491718_16.001')

##### Normalization with edgeR package ####

# We use TMM method , which is a normalization method intra and inter-sample and we create CPM matrices 

# Let`s check how many human specific genes we have in our dataset
HSgenes_tumor <- tumor_adjusted %>% dplyr::filter(tumor_adjusted$ensembl_gene_id %in% Human_genes$`Ensembl ID`)  
HSgenes_control <- control_adjusted %>% dplyr::filter(control_adjusted$ensembl_gene_id %in% Human_genes$`Ensembl ID`) 

# Previously the result is that of 873 human specific genes there are present 498 in control 
# But after adding controls ans tumors samples ( 88 control total and 705 tumors total ) in control we have 478 and tumor 603 
# After adding 377 tumor samples results: 478 controls and 598 tumors

# We don't have a precise threshold to eliminate the low expressed genes, so we know that DESeq2 set a 
#threshold based on the expression of our data by doing result()

# set the dataframe more easier for us to use 
control_adjusted1 <- control_adjusted %>% column_to_rownames('ensembl_gene_id')
tumor_adjusted1 <- tumor_adjusted %>% column_to_rownames('ensembl_gene_id')

#let's have a look at the data using a boxplot, we will make a comparison after the normalization
Pre_control_df <- tidyr::gather(control_adjusted1,key = 'sample',value = 'read_number')
Pre_tumor_df_subset <- tidyr::gather(tumor_adjusted1[1:20],key = 'sample',value = 'read_number')

jpeg(filename = '../images/control_Pre_TMM_boxplot.jpeg')
ggplot() +
  geom_boxplot(colour = 'darkgreen',fill = 'olivedrab',alpha = 0.7, mapping = aes(sample,read_number+1),data = Pre_control_df)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_y_log10()
dev.off()

# in this case there are too many samples so we are gonging to plot 20 samples instead of 1083 
jpeg(filename = '../images/tumor_Pre_TMM_boxplot_subset.jpeg')
ggplot() +
  geom_boxplot(colour = 'darkred',fill = 'indianred',alpha = 0.5, mapping = aes(sample,read_number+1),data = Pre_tumor_df_subset, width = 0.5)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_y_log10()
dev.off()

# first we find the DGEList object 
edge_c_control <- DGEList(counts = control_adjusted1)
edge_c_tumor <- DGEList(counts = tumor_adjusted1)

# normalize with the edgeR package using the TMM method, which apply and inter and intra normalization of the data, both for the controls and the tumor 
edge_n_control <- calcNormFactors(edge_c_control,method = 'TMM') 
edge_n_tumor <- calcNormFactors(edge_c_tumor,method = 'TMM')

# from that we create a CPM table with normalized expression values 
CPM_control <- as.data.frame(round(cpm(edge_n_control),2)) 
CPM_tumor <-  as.data.frame(round(cpm(edge_n_tumor),2))

# For what we created this datraframes?? 
CPM_control_df <- tidyr::gather(CPM_control,key = 'sample',value = 'CPM')
CPM_tumor_df <- tidyr::gather(CPM_tumor,key = 'sample',value = 'CPM')

# Save the CPM table 
# write.csv(CPM_tumor,file = 'CPM_Tumor_dataframe.csv',row.names = T)
# write.csv(CPM_control,file = 'CPM_Control_dataframe.csv',row.names = T)

CPM_tumor_df_toplot <- tidyr::gather(CPM_tumor[1:20],key = 'sample',value = 'CPM')

jpeg(filename = '../images/control_TMM_boxplot.jpeg')
ggplot() +
         geom_boxplot(colour = 'darkgreen',fill = 'olivedrab',alpha = 0.7, mapping = aes(sample,CPM+1),data = CPM_control_df)+
         theme_bw()+
         scale_x_discrete(guide = guide_axis(angle = 90))+
         scale_y_log10()
dev.off()

# in this case there are too many samples so we are gonging to plot 20 samples instead of 640 
jpeg(filename = '../images/tumor_TMM_boxplot_subset.jpeg')
ggplot() +
  geom_boxplot(colour = 'darkred',fill = 'indianred',alpha = 0.5, mapping = aes(sample,CPM+1),data = CPM_tumor_df_toplot, width = 0.5)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  scale_y_log10()
dev.off()

# -> noramlizzazione non un granche', come mai? Forse perche' adesso sono tanti? 

##### Differential gene expression analysis ####

total_adjusted <- merge(control_adjusted,tumor_adjusted,by='ensembl_gene_id')
total_adjusted1 <- total_adjusted %>% column_to_rownames('ensembl_gene_id')

# creating a dataframe containing the info on the samples, this is needed to be able to perform the DGE, we set the conditions of the samples as healty (H) and tumoral (T)
info_sample_1<-data.frame("sample"=colnames(total_adjusted1))
rownames(info_sample_1)<-info_sample_1$sample
info_sample_2<-as.data.frame(str_split(string=info_sample_1$sample, pattern="_", simplify=T)) #? serve? 
colnames(info_sample_2)<-c("condition","replicate")
info_samples<-cbind(info_sample_1, info_sample_2[1:2])
info_samples$condition<-c(rep("H",88),rep("T",1082)) # which are 88 healthy and 705 tumors 
info_samples$replicate<-c(rep(1,1170))


# Eliminating the samples and control base on the median and tharshold 

# new dataframe for the info samples to use to filter the 0 from our dataset
# The conditions are AT = adult tumor, PH = pediatric healthy, AH = adult healthy, PT = pediatric tumor
# info_samples_new_cond<-info_samples
# 
# ## Setting controls 
# info_samples_new_cond$condition<-"AT" 
# #GSE84445 -> first 20 AH 
# info_samples_new_cond$condition[1:20]<-"AH"
# #GSE227832 -> 10 PH
# info_samples_new_cond$condition[21:30]<-"PH"
# #GSE139073 -> 40 AH 
# info_samples_new_cond$condition[31:70] <- 'AH'
# #GSE115736 -> 18 unkonw +UK
# 
# 
# info_samples_new_cond$condition[31:562] <- 'PT'
# info_samples_new_cond$condition[info_samples_new_cond$sample %in% c('CMUTALLS4','T59','T91','T89','T87','T82','T81','T74','T59','T112','T102','SIHTALLS32','SIHTALLS25','SIHTALLS12','H301TALLS3','H301TALLS13','H301TALLS11','CMUTALLS9','CMUTALLS13','T67','T77','T103')] <-"PT"

# let's filter the dataset and setting the threshold definition
median_thr<-5
cond_tresh<-0.5
# filter_vec<-apply(total_adjusted1, 1, function(y) max(by(y,info_samples_new_cond$condition, function(x) median(x>=median_thr))) )
filter_vec<-apply(total_adjusted1, 1, function(y) max(by(y,info_samples$condition, function(x) median(x>=median_thr))) )
filter_counts_df <- total_adjusted1[filter_vec>=cond_tresh,]


# Now we can create the DGEList object
# edge_c_total <- DGEList(counts = total_adjusted1, group=info_samples$condition, samples=info_samples, genes=total_adjusted1)
# edge_n_total <- calcNormFactors(edge_c_total,method = 'TMM')


edge_c_total <- DGEList(counts = filter_counts_df, group=info_samples$condition, samples=info_samples, genes=filter_counts_df)
edge_n_total <- calcNormFactors(edge_c_total,method = 'TMM')

# We create the cpm table
cpm_table <-as.data.frame(round(cpm(edge_n_total),2)) # the library size is scaled by the normalization factor

# Here we define the experimental design matrix, we build a model with no intercept also we have two varaibles, one for each condition 
# 1 for control and 2 for tumor 
design <- model.matrix(~0+group, data = edge_n_total$samples)
colnames(design) <- levels(edge_n_total$samples$group)
rownames(design) <- edge_n_total$samples$sample

# Calculate dispersion and fit the result with edgeR (necessary for differential expression analysis)
edge_d_total <- estimateDisp(edge_n_total,design)

# Fit the data we model the data using a negative binomial distribution
edge_f<-glmQLFit(edge_d_total, design)

# Definition of the contrast (conditions to be compared)
contro <- makeContrasts("T-H", levels=design)

# Fit the model with generalized linear models
edge_t <- glmQLFTest(edge_f,contrast=contro)

# edge_t contains the results of the DE analysis
# -> we can extract the data using the function topTags -> extract the top20, using a cut off and sorting by fold-change
# -> we get the top 20 DE genes
#  We sort for the fold change
DEGs <- as.data.frame(topTags(edge_t,n=12969,p.value = 0.01,sort.by = "logFC"))

# We add a new column to the DEGs dataframe called class.
# Used to express the values of the fold change of the transcripts.
# The selection is based on the log fold change ratio (>1.5 for up-regulated genes and < (-1.5) for down-regulated genes)
# and a log CPM (>1 for both cases). From the contingency table of our DEGs we can see that the up regulated genes
# correspond to the 3.7% of the total and the down regulated are the 16% of the total.

DEGs$class <- '='
DEGs$class[which(DEGs$logCPM > 1 & DEGs$logFC > 1.5)] = '+'
DEGs$class[which(DEGs$logCPM > 1 & DEGs$logFC < (-1.5))] = '-'
DEGs <- DEGs[order(DEGs$logFC, decreasing = T),] # we order based on the fold change

table(DEGs$class)
# AFTER ADDING NEW DATA
#  - 736   + 2161   = 7555

#after adding the new 377 B-ALl samples 
# -  706  + 2012   = 7528
  
# Let`s check how many human specific genes we have in the up regulated and down regulated genes

DEGs_Hsgenes <- DEGs %>% dplyr::filter(rownames(DEGs) %in% Human_genes$`Ensembl ID`)
Up_HSgenes <- DEGs[DEGs$class=='+',] %>% dplyr::filter(rownames(DEGs[DEGs$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HSgenes <- DEGs[DEGs$class=='-',] %>% dplyr::filter(rownames(DEGs[DEGs$class=='-',]) %in% Human_genes$`Ensembl ID`) 

table(DEGs_Hsgenes$class)
# AFTER ADDING 12 - and 85 +  163 =HS DEGs

#after adding the new 377 B-ALl samples 
# -14   + 72  = 160
     


# Display the results using a volcano plot (x-axes: log FoldChange, y-axes: inverse function of the p-value).
# We can see the most significant DEGs colored in green, which are genes that surpass a threshold set on both the p-value
# and the Fold Change.

######### vulcano all gene ####
jpeg(filename = '../images/Vulcano_plot_DEGs.jpeg')
input_df<-DEGs
xlabel<- "log2 FC control vs case"
ylabel<-"-log10 p-value"
par(fig=c(0,1, 0,1), mar=c(4,4,1,2), mgp=c(2, 0.75,0))
plot(DEGs$logFC,-log(DEGs$PValue, base=10), xlab=xlabel,ylab = ylabel, col=ifelse(DEGs$class=="=", "grey70", "olivedrab4"), pch=20, frame.plot=TRUE, cex=0.8, main="Volcano plot") %>% 
abline(v = 0, lty = 2, col="grey20")
dev.off()

######### vulcano hs ####
jpeg(filename = '../images/Vulcano_plot_DEGsHS.jpeg')
input_df<-DEGs_Hsgenes
xlabel<- "log2 FC control vs case"
ylabel<-"-log10 p-value"
par(fig=c(0,1, 0,1), mar=c(4,4,1,2), mgp=c(2, 0.75,0))
plot(DEGs_Hsgenes$logFC,-log(DEGs_Hsgenes$PValue, base=10), xlab=xlabel,ylab = ylabel, col=ifelse(DEGs_Hsgenes$class=="=", "grey70", "olivedrab4"), pch=20, frame.plot=TRUE, cex=0.8, main="Volcano plot") %>% 
  abline(v = 0, lty = 2, col="grey20")
dev.off()

# We can also represent the genes using a heatmap.
# A clustering process is operated. We plot only up or down expressed genes using data from both the normalized CPM and
# the log transformation of the CPM table.

############ heatmap all genes ####

col <- rep('chartreuse4', 1170) # 1170 numebr of sampels 
col[which(info_samples$condition == 'T')] <- 'burlywood3' 
pal <- c('blue','white','red')
pal <- colorRampPalette(pal)(1170) # 1170 numebr of sampels 
DEGs_selected <- DEGs %>% dplyr::filter(DEGs$class != '=')
jpeg(filename = '../images/Heatmap_plot_DEGs.jpeg')
heatmap(as.matrix(cpm_table[which(rownames(cpm_table) %in% rownames(DEGs_selected)),]),ColSideColors = col, cexCol = 0.5,margins = c(4,4), col = pal, cexRow = 0.2)
dev.off()

#for improving the clusterization we set the cpm table as logarithmic

cpm_table_log <- as.data.frame(round(log10(cpm(edge_n_total)+1),2))
# write.csv(cpm_table_log,file ='Acute_lymphoide_leukemia_cpm_log_expression_table.csv',row.names = T )

# moment1 <- read.csv('Acute_lymphoide_leukemia_cpm_log_expression_table.csv')
# 
# moment1 <- moment1 %>% column_to_rownames('gene_names')
# write.csv(moment1,file ='Acute_lymphoide_leukemia_cpm_log_expression_table.csv')


jpeg(filename = '../images/Heatmap_plot_DEGs_log.jpeg')
heatmap(as.matrix(cpm_table_log[which(rownames(cpm_table_log) %in% rownames(DEGs_selected)),]),ColSideColors = col, cexCol = 0.5,margins = c(4,4), col = pal, cexRow = 0.2)
dev.off()

############ heatmap human specific ####

col <- rep('chartreuse4', 1170) # 1170 numebr of sampels 
col[which(info_samples$condition == 'T')] <- 'burlywood3' 
pal <- c('blue','white','red')
pal <- colorRampPalette(pal)(1170) # 1170 numebr of sampels 
DEGs_selected <- DEGs_Hsgenes %>% dplyr::filter(DEGs_Hsgenes$class != '=')
jpeg(filename = '../images/Heatmap_plot_DEGsHS.jpeg')
heatmap(as.matrix(cpm_table[which(rownames(cpm_table) %in% rownames(DEGs_selected)),]),ColSideColors = col, cexCol = 0.5,margins = c(4,4), col = pal, cexRow = 0.2)
dev.off()

#for improving the clusterization we set the cpm table as logarithmic
cpm_table_log <- as.data.frame(round(log10(cpm(edge_n_total)+1),2))
jpeg(filename = '../images/Heatmap_plot_DEGs_logHS.jpeg')
heatmap(as.matrix(cpm_table_log[which(rownames(cpm_table_log) %in% rownames(DEGs_selected)),]),ColSideColors = col, cexCol = 0.5,margins = c(4,4), col = pal, cexRow = 0.2)
dev.off()



##### PCA analysis ####
Diff_expressed <- DEGs[which(DEGs$class != '='),]
PCA_cpm_log_nonHS <- cpm_table_log[which(rownames(cpm_table_log) %in% rownames(Diff_expressed)),] # all genes Diff expressed
PCA_cpm_log <- cpm_table_log[which(rownames(cpm_table_log) %in% rownames(DEGs_selected)),] #only HS not diff expressed

PCA_cpm_log_woHS <- cpm_table_log[which(rownames(cpm_table_log) %in% rownames(Diff_expressed) & rownames(cpm_table_log) %nin% rownames(DEGs_selected)),] # Diff expressed genes without HS

# # we need to have both for the columns and the row a variance different from zero (because divide for the varaince )
PCA_cpm_log_filtered<-PCA_cpm_log[,which(apply(PCA_cpm_log, 2, var) != 0)]
PCA_cpm_log_filtered<- PCA_cpm_log_filtered[which(apply(PCA_cpm_log_filtered, 1, var) != 0),]
color<- c(rep('darkgreen',88),rep('indianred',1082))

PCA_cpm_log_nonHS_filtered <- PCA_cpm_log_nonHS[,which(apply(PCA_cpm_log_nonHS, 2, var) != 0)]
PCA_cpm_log_nonHS_filtered <- PCA_cpm_log_nonHS[which(apply(PCA_cpm_log_nonHS, 1, var) != 0),]


PCA_cpm_log_woHS_filtered <- PCA_cpm_log_woHS[,which(apply(PCA_cpm_log_woHS, 2, var) != 0)]
PCA_cpm_log_woHS_filtered <- PCA_cpm_log_woHS[which(apply(PCA_cpm_log_woHS, 1, var) != 0),]

#### PCA plot HS ####
data.PC_HS <- prcomp(t(PCA_cpm_log_filtered),scale. = T)
jpeg(filename = '../images/PCA_plot_DEGs_log_HS.jpeg')
plot(data.PC_HS$x[,1:2],col=color,pch = 19) 
dev.off()

#### PCA plot all genes####
data.PC_nonHS <- prcomp(t(PCA_cpm_log_nonHS_filtered),scale. = T)
plot(data.PC_nonHS$x[,1:2],col=color,pch = 19) 

#### PCA without HS ####
data.PC_woHS <- prcomp(t(PCA_cpm_log_woHS_filtered),scale. = T)
plot(data.PC_woHS$x[,1:2],col=color,pch=19)


#### PCA plot of tumor only ####
# only HS
data.PC.tumor_HS <- prcomp(t(PCA_cpm_log_filtered[89:1170]),scale. = T )
# jpeg(filename = '../images/PCA_plot_DEGs_log_tumor.jpeg')
plot(data.PC.tumor_HS$x[,1:2],pch = 19,col = c(rep('indianred',1082)))
# dev.off()

# all together 
data.PC_nonHS_tumor <- prcomp(t(PCA_cpm_log_nonHS_filtered[89:1170]),scale. = T )
plot(data.PC_nonHS_tumor$x[,1:2],pch = 19,col = c(rep('indianred',1082)))

#without HS 
data.PC_woHS_tumor <- prcomp(t(PCA_cpm_log_woHS_filtered[89:1170]),scale. = T )
plot(data.PC_woHS_tumor$x[,1:2],pch=19, col = c(rep('indianred',1082)))


#### PCA plot of control only ####

# only HS
data.PC.control_HS <- prcomp(t(PCA_cpm_log_filtered[1:88]),scale. = T )
# jpeg(filename = '../images/PCA_plot_DEGs_log_tumor.jpeg')
plot(data.PC.control_HS$x[,1:2],pch = 19, col= c(rep('darkgreen',88)))
# dev.off()

# all together 
data.PC_nonHS_control <- prcomp(t(PCA_cpm_log_nonHS_filtered[1:88]),scale. = T )
plot(data.PC_nonHS_control$x[,1:2],pch = 19,col= c(rep('darkgreen',88)))

#without HS
data.PC_woHS_control <- prcomp(t(PCA_cpm_log_woHS_filtered[1:88]),scale. = T )
plot(data.PC_woHS_control$x[,1:2],col= c(rep('darkgreen',88)),pch=19)



##### Partitioning around medoids, need to also to install cmake ####

#### I SEGUENTI SERVONO?####

# # for control-tumor -> 2 clusters as seen from graphs under
# fviz_nbclust(data.PC$x,FUNcluster = cluster::pam,k.max = 10)
# fviz_nbclust(data.PC$x,FUNcluster = cluster::pam,k.max = 10,method = 'gap_stat')+ theme_classic()
# fviz_nbclust(data.PC$x,FUNcluster = cluster::pam,k.max = 10, method = "wss")
# 
# # For subtypes of tumors
# fviz_nbclust(data.PC.tumor$x,FUNcluster = cluster::pam,k.max = 15)
# fviz_nbclust(data.PC.tumor$x,FUNcluster = cluster::pam,k.max = 15,method = 'gap_stat')+ theme_classic()
# fviz_nbclust(data.PC.tumor$x,FUNcluster = cluster::pam,k.max = 15, method = "wss")
# 
# # for control-tumor -> 2 clusters as seen from graphs under no HS
# fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 10)
# fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 10,method = 'gap_stat')+ theme_classic()
# fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 10, method = "wss")
# 
# # For subtypes of tumors no HS
# fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 15)
# fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 15,method = 'gap_stat')+ theme_classic()
# fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 15, method = "wss")


##### Questi servono mi sa #######
# # PAM on control-tumor
# # pam1<-eclust(data.PC_HS$x, "pam", k=9)
# pam1_nonHS <- eclust(data.PC_nonHS$x,'pam',k=9)
# 
# # PAM tumors subtypes
# pam2<-eclust(data.PC.tumor_HS$x, "pam", k=8)
# pam2_nonHS <- eclust(data.PC_nonHS_tumor$x,'pam',k=8)
# 
# 
# ##### hierarchical clustering ####
# #fanno schifo  
# 
# #calculate distances between observations and create a simple dendogram
# # dm <- dist(data.PC$x)
# # hc <- hclust(dm,method = 'average')
# # plot(hc,hang =-1)
# # rect.hclust(hc,k=2,border = 'red')
# # clust.vec.2<-cutree(hc,k=2)
# # fviz_cluster(list(data=data.PC$x, cluster=clust.vec.2))
# # 
# # # same for tumor subtypes
# # dm2 <- dist(data.PC.tumor$x)
# # hc2 <- hclust(dm2,method = 'average')
# # plot(hc2,hang =-1)
# # rect.hclust(hc2,k=2,border = 'red')
# # clust.vec<-cutree(hc2,k=8)
# # fviz_cluster(list(data=data.PC.tumor$x, cluster=clust.vec))
# # 
# # # clusters <- mutate(cpm_table_log[31:670],cluster =clust.vec
# 
# clusterino_pam2<-as.data.frame((pam2$clustering))
# components<-data.PC.tumor_HS[["x"]]
# components<-data.frame(components)
# components<-cbind(components, clusterino_pam2)
# components$PC2<- -components$PC2 # non capisco perchè mettiam positivo pc2
# 
# 
# components$`(pam2$clustering)` <- as.factor(components$`(pam2$clustering)`)
# 
# fig<-plot_ly(components, x=~PC1, y=~PC2, color=clusterino_pam2$`(pam2$clustering)`,colors=c('cadetblue1', 'red', 'chartreuse3','blueviolet','blue4','darkgoldenrod2','darksalmon','seagreen4') ,type='scatter',mode='markers') #  %>%
# # layout(legend = list(title = list(text = 'color')))
# fig
# # 
# 
# fig2<-plot_ly(components, x=~PC1, y=~PC2, z=~PC3,color=components$`(pam2$clustering)`,colors=brewer.pal(n = 8, name = "RdBu") ,mode='markers',marker = list(size = 4)) #  %>%
# fig2
# 
# 
# clusterino_pam2_nonHS<-as.data.frame((pam2_nonHS$clustering))
# components_nonHS<-data.PC_nonHS_tumor$x
# # components<-data.frame(components)
# components_nonHS<-cbind(components_nonHS, clusterino_pam2_nonHS)
# components_nonHS$PC2<- -components_nonHS$PC2
# components_nonHS$`(pam2_nonHS$clustering)`<- as.factor(components_nonHS$`(pam2_nonHS$clustering)`)
# 
# 
# fig_nonHS<-plot_ly(components_nonHS, x=~PC1, y=~PC2, color=clusterino_pam2_nonHS$`(pam2_nonHS$clustering)`,colors=c('cadetblue1', 'red', 
#                                                                                             'chartreuse3','blueviolet','blue4','darkgoldenrod2','darksalmon','seagreen4') ,type='scatter',mode='markers') #  %>%
# # layout(legend = list(title = list(text = 'color')))
# 
# fig_nonHS
# 
# fig2_nonHS<-plot_ly(components_nonHS, x=~PC1, y=~PC2, z=~PC3, color =components_nonHS$`(pam2_nonHS$clustering)`, colors=brewer.pal(n = 8, name = "RdBu") ,mode='markers',marker = list(size = 4)) #  %>%
# fig2_nonHS
# layout(legend = list(title = list(text = 'color')))

# GSE181157, 173 pediatric tumor samples
# GSE227832, 321 tumor samples, all pediatric from the metadata
# GSE133499, 38 tumor samples pediatric
# T_all, 107 tumor samples, both pediatric and adults
# GSE228632, 65 samples all tumor samples, bone marrow, all pediatric

# SO, we have: 619 sample pediatrici
# and 85 sample adult (only T!)


### Setting pediatric and adults 
# clusterino_pam2<-data.frame(data.PC.tumor_HS)


''' Allora si può mettere la stesso tipo di dataframe come usando dati pca$x gia presenti, 
l unica cosa è cambiare il vaore di pC2 a positivo, però non capisco il motivo, per averlo più comodo nel plot??  '''


clusterino_pam2 <- as.data.frame(data.PC.tumor_HS$x)
clusterino_pam2$PC2 <- - clusterino_pam2$PC2
clusterino_pam2$type <- 'pediatric'
clusterino_pam2$age <- 0
metadataGSE181157<-  readxl::read_xlsx('../Tumors/GSE181157_SampleMetadata.xlsx')
metadataGSE181157$`DFCI ID`<- rownames(clusterino_pam2)[1:173]# HERE I KEEP PRE-B AND PRE-T TYPES!
for (row in 1:dim(metadataGSE181157)[1]){
  Age <- metadataGSE181157$`Age at Dx (years)`[row]  
  if (Age<=18){
    clusterino_pam2$type[rownames(clusterino_pam2) == metadataGSE181157$`DFCI ID`[row]] <- 'pediatric'
    clusterino_pam2$age[rownames(clusterino_pam2) == metadataGSE181157$`DFCI ID`[row]] <- round(Age,2)} else{
      clusterino_pam2$type[rownames(clusterino_pam2) == metadataGSE181157$`DFCI ID`[row]] <- 'adult' 
      clusterino_pam2$age[rownames(clusterino_pam2) == metadataGSE181157$`DFCI ID`[row]] <- round(Age,2)}
}

metadata_choort_7_8 <- readxl::read_xlsx('../Tumors/Metadata_choort_7_8.xlsx', skip=1)

for (row in 1:dim(metadata_choort_7_8)[1]){  Age <- metadata_choort_7_8$`Age (year)`[row]
            if (Age=='Not available'){    
              clusterino_pam2$type[rownames(clusterino_pam2) == metadata_choort_7_8$ID[row]] <- 'Unknown'
              clusterino_pam2$age[rownames(clusterino_pam2) == metadata_choort_7_8$ID[row]] <- -1} 
            else if (as.numeric(Age)<=18){    
              clusterino_pam2$type[rownames(clusterino_pam2) == metadata_choort_7_8$ID[row]] <- 'pediatric'
              clusterino_pam2$age[rownames(clusterino_pam2) == metadata_choort_7_8$ID[row]] <- as.numeric(Age)
              } 
            else {
              clusterino_pam2$type[rownames(clusterino_pam2) == metadata_choort_7_8$ID[row]] <- 'adult'
              clusterino_pam2$age[rownames(clusterino_pam2) == metadata_choort_7_8$ID[row]] <- as.numeric(Age)
}}






components2 <- as.data.frame(data.PC.tumor$x)
components2<-cbind(components2, clusterino_pam2)
components2$PC2 <- -components2$PC2
fig3<-plot_ly(components2, x=~PC1, y=~PC2, color=clusterino_pam2$type,colors=c('red2', 'blue4') ,type='scatter',mode='markers')

fig3

fig4<-plot_ly(components2, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2$type,colors=c('darkred', 'blue4') ,mode='markers')
fig4

components2_nonHS <- as.data.frame(data.PC_nonHG_tumor$x)
components2_nonHS<-cbind(components2_nonHS, clusterino_pam2_nonHS)
components2_nonHS$PC2 <- -components2_nonHS$PC2
fig3_nonHS<-plot_ly(components2_nonHS, x=~PC1, y=~PC2, color=clusterino_pam2_nonHS$type,colors=c('red2', 'blue4') ,type='scatter',mode='markers')
fig3_nonHS

fig4_nonHS<-plot_ly(components2_nonHS, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2_nonHS$type,colors=c('darkred', 'blue4') ,mode='markers')
fig4_nonHS

# #### QUESTO DIREI CHE NON SERVE ####
# setwd("~/Desktop/magistrale_Qcb/3master_QCB_first_semester_second_year/biological_data_mining_blanzieri/Laboratory_Biological_Data_Mining/Dataset/GSE181157")
# 
#   
# setwd("~/Desktop/magistrale_Qcb/3master_QCB_first_semester_second_year/biological_data_mining_blanzieri/Laboratory_Biological_Data_Mining/Datasets_finals")

# metadata$`Final Risk`<- replace(metadata$`Final Risk`, metadata$`Final Risk` == 'Not Available', NA) 

# clusterino_pam2$risk <- 'Not Available'
# clusterino_pam2$risk[rownames(clusterino_pam2) %in% metadata$`DFCI ID`] <- metadata$`Final Risk`
# componet3 <- data.PC.tumor$x
# componet3 <- cbind(componet3,clusterino_pam2)
# componet3$PC2 <- -componet3$PC2 
# 
# fig5<-plot_ly(componet3, x=~PC1, y=~PC2, color=clusterino_pam2$risk,colors=c('red2', 'blue4') ,type='scatter',mode='markers')
# fig5
# 
# fig6<-plot_ly(componet3, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2$risk,colors=c('darkred', 'blue4') ,mode='markers')
# fig6
# 

# setwd("~/Desktop/magistrale_Qcb/3master_QCB_first_semester_second_year/biological_data_mining_blanzieri/Laboratory_Biological_Data_Mining/Dataset/GSE181157")
# 
# 
# setwd("~/Desktop/magistrale_Qcb/3master_QCB_first_semester_second_year/biological_data_mining_blanzieri/Laboratory_Biological_Data_Mining/Datasets_finals")
# 
# # metadata$`Final Risk`<- replace(metadata$`Final Risk`, metadata$`Final Risk` == 'Not Available', NA) 
# 
# clusterino_pam2_nonHS$risk <- 'Not Available'
# clusterino_pam2_nonHS$risk[rownames(clusterino_pam2_nonHS) %in% metadata_nonHS$`DFCI ID`] <- metadata_nonHS$`Final Risk`
# componet3_nonHS <- data.PC_nonHG_tumor$x
# componet3_nonHS <- cbind(componet3_nonHS,clusterino_pam2_nonHS)
# componet3_nonHS$PC2 <- -componet3_nonHS$PC2 
# 
# fig5_nonHS<-plot_ly(componet3_nonHS, x=~PC1, y=~PC2, color=clusterino_pam2_nonHS$risk,colors=c('red2', 'blue4') ,type='scatter',mode='markers')
# fig5_nonHS
# 
# fig6_nonHS<-plot_ly(componet3, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2$risk,colors=c('darkred', 'blue4') ,mode='markers')
# fig6_nonHS




#### DA QUI SERVE E FACCIAMO SOTTOTIPI CLASSIFICATION ####


# GSE181157,GSE227832,GSE133499,T_all,GSE228632
# GSE181157 is all Pre B and Pre T
# GSE227832 and GSE228632 is mixed
# GSE133499 is mixed


# NB SET WORKING DIRECTORY TO WHERE YPU HAVE THE METADATA INFO #

# set all the subtypes as B cells
clusterino_pam2$Cell_type <- 'B'

# GSE181157 are the first 173 samples
metadata<-  readxl::read_xlsx('../Tumors/GSE181157_SampleMetadata.xlsx')
#metadata_nonHS<-  readxl::read_xlsx('GSE181157_SampleMetadata.xlsx')
#metadata_nonHS$`DFCI ID` <- rownames(clusterino_pam2_nonHS)[1:173]
metadata$`DFCI ID`<- rownames(clusterino_pam2)[1:173]
# HERE I KEEP PRE-B AND PRE-T TYPES!

for (row in 1:dim(metadata)[1]){
  Diagnosis <- metadata$Diagnosis[row]
  if (Diagnosis=='9836/3 - Pre-B ALL'){
    clusterino_pam2$Cell_type[rownames(clusterino_pam2) == metadata$`DFCI ID`[row]] <- 'PreB'
  } else{
    clusterino_pam2$Cell_type[rownames(clusterino_pam2) == metadata$`DFCI ID`[row]] <- 'PreT'
  }
}



# Dataset T-ALL is all T subtype
clusterino_pam2$Cell_type[533:640] <- 'T' #by literature T_ALL is of only T cells 

# GSE227832 and GSE228632 is mixed
metadata_GSE227832<-  readxl::read_xlsx('../Tumors/Metadata_GSE227832_GSE228632/NEW_Metadata_GSE227832_GSE228632.xlsx',skip=1, col_names=T) # non ho questa cartella, manco su github

table(metadata_GSE227832$`Subtype at ALL diagnosis`)

for (row in 1:dim(metadata_GSE227832)[1]){
  Diagnosis<-metadata_GSE227832$`Subtype at ALL diagnosis`[row]
  if (Diagnosis =='T-ALL'){
    clusterino_pam2$Cell_type[rownames(clusterino_pam2) == metadata_GSE227832$public_id[row]] <- 'T'
  } else{
    clusterino_pam2$Cell_type[rownames(clusterino_pam2) == metadata_GSE227832$public_id[row]] <- "B"
  }
}

# Dataset GSE133499 is mixed 
metadata_GSE133499<-  readxl::read_xlsx('../Tumors/Metadata_GSE133499.xlsx', col_names=T) # the column IPT stand for immunephenotype -> ergo: commoni (Btype), pre B and T or unknown

for (row in 1:dim(metadata_GSE133499)[1]){
  Diagnosis<-metadata_GSE133499$IPT[row]
  if (Diagnosis=="T"){
    clusterino_pam2$Cell_type[rownames(clusterino_pam2) == metadata_GSE133499$`Anonym ID`[row]] <- "T"
  } else if (Diagnosis=="B"){
    clusterino_pam2$Cell_type[rownames(clusterino_pam2) == metadata_GSE133499$`Anonym ID`[row]] <- "B"
  } else if (Diagnosis=="pre-B"){
    clusterino_pam2$Cell_type[rownames(clusterino_pam2) == metadata_GSE133499$`Anonym ID`[row]] <- "PreB"
  } else{
    clusterino_pam2$Cell_type[rownames(clusterino_pam2) == metadata_GSE133499$`Anonym ID`[row]] <- "Unknown"
  }
  
}


table(clusterino_pam2$Cell_type)
#  B-ALL   Pre-B   Pre-T   T-ALL Unknown 
#  367     144      37     133      24 

componet4 <- data.PC.tumor$x
componet4 <- cbind(componet4,clusterino_pam2)
componet4$PC2 <- -componet4$PC2

fig7<-plot_ly(componet4, x=~PC1, y=~PC2, color=clusterino_pam2$Cell_type,colors=c('red', 'blue','green','orange', 'grey') ,type='scatter',mode='markers')
fig7

fig8<-plot_ly(componet4, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2$Cell_type,colors=c('darkred', 'blue4','green','orange', 'grey'),type='scatter3d', size=10, mode='markers')
fig8

#### All genes in only tumor ####

data.PC_nonHG_tumor

componet7 <- data.PC_nonHG_tumor$x
componet7 <- cbind(componet7,clusterino_pam2)
componet7$PC2 <- -componet7$PC2

fig7<-plot_ly(componet7, x=~PC1, y=~PC2, color=clusterino_pam2$Cell_type,colors=c('red', 'blue','green','orange', 'grey') ,type='scatter',mode='markers')
fig7

fig8<-plot_ly(componet7, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2$Cell_type,colors=c('darkred', 'blue4','green','orange', 'grey'),type='scatter3d', size=10, mode='markers')
fig8


# 
# 
# #HS
# clusterino_pam1<-as.data.frame((pam1$clustering))
# clusterino_pam1$C_T <- "Tumor"
# clusterino_pam1$C_T[rownames(clusterino_pam1) %in% c("X817_T","X845_B","X845_T","X858_B","X858_T","X867_B","X867_T","X899_B","X899_T","X817_B","TU0049_CD4_HC","TU0049_CD8_HC",
#                                                      "TU0051_CD4_HC","TU0051_CD8_HC","TU0054_CD4_HC","TU0054_CD8_HC","XT0130_CD4_HC","XT0130_CD8_HC","XT0133_CD4_HC","XT0133_CD8_HC",
#                                                      "XT0108_CD4_HC","XT0108_CD8_HC","XT0115_CD4_HC","XT0115_CD8_HC","XT0127_CD4_HC","XT0127_CD8_HC","XT0131_CD4_HC","XT0131_CD8_HC",
#                                                      "XT0141_CD4_HC","XT0141_CD8_HC")] <- 'Control'
# clusterino_pam1$type <- 'pediatric'
# clusterino_pam1$type[563:670] <- 'adult'
# clusterino_pam1$type[rownames(clusterino_pam1) %in% c('CMUTALLS4','T59','T91','T89','T87','T82','T81','T74','T59','T112','T102','SIHTALLS32','SIHTALLS25','SIHTALLS12','H301TALLS3','H301TALLS13','H301TALLS11','CMUTALLS9','CMUTALLS13','T67','T77','T103')] <- 'pediatric'
# clusterino_pam1$type[11:30] <- 'adult'
# clusterino_pam1$risk <- 'Not Available'
# clusterino_pam1$risk[rownames(clusterino_pam1) %in% metadata$`DFCI ID`] <- metadata$`Final Risk`
# clusterino_pam1$Cell_type <- 'Unkown'
# clusterino_pam1$Cell_type[1:30] <- 'Control'
# clusterino_pam1$Cell_type[563:670] <- 'T Cell' #by letaruet of only T cells
# clusterino_pam1$Cell_type[rownames(clusterino_pam1) %in% metadata$`DFCI ID`] <- metadata$Diagnosis
# component5 <- data.PC$x
# component5<- cbind(component5,clusterino_pam1)
# component5$PC2 <- -component5$PC2
# 
# ### sistemare
# fig9<-plot_ly(component5, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam1$Cell_type,colors=brewer.pal(n = 4, name = "RdBu"),  symbol = clusterino_pam1$type, symbols = c('diamond','circle'), mode='markers',marker = list(size = 4))
# fig9
# 
# #write.csv(component5,file='ML_HS.csv',row.names = T)
# #### Non HS
# clusterino_pam1_nonHS<-as.data.frame((pam1_nonHS$clustering))
# clusterino_pam1_nonHS$C_T <- "Tumor"
# clusterino_pam1_nonHS$C_T[rownames(clusterino_pam1_nonHS) %in% c("X817_T","X845_B","X845_T","X858_B","X858_T","X867_B","X867_T","X899_B","X899_T","X817_B","TU0049_CD4_HC","TU0049_CD8_HC",
#                                                                  "TU0051_CD4_HC","TU0051_CD8_HC","TU0054_CD4_HC","TU0054_CD8_HC","XT0130_CD4_HC","XT0130_CD8_HC","XT0133_CD4_HC","XT0133_CD8_HC",
#                                                                  "XT0108_CD4_HC","XT0108_CD8_HC","XT0115_CD4_HC","XT0115_CD8_HC","XT0127_CD4_HC","XT0127_CD8_HC","XT0131_CD4_HC","XT0131_CD8_HC",
#                                                                  "XT0141_CD4_HC","XT0141_CD8_HC")] <- 'Control'
# clusterino_pam1_nonHS$type <- 'pediatric'
# clusterino_pam1_nonHS$type[563:670] <- 'adult'
# clusterino_pam1_nonHS$type[rownames(clusterino_pam1_nonHS) %in% c('CMUTALLS4','T59','T91','T89','T87','T82','T81','T74','T59','T112','T102','SIHTALLS32','SIHTALLS25','SIHTALLS12','H301TALLS3','H301TALLS13','H301TALLS11','CMUTALLS9','CMUTALLS13','T67','T77','T103')] <- 'pediatric'
# clusterino_pam1_nonHS$type[11:30] <- 'adult'
# clusterino_pam1_nonHS$risk <- 'Not Available'
# clusterino_pam1_nonHS$risk[rownames(clusterino_pam1_nonHS) %in% metadata$`DFCI ID`] <- metadata$`Final Risk`
# clusterino_pam1_nonHS$Cell_type <- 'Unkown'
# clusterino_pam1_nonHS$Cell_type[1:30] <- 'Control'
# clusterino_pam1_nonHS$Cell_type[563:670] <- 'T Cell' #by letaruet of only T cells
# clusterino_pam1_nonHS$Cell_type[rownames(clusterino_pam1_nonHS) %in% metadata$`DFCI ID`] <- metadata$Diagnosis
# component6 <- data.PC_nonHG$x
# component6<- cbind(component6,clusterino_pam2)
# component6$PC2 <- -component6$PC2
# #write.csv(component6,file='ML_nonHS.csv',row.names = T)
# 
# ##### sistemare
# fig9_nonHS<-plot_ly(component6, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2$Cell_type,colors=brewer.pal(n = 4, name = "RdBu"), mode='markers',marker = list(size = 4))
# fig9_nonHS


####### Categories of HS genes found UP and Down between Tumor and Control
## Above we defined Diff_expressed <- DEGs[which(DEGs$class != '='),]

Categories_HS <- Human_genes %>% dplyr::filter(Human_genes$`Ensembl ID` %in% rownames(Diff_expressed ))
table(Categories_HS$`General Mechanism of Origin`)
barplot(table(Categories_HS$`General Mechanism of Origin`), horiz=T, names.arg=c("1","2","3","4","5","6","7","8","9"),col=brewer.pal(n = 9, name = "Set1"))
legend(x="topright", inset=.02,y.intersp = 1,title="Mechanism of origin", legend=c("de novo origin","amplification","loss","sequence alteration","structure alteration","undefined feature","lost in chimpanzee","new non-coding gene","regulatory region alteration"), fill=brewer.pal(n = 9, name = "Set1"), cex=0.7)

slices<-c(1,63,4,21,10,11,3,1,8)
pct<-round(slices/sum(slices)*100)
lbs<-c("","","","","","","","","")
lbs<-paste(lbs,pct)
lbs<-paste(lbs,"%",sep=" ")
pie(slices, lbs,col=brewer.pal(n = 9, name = "Set1"))
legend(x="topright", inset=.02,y.intersp = 1,title="Mechanism of origin", legend=c("de novo origin","amplification","loss","sequence alteration","structure alteration","human specific","lost in chimpanzee","new non-coding gene","regulatory region alteration"), fill=brewer.pal(n = 9, name = "Set1"), cex=0.7)

#### GSEA Tumor-Control DEGS ####

ensmebl <- useMart(biomart = 'ensembl',dataset = 'hsapiens_gene_ensembl')

############ GSEA non HS ####
convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = rownames(DEGs), mart = ensmebl)

DEGs_2 <- tibble::rownames_to_column(DEGs, var = 'ensembl_gene_id')

DEGs_merge_convert <- merge(DEGs_2, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

DEGs_merge_convert<-DEGs_merge_convert[which(!is.na(DEGs_merge_convert$entrezgene_id)),] 

DEGs_merge_convert<-DEGs_merge_convert[-which(duplicated(DEGs_merge_convert$entrezgene_id)),]

Up_DEGs_merge_convert <- DEGs_merge_convert %>% dplyr::filter(DEGs_merge_convert$class == '+')

Down_DEGs_merge_convert <- DEGs_merge_convert %>% dplyr::filter(DEGs_merge_convert$class == '-')
#### Gene Ontology enrichment analysis (biological process) UP####
ego_BP_UP <- enrichGO(gene = Up_DEGs_merge_convert$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05)

barplot(ego_BP_UP, showCategory = 15)

dotplot(ego_BP_UP, showCategory=15)

heatplot(ego_BP_UP, showCategory = 2)
#### We perform KEGG enrichment analysis. UP####
#We use function enrichWP to retrieve the list of genes from the wiki pathways, we can see which pathways are more expressed.

eWP_BP_UP <- enrichWP(gene =Up_DEGs_merge_convert$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.05, qvalueCutoff = 0.1 )

head(eWP_BP_UP,10)

# Resutls 
# WP2446 WP2446                                  Retinoblastoma gene in cancer   40/1038 89/8421 1.7e-14  1.1e-11 9.9e-12
# WP5218 WP5218 Extrafollicular and follicular B cell activation by SARS CoV 2   32/1038 76/8421 6.2e-11  2.0e-08 1.8e-08
# WP2849 WP2849                        Hematopoietic stem cell differentiation   27/1038 58/8421 1.2e-10  2.5e-08 2.3e-08
# WP3937 WP3937                        Microglia pathogen phagocytosis pathway   21/1038 40/8421 8.6e-10  1.4e-07 1.2e-07
# WP45     WP45                                     G1 to S cell cycle control   26/1038 64/8421 9.7e-09  1.2e-06 1.1e-06
# WP3945 WP3945                             TYROBP causal network in microglia   25/1038 62/8421 2.2e-08  2.4e-06 2.2e-06
# WP466   WP466                                                DNA replication   19/1038 42/8421 1.2e-07  1.1e-05 1.0e-05
# WP23     WP23                              B cell receptor signaling pathway   30/1038 92/8421 2.6e-07  2.1e-05 1.9e-05
# WP4016 WP4016                    DNA IR damage and cellular response via ATR   27/1038 81/8421 6.3e-07  4.4e-05 4.0e-05
# WP4884 WP4884      Pathogenesis of SARS CoV 2 mediated by nsp9 nsp10 complex   12/1038 21/8421 1.2e-06  7.4e-05 6.7e-05



#### Gene Ontology enrichment analysis (biological process) Down ####
ego_BP_DW <- enrichGO(gene = Down_DEGs_merge_convert$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05)

barplot(ego_BP_DW, showCategory = 15)

dotplot(ego_BP_DW, showCategory=15)

heatplot(ego_BP_DW, showCategory = 2)

#### We perform KEGG enrichment analysis. Down ####

eWP_BP_DW <- enrichWP(gene =Down_DEGs_merge_convert$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.05, qvalueCutoff = 0.1 )

# head(eWP_BP_DW@result[["Description"]],10)
head(eWP_BP_DW,10)

# Result
# WP306   WP306                                              Focal adhesion    35/405 199/8421 1.5e-11  7.5e-09 6.7e-09
# WP2911 WP2911                 miRNA targets in ECM and membrane receptors    12/405  38/8421 1.1e-07  2.9e-05 2.5e-05
# WP5055 WP5055                                          Burn wound healing    16/405  75/8421 3.8e-07  6.5e-05 5.7e-05
# WP2572 WP2572            Primary focal segmental glomerulosclerosis FSGS     14/405  72/8421 6.6e-06  8.3e-04 7.4e-04
# WP5284 WP5284 Cell interactions of the pancreatic cancer microenvironment     8/405  26/8421 1.9e-05  2.0e-03 1.7e-03
# WP3888 WP3888                                      VEGFA VEGFR2 signaling    40/405 433/8421 4.7e-05  4.0e-03 3.5e-03
# WP5348 WP5348                     11p11 2 copy number variation syndrome     11/405  56/8421 5.8e-05  4.2e-03 3.8e-03
# WP3874 WP3874                 Canonical and non canonical TGF B signaling     6/405  17/8421 9.4e-05  6.0e-03 5.3e-03
# WP185   WP185                             Integrin mediated cell adhesion    14/405 102/8421 3.5e-04  2.0e-02 1.7e-02
# WP453   WP453                               Inflammatory response pathway     7/405  30/8421 4.4e-04  2.2e-02 2.0e-02

##### GSEA Tumor-Control HS ####

convert_HS <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = rownames(DEGs_Hsgenes), mart = ensmebl)

DEGs_Hsgenes_2 <- rownames_to_column(DEGs_Hsgenes, var = 'ensembl_gene_id')

DEGs_merge_convert_HS <- merge(DEGs_Hsgenes_2, convert_HS, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

DEGs_merge_convert_HS<-DEGs_merge_convert_HS[which(!is.na(DEGs_merge_convert_HS$entrezgene_id)),] 

DEGs_merge_convert_HS<-DEGs_merge_convert_HS[-which(duplicated(DEGs_merge_convert_HS$entrezgene_id)),]

Up_DEGs_merge_convert_HS<- DEGs_merge_convert_HS %>% dplyr::filter(DEGs_merge_convert_HS$class == '+')

Down_DEGs_merge_convert_HS <- DEGs_merge_convert_HS %>% dplyr::filter(DEGs_merge_convert_HS$class == '-')

#### Gene Ontology enrichment analysis (biological process) HS UP ####

ego_BP_UP_HS <- enrichGO(gene = Up_DEGs_merge_convert_HS$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.2, qvalueCutoff =0.2)

barplot(ego_BP_UP_HS)

dotplot(ego_BP_UP_HS, showCategory=15)

heatplot(ego_BP_UP_HS, showCategory = 2)

#### We perform KEGG enrichment analysis. HS UP ####

# pathway annotation will be done after the gene expansion, due to lack of informations, to low numbers of terms
eWP_BP_UP_HS <- enrichWP(gene =Up_DEGs_merge_convert_HS$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 1, qvalueCutoff = 0.2)

head(eWP_BP_UP_HS,10)

#Result
# WP2431 WP2431                                        Spinal cord injury      5/34 120/8421 0.00011     0.01 0.0099 
# WP5092 WP5092 Interactions of natural killer cells in pancreatic cancer      3/34  28/8421 0.00018     0.01 0.0099


##### Gene Ontology enrichment analysis (biological process) HS Down ####

ego_BP_DW_HS <- enrichGO(gene = Down_DEGs_merge_convert_HS$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.2, qvalueCutoff = 0.2)

barplot(ego_BP_DW_HS)

dotplot(ego_BP_DW_HS)

heatplot(ego_BP_DW_HS, showCategory = 2)


#### DEG subtype vs subtype without controls ####

# creating a dataframe containing the info on the samples, this is needed to be able to perform the DGE

# clusterino_pam2 contains in Cell type columns the info on the subtypes
info_subtypes<-clusterino_pam2
info_subtypes$sample<-rownames(info_subtypes)
info_subtypes<-info_subtypes[info_subtypes$Cell_type!="Unknown",]
tumors_subtype<-tumor_adjusted1[colnames(tumor_adjusted1) %in% rownames(info_subtypes)]


edge_c_subtypes <- DGEList(counts = tumors_subtype, group=info_subtypes$Cell_type, samples=info_subtypes, genes=tumors_subtype)
edge_n_subtypes <- calcNormFactors(edge_c_subtypes,method = 'TMM')
# We create the cpm table
cpm_table_subtypes <-as.data.frame(round(cpm(edge_n_subtypes),2)) # the library size is scaled by the normalization factor

# Here we define the experimental design matrix, we build a model with no intercept also we have two varaibles, one for each condition 
# 1 for control and 2 for tumor 
design_subtype <- model.matrix(~0+group, data = edge_n_subtypes$samples, contrast.arg=list(group='contr.sum'))
colnames(design_subtype) <- levels(edge_n_subtypes$samples$group)
rownames(design_subtype) <- edge_n_subtypes$samples$sample


colnames(design_subtype)<-c("B","PreB", "PreT", "T")


# Calculate dispersion and fit the result with edgeR (necessary for differential expression analysis)
edge_d_subtype <- estimateDisp(edge_n_subtypes,design_subtype)

# Fit the data we model the data using a negative binomial distribution
edge_f_subtype<-glmQLFit(edge_d_subtype, design_subtype)

####PreB vs all other subtype####

# Definition of the contrast (conditions to be compared)
contro_subtype_PreB <- makeContrasts("PreB-(PreT+T+B)/3", levels=design_subtype)
#contro_subtype[,"PreB-PreT-T"]<-c(1, -0.5, -0.5)

# Fit the model with generalized linear models
edge_t_subtype_PreB <- glmQLFTest(edge_f_subtype,contrast=contro_subtype_PreB)

# edge_t contains the results of the DE analysis
# -> we can extract the data using the function topTags -> extract the top20, using a cut off and sorting by fold-change
# -> we get the top 20 DE genes
#  We sort for the fold change

# the n are the number of row. this case 21376
DEGs_subtype_PreB <- as.data.frame(topTags(edge_t_subtype_PreB,n=21376,p.value = 0.01,sort.by = "logFC"))

# We add a new column to the DEGs dataframe called class.
# Used to express the values of the fold change of the transcripts.
# The selection is based on the log fold change ratio (>1.5 for up-regulated genes and < (-1.5) for down-regulated genes)
# and a log CPM (>1 for both cases). From the contingency table of our DEGs we can see that the up regulated genes
# correspond to the 3.7% of the total and the down regulated are the 16% of the total.

DEGs_subtype_PreB$class <- '='
DEGs_subtype_PreB$class[which(DEGs_subtype_PreB$logCPM > 1 & DEGs_subtype_PreB$logFC > 1.5)] = '+'
DEGs_subtype_PreB$class[which(DEGs_subtype_PreB$logCPM > 1 & DEGs_subtype_PreB$logFC < (-1.5))] = '-'
DEGs_subtype_PreB <- DEGs_subtype_PreB[order(DEGs_subtype_PreB$logFC, decreasing = T),] # we order based on the fold change

table(DEGs_subtype_PreB$class)
# with new samples -327    +276    =8394 
   

# Let`s check how many human specific genes we have in the up regulated and down regulated genes in Pre B type

DEGs_subtype_PreB_HS <- DEGs_subtype_PreB %>% dplyr::filter(rownames(DEGs_subtype_PreB) %in% Human_genes$`Ensembl ID`)
Up_HS_PreB <- DEGs_subtype_PreB[DEGs_subtype_PreB$class=='+',] %>% dplyr::filter(rownames(DEGs_subtype_PreB[DEGs_subtype_PreB$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HS_PreB<- DEGs_subtype_PreB[DEGs_subtype_PreB$class=='-',] %>% dplyr::filter(rownames(DEGs_subtype_PreB[DEGs_subtype_PreB$class=='-',]) %in% Human_genes$`Ensembl ID`) 

table(DEGs_subtype_PreB_HS$class)
# with new samples -24   +8   =237
     


#### Subtype PreT vs all the others ####

# Definition of the contrast (conditions to be compared)
contro_subtype_PT <- makeContrasts("PreT-(PreB+T+B)/3", levels=design_subtype)
#contro_subtype[,"PreB-PreT-T"]<-c(1, -0.5, -0.5)

# Fit the model with generalized linear models
edge_t_subtype_PT<- glmQLFTest(edge_f_subtype,contrast=contro_subtype_PT)

# edge_t contains the results of the DE analysis
# -> we can extract the data using the function topTags -> extract the top20, using a cut off and sorting by fold-change
# -> we get the top 20 DE genes
#  We sort for the fold change

# the n are the number of row. this case 21376
DEGs_subtype_PT <- as.data.frame(topTags(edge_t_subtype_PT,n=21376,p.value = 0.01,sort.by = "logFC"))

# We add a new column to the DEGs dataframe called class.
# Used to express the values of the fold change of the transcripts.
# The selection is based on the log fold change ratio (>1.5 for up-regulated genes and < (-1.5) for down-regulated genes)
# and a log CPM (>1 for both cases). From the contingency table of our DEGs we can see that the up regulated genes
# correspond to the 3.7% of the total and the down regulated are the 16% of the total.

DEGs_subtype_PT$class <- '='
DEGs_subtype_PT$class[which(DEGs_subtype_PT$logCPM > 1 & DEGs_subtype_PT$logFC > 1.5)] = '+'
DEGs_subtype_PT$class[which(DEGs_subtype_PT$logCPM > 1 & DEGs_subtype_PT$logFC < (-1.5))] = '-'
DEGs_subtype_PT <- DEGs_subtype_PT[order(DEGs_subtype_PT$logFC, decreasing = T),] # we order based on the fold change

table(DEGs_subtype_PT$class)
# with new samples  - 868    + 371   = 6809
    

# Let`s check how many human specific genes we have in the up regulated and down regulated genes in Pre B type
DEGs_subtype_PT_HS <- DEGs_subtype_PT %>% dplyr::filter(rownames(DEGs_subtype_PT) %in% Human_genes$`Ensembl ID`)
Up_HS_PreT <- DEGs_subtype_PT[DEGs_subtype_PT$class=='+',] %>% dplyr::filter(rownames(DEGs_subtype_PT[DEGs_subtype_PT$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HS_PreT<- DEGs_subtype_PT[DEGs_subtype_PT$class=='-',] %>% dplyr::filter(rownames(DEGs_subtype_PT[DEGs_subtype_PT$class=='-',]) %in% Human_genes$`Ensembl ID`) 

table(DEGs_subtype_PT_HS$class)
# with new samples   - 45  + 14   = 186 HS DEGs
    

#### T subtype vs all other subtypes ####

# Definition of the contrast (conditions to be compared)
contro_subtype_T <- makeContrasts("T-(PreB+PreT+B)/3", levels=design_subtype)
#contro_subtype[,"PreB-PreT-T"]<-c(1, -0.5, -0.5)

# Fit the model with generalized linear models
edge_t_subtype_T<- glmQLFTest(edge_f_subtype,contrast=contro_subtype_T)

# edge_t contains the results of the DE analysis
# -> we can extract the data using the function topTags -> extract the top20, using a cut off and sorting by fold-change
# -> we get the top 20 DE genes
#  We sort for the fold change

# the n are the number of row. this case 21376

DEGs_subtype_T <- as.data.frame(topTags(edge_t_subtype_T,n=21376,p.value = 0.01,sort.by = "logFC"))

# We add a new column to the DEGs dataframe called class.
# Used to express the values of the fold change of the transcripts.
# The selection is based on the log fold change ratio (>1.5 for up-regulated genes and < (-1.5) for down-regulated genes)
# and a log CPM (>1 for both cases). From the contingency table of our DEGs we can see that the up regulated genes
# correspond to the 3.7% of the total and the down regulated are the 16% of the total.

DEGs_subtype_T$class <- '='
DEGs_subtype_T$class[which(DEGs_subtype_T$logCPM > 1 & DEGs_subtype_T$logFC > 1.5)] = '+'
DEGs_subtype_T$class[which(DEGs_subtype_T$logCPM > 1 & DEGs_subtype_T$logFC < (-1.5))] = '-'
DEGs_subtype_T <- DEGs_subtype_T[order(DEGs_subtype_T$logFC, decreasing = T),] # we order based on the fold change

table(DEGs_subtype_T$class)
# with new samples    - 16   + 285   = 3417
    

# Let`s check how many human specific genes we have in the up regulated and down regulated genes in Pre B type
DEGs_subtype_T_HS <- DEGs_subtype_T %>% dplyr::filter(rownames(DEGs_subtype_T) %in% Human_genes$`Ensembl ID`)
Up_HS_T <- DEGs_subtype_T[DEGs_subtype_T$class=='+',] %>% dplyr::filter(rownames(DEGs_subtype_T[DEGs_subtype_T$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HS_T<- DEGs_subtype_T[DEGs_subtype_T$class=='-',] %>% dplyr::filter(rownames(DEGs_subtype_T[DEGs_subtype_T$class=='-',]) %in% Human_genes$`Ensembl ID`) 

table(DEGs_subtype_T_HS$class)
# with new samples +34   = 111
 


#### B subtype vs all other subtypes ####

# Definition of the contrast (conditions to be compared)
contro_subtype_B <- makeContrasts("B-(PreB+PreT+T)/3", levels=design_subtype)
#contro_subtype[,"PreB-PreT-T"]<-c(1, -0.5, -0.5)

# Fit the model with generalized linear models
edge_t_subtype_B<- glmQLFTest(edge_f_subtype,contrast=contro_subtype_B)

# edge_t contains the results of the DE analysis
# -> we can extract the data using the function topTags -> extract the top20, using a cut off and sorting by fold-change
# -> we get the top 20 DE genes
#  We sort for the fold change

# the n are the number of row. this case 21376

DEGs_subtype_B <- as.data.frame(topTags(edge_t_subtype_B,n=21376,p.value = 0.01,sort.by = "logFC"))

# We add a new column to the DEGs dataframe called class.
# Used to express the values of the fold change of the transcripts.
# The selection is based on the log fold change ratio (>1.5 for up-regulated genes and < (-1.5) for down-regulated genes)
# and a log CPM (>1 for both cases). From the contingency table of our DEGs we can see that the up regulated genes
# correspond to the 3.7% of the total and the down regulated are the 16% of the total.

DEGs_subtype_B$class <- '='
DEGs_subtype_B$class[which(DEGs_subtype_B$logCPM > 1 & DEGs_subtype_B$logFC > 1.5)] = '+'
DEGs_subtype_B$class[which(DEGs_subtype_B$logCPM > 1 & DEGs_subtype_B$logFC < (-1.5))] = '-'
DEGs_subtype_B <- DEGs_subtype_B[order(DEGs_subtype_B$logFC, decreasing = T),] # we order based on the fold change

table(DEGs_subtype_B$class)
# +: 293  -:51   =: 6928 


# Let`s check how many human specific genes we have in the up regulated and down regulated genes in Pre B type
DEGs_subtype_B_HS <- DEGs_subtype_B %>% dplyr::filter(rownames(DEGs_subtype_B) %in% Human_genes$`Ensembl ID`)
Up_HS_B <- DEGs_subtype_B[DEGs_subtype_B$class=='+',] %>% dplyr::filter(rownames(DEGs_subtype_B[DEGs_subtype_B$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HS_B<- DEGs_subtype_B[DEGs_subtype_B$class=='-',] %>% dplyr::filter(rownames(DEGs_subtype_B[DEGs_subtype_B$class=='-',]) %in% Human_genes$`Ensembl ID`) 

table(DEGs_subtype_B_HS$class)
# with new samples -2 +31 = 233




#### COMPARISON BETWEEN THE HS IN THE 4 SUBTYPES ####

# NB: done in this way otherwise by using for example DEGs HS of T vs DEGs HS of PreT then can be some matches that are not supposed to be done 

## Up vs Up 
Up_HS_T_unique <- Up_HS_T[rownames(Up_HS_T) %nin% c(rownames(Up_HS_PreT),rownames(Up_HS_PreB),rownames(Up_HS_B)),]
table(Up_HS_T_unique$class) # 3+ 

Up_HS_PreT_unique <-  Up_HS_PreT[rownames(Up_HS_PreT) %nin% c(rownames(Up_HS_T),rownames(Up_HS_PreB),rownames(Up_HS_B)),]
table(Up_HS_PreT_unique$class) #9+ 

Up_HS_PreB_unique <- Up_HS_PreB[rownames(Up_HS_PreB) %nin% c(rownames(Up_HS_T),rownames(Up_HS_PreT),rownames(Up_HS_B)),]
table(Up_HS_PreB_unique$class) # 4+ 

Up_HS_B_unique <- Up_HS_B[rownames(Up_HS_B) %nin% c(rownames(Up_HS_PreT),rownames(Up_HS_PreB),rownames(Up_HS_T)),]
table(Up_HS_B_unique$class) # 3+ 

## Down vs Down
Down_HS_T_unique <- Down_HS_T[rownames(Down_HS_T) %nin% c(rownames(Down_HS_PreT),rownames(Down_HS_PreB),rownames(Down_HS_B)),]
table(Down_HS_T_unique$class) # 0! becasue we didn't have any down reg 

Down_HS_B_unique <- Down_HS_B[rownames(Down_HS_B) %nin% c(rownames(Down_HS_PreT),rownames(Down_HS_PreB),rownames(Down_HS_T)),]
table(Down_HS_B_unique$class) # 0! 

Down_HS_PreB_unique <-Down_HS_PreB[rownames(Down_HS_PreB) %nin% c(rownames(Down_HS_PreT),rownames(Down_HS_B),rownames(Down_HS_T)),]
table(Down_HS_PreB_unique$class) # 11 - 

Down_HS_PreT_unique <-Down_HS_PreT[rownames(Down_HS_PreT) %nin% c(rownames(Down_HS_PreB),rownames(Down_HS_B),rownames(Down_HS_T)),]
table(Down_HS_PreT_unique$class) # 34 - 

# Results of DEGS of HS without commons 

# Pre B with new samples 4+ 11-
# Pre T # with new samples 34- 9+ 
# T with new samples 3+ 0- 
# B with new samples 3+ 0- 


Total_HS_unqiue <- c(rownames(Up_HS_B_unique),rownames(Down_HS_B_unique),rownames(Up_HS_T_unique),rownames(Down_HS_T_unique),rownames(Up_HS_PreB_unique),rownames(Down_HS_PreB_unique),
                     rownames(Up_HS_PreT_unique),rownames(Down_HS_PreT_unique))



#### UMAP DEGs tumors only HS  #### 

cpm_table_subtypes_log <- as.data.frame(round(log10(cpm(edge_n_subtypes)+1),2))
Subtype_cpm_log_only_HS <- cpm_table_subtypes_log[which(rownames(cpm_table_subtypes_log) %in% Total_HS_unqiue),]
Subtype_cpm_log_only_HS_filtered <- Subtype_cpm_log_only_HS[,which(apply(Subtype_cpm_log_only_HS, 2, var) != 0)]
Subtype_cpm_log_only_HS_filtered <- Subtype_cpm_log_only_HS_filtered[which(apply(Subtype_cpm_log_only_HS_filtered, 1, var) != 0),]



# Perform UMAP dimensionality reduction
umap_result4 <- umap(t(Subtype_cpm_log_only_HS_filtered),
                     n_neighbors = 5,         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)



umap_result4 <- umap(t(Subtype_cpm_log_only_HS_filtered),
                     n_neighbors = 15,         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)

umap_result4 <- umap(t(Subtype_cpm_log_only_HS_filtered),
                    n_neighbors = sqrt(dim(t(Subtype_cpm_log_only_HS_filtered))[1]),         #or the square root of the rows 
                    min_dist = 0.1,          #
                    metric = "euclidean",    #you can change it
                    n_components = 3
)

clusterino_umap_subtypes <- clusterino_pam2[clusterino_pam2$Cell_type != 'Unknown',]
umap_df4 <- data.frame(umap_result4$layout,
                       Cell_type = clusterino_umap_subtypes$Cell_type)

colnames(umap_df4) <- c("umap_1","umap_2","umap_3", "Cell_type")
rownames(umap_df4) <-  rownames(t(Subtype_cpm_log_only_HS_filtered))
# Print UMAP
print(umap_result4)

#umap_df <- cbind(umap_df,disease_state)
# Create a 2D scatter plot using Plotly
fig2U <- plot_ly(umap_df4, 
                 x = ~umap_1, y = ~umap_2, z = ~umap_3,
                 color = umap_df4$Cell_type,
                 colors = c("blue","red","green","orange", "grey"),   
                 mode = 'markers',
                 size=10) %>% layout(title = 'Tumor subtypes HS, metric euclidian, neighbors = square')

# Display the 2D scatter plot
fig2U

#### COMPARISON BETWEEN THE DEGS not only HS IN THE 4 SUBTYPES ####

## subtype B 

Up_subtype_B_unique <- DEGs_subtype_B[DEGs_subtype_B$class == '+',] # 293 + 
Down_subtype_B_unique <- DEGs_subtype_B[DEGs_subtype_B$class == '-',]# 51 -

Up_subtype_B_unique <- Up_subtype_B_unique[rownames(Up_subtype_B_unique) %nin% c(rownames(DEGs_subtype_PreB[DEGs_subtype_PreB$class == '+',]),rownames(DEGs_subtype_PT[DEGs_subtype_PT$class == '+',]),
                                                                                      rownames(DEGs_subtype_T[DEGs_subtype_T$class == '+',])),]

Down_subtype_B_unique <- Down_subtype_B_unique[rownames(Down_subtype_B_unique) %nin% c(rownames(DEGs_subtype_PreB[DEGs_subtype_PreB$class == '-',]),rownames(DEGs_subtype_PT[DEGs_subtype_PT$class == '-',]),
                                                                                      rownames(DEGs_subtype_T[DEGs_subtype_T$class == '-',])),]
table(Up_subtype_B_unique$class) # 35+ 
table(Down_subtype_B_unique$class) # 4 -


## Subtype T 

Up_subtype_T_unique <- DEGs_subtype_T[DEGs_subtype_T$class == '+',] # 285 + 
Down_subtype_T_unique <- DEGs_subtype_T[DEGs_subtype_T$class == '-',]# 16 -

Up_subtype_T_unique <- Up_subtype_T_unique[rownames(Up_subtype_T_unique) %nin% c(rownames(DEGs_subtype_PreB[DEGs_subtype_PreB$class == '+',]),rownames(DEGs_subtype_PT[DEGs_subtype_PT$class == '+',]),
                                                                                 rownames(DEGs_subtype_B[DEGs_subtype_B$class == '+',])),]

Down_subtype_T_unique <- Down_subtype_T_unique[rownames(Down_subtype_T_unique) %nin% c(rownames(DEGs_subtype_PreB[DEGs_subtype_PreB$class == '-',]),rownames(DEGs_subtype_PT[DEGs_subtype_PT$class == '-',]),
                                                                                       rownames(DEGs_subtype_B[DEGs_subtype_B$class == '-',])),]
table(Up_subtype_T_unique$class) # 53 + 
table(Down_subtype_T_unique$class) # 11 -

## Subtype Pre B 

Up_subtype_PreB_unique <- DEGs_subtype_PreB[DEGs_subtype_PreB$class == '+',] # 276 + 
Down_subtype_PreB_unique <- DEGs_subtype_PreB[DEGs_subtype_PreB$class == '-',]# 327 -

Up_subtype_PreB_unique <- Up_subtype_PreB_unique[rownames(Up_subtype_PreB_unique) %nin% c(rownames(DEGs_subtype_T[DEGs_subtype_T$class == '+',]),rownames(DEGs_subtype_PT[DEGs_subtype_PT$class == '+',]),
                                                                                 rownames(DEGs_subtype_B[DEGs_subtype_B$class == '+',])),]

Down_subtype_PreB_unique <- Down_subtype_PreB_unique[rownames(Down_subtype_PreB_unique) %nin% c(rownames(DEGs_subtype_T[DEGs_subtype_T$class == '-',]),rownames(DEGs_subtype_PT[DEGs_subtype_PT$class == '-',]),
                                                                                       rownames(DEGs_subtype_B[DEGs_subtype_B$class == '-',])),]
table(Up_subtype_PreB_unique$class) # 127 + 
table(Down_subtype_PreB_unique$class) # 241 -

##Subtype Pre T

Up_subtype_PreT_unique <- DEGs_subtype_PT[DEGs_subtype_PT$class == '+',] # 371 + 
Down_subtype_PreT_unique <- DEGs_subtype_PT[DEGs_subtype_PT$class == '-',]# 868 -

Up_subtype_PreT_unique <- Up_subtype_PreT_unique[rownames(Up_subtype_PreT_unique) %nin% c(rownames(DEGs_subtype_T[DEGs_subtype_T$class == '+',]),rownames(DEGs_subtype_PreB[DEGs_subtype_PreB$class == '+',]),
                                                                                          rownames(DEGs_subtype_B[DEGs_subtype_B$class == '+',])),]

Down_subtype_PreT_unique <- Down_subtype_PreT_unique[rownames(Down_subtype_PreT_unique) %nin% c(rownames(DEGs_subtype_T[DEGs_subtype_T$class == '-',]),rownames(DEGs_subtype_PreB[DEGs_subtype_PreB$class == '-',]),
                                                                                                rownames(DEGs_subtype_B[DEGs_subtype_B$class == '-',])),]
table(Up_subtype_PreT_unique$class) # 293 + 
table(Down_subtype_PreT_unique$class) # 811 -


# Results of DEGS of the subtypes by taking in consideration all the genes without commons 

# Pre B with new samples 127+ 241-
# Pre T with new samples 811- 293+ 
# T with new samples 53+ 11- 
# B with new samples 35+ 4-


#### UMAP without HS genes and commons ####

Up_subtype_B_unique_withoutHS <- Up_subtype_B_unique[rownames(Up_subtype_B_unique) %nin% Human_genes$`Ensembl ID`,]
Down_subtype_B_unique_withoutHS <- Down_subtype_B_unique[rownames(Down_subtype_B_unique) %nin% Human_genes$`Ensembl ID`,]

Up_subtype_T_unique_withoutHS <-  Up_subtype_T_unique[rownames(Up_subtype_T_unique) %nin% Human_genes$`Ensembl ID`,]
Down_subtype_T_unique_withoutHS <- Down_subtype_T_unique[rownames(Down_subtype_T_unique) %nin% Human_genes$`Ensembl ID`,]

Up_subtype_PreB_unique_withoutHS <-  Up_subtype_PreB_unique[rownames(Up_subtype_PreB_unique) %nin% Human_genes$`Ensembl ID`,]
Down_subtype_PreB_unique_withoutHS <-  Down_subtype_PreB_unique[rownames(Down_subtype_PreB_unique) %nin% Human_genes$`Ensembl ID`,]

Up_subtype_PreT_unique_withoutHS <- Up_subtype_PreT_unique[rownames(Up_subtype_PreT_unique) %nin% Human_genes$`Ensembl ID`,]
Down_subtype_PreT_unique_withoutHS <- Down_subtype_PreT_unique[rownames(Down_subtype_PreT_unique) %nin% Human_genes$`Ensembl ID`,]
# DEGs_subtype_B_without_HS<-DEGs_subtype_B[DEGs_subtype_B$class != '=',]
# DEGs_subtype_B_without_HS<-DEGs_subtype_B_without_HS[rownames(DEGs_subtype_B_without_HS) %nin% Human_genes$`Ensembl ID`,]
# 
# DEGs_subtype_T_without_HS<-DEGs_subtype_T[DEGs_subtype_T$class != '=',]
# DEGs_subtype_T_without_HS<-DEGs_subtype_T_without_HS[rownames(DEGs_subtype_T_without_HS) %nin% Human_genes$`Ensembl ID`,]
# 
# DEGs_subtype_PreT_without_HS<-DEGs_subtype_PT[DEGs_subtype_PT$class != '=',]
# DEGs_subtype_PreT_without_HS<-DEGs_subtype_PreT_without_HS[rownames(DEGs_subtype_PreT_without_HS) %nin% Human_genes$`Ensembl ID`,]
# 
# DEGs_subtype_PreB_without_HS<-DEGs_subtype_PreB[DEGs_subtype_PreB$class != '=',]
# DEGs_subtype_PreB_without_HS<-DEGs_subtype_PreB_without_HS[rownames(DEGs_subtype_PreB_without_HS) %nin% Human_genes$`Ensembl ID`,]

Total_without_HS_unique <- c(rownames(Up_subtype_B_unique_withoutHS),rownames(Down_subtype_B_unique_withoutHS),rownames(Up_subtype_T_unique_withoutHS),rownames(Down_subtype_T_unique_withoutHS), 
                             rownames(Up_subtype_PreB_unique_withoutHS),rownames(Down_subtype_PreB_unique_withoutHS),rownames(Up_subtype_PreT_unique_withoutHS),rownames(Down_subtype_PreT_unique_withoutHS))


Subtype_cpm_log_without_HS <- cpm_table_subtypes_log[which(rownames(cpm_table_subtypes_log) %in% Total_without_HS_unique),]
Subtype_cpm_log_without_HS_filtered <- Subtype_cpm_log_without_HS[,which(apply(Subtype_cpm_log_without_HS, 2, var) != 0)]
Subtype_cpm_log_without_HS_filtered <- Subtype_cpm_log_without_HS_filtered[which(apply(Subtype_cpm_log_without_HS_filtered, 1, var) != 0),]



# Perform UMAP dimensionality reduction
umap_result5 <- umap(t(Subtype_cpm_log_without_HS_filtered),
                     n_neighbors = 5,         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)



umap_result5 <- umap(t(Subtype_cpm_log_without_HS_filtered),
                     n_neighbors = 15,         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)

umap_result5 <- umap(t(Subtype_cpm_log_without_HS_filtered),
                     n_neighbors = sqrt(dim(t(Subtype_cpm_log_without_HS_filtered))[1]),         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)

#clusterino_umap_subtypes <- clusterino_pam2[clusterino_pam2$Cell_type != 'Unknown',]
umap_df5 <- data.frame(umap_result5$layout,
                       Cell_type = clusterino_umap_subtypes$Cell_type)

colnames(umap_df5) <- c("umap_1","umap_2","umap_3", "Cell_type")
rownames(umap_df5) <-  rownames(t(Subtype_cpm_log_without_HS_filtered))
# Print UMAP
print(umap_result5)

#umap_df <- cbind(umap_df,disease_state)
# Create a 2D scatter plot using Plotly
fig3U <- plot_ly(umap_df5, 
                 x = ~umap_1, y = ~umap_2, z = ~umap_3,
                 color = umap_df4$Cell_type,
                 colors = c("blue","red","green","orange", "grey"),   
                 mode = 'markers',
                 size=10) %>%  layout(title = 'Tumor subtypes without HS and commons, metric euclidian, neighbors = square')

# Display the 2D scatter plot
fig3U

## FIle creation
#write.csv(DEGs_subtype_T,file = 'DEGs_subtype_T.csv',row.names = T, col.names = T)

#write.csv(DEGs_subtype_PT,file = 'DEGs_subtype_PT.csv',row.names = T, col.names = T)

#write.csv(DEGs_subtype_B,file = 'DEGs_subtype_B.csv',row.names = T, col.names = T)


#### GSEA subtypes HS ####

### Subtype pre-B ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_HS_PreB_unique), mart = ensmebl)
Up_HS_PreB_unique$ensembl_gene_id<-row.names(Up_HS_PreB_unique)

merged <- merge(Up_HS_PreB_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ###

up_PreB_hs <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(up_PreB_hs, showCategory = 15)

dotplot(up_PreB_hs, showCategory=15)

heatplot(up_PreB_hs, showCategory = 5)

###We perform KEGG enrichment analysis. HS UP ###

up_PreB_hs_kegg <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(up_PreB_hs_kegg,10)

#Results 
# WP2431 WP2431                            Spinal cord injury       1/1 120/8421  0.014    0.016     NA   1464     1
# WP5417 WP5417 Cell lineage map for neuronal differentiation       1/1 132/8421  0.016    0.016     NA   1464     1

###Gene Ontology enrichment analysis (biological process) HS Down ###

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Down_HS_PreB_unique), mart = ensmebl)
Down_HS_PreB_unique$ensembl_gene_id<-row.names(Down_HS_PreB_unique)

merged <- merge(Down_HS_PreB_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

down_PreB_hs <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(down_PreB_hs, showCategory = 15)

dotplot(down_PreB_hs, showCategory=15)

heatplot(down_PreB_hs, showCategory = 5)

### We perform KEGG enrichment analysis. HS Down ###

down_PreB_hs_wp <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(down_PreB_hs_wp, 10)

# WP1425 WP1425         Bone morphogenic protein signaling and regulation       1/3 12/8421 0.0043    0.048  0.011    657     1
# WP3874 WP3874               Canonical and non canonical TGF B signaling       1/3 17/8421 0.0060    0.048  0.011    657     1
# WP1591 WP1591                                         Heart development       1/3 47/8421 0.0167    0.048  0.011    657     1
# WP4917 WP4917                                 Proximal tubule transport       1/3 57/8421 0.0202    0.048  0.011    486     1
# WP5053 WP5053                 Development of ureteric collection system       1/3 60/8421 0.0212    0.048  0.011    657     1
# WP5352 WP5352             10q11 21q11 23 copy number variation syndrome       1/3 61/8421 0.0216    0.048  0.011    657     1
# WP474   WP474                                 Endochondral ossification       1/3 63/8421 0.0223    0.048  0.011    657     1
# WP4808 WP4808        Endochondral ossification with skeletal dysplasias       1/3 63/8421 0.0223    0.048  0.011    657     1
# WP5402 WP5402                            10q22q23 copy number variation       1/3 64/8421 0.0226    0.048  0.011    657     1
# WP2840 WP2840 Hair follicle development cytodifferentiation part 3 of 3       1/3 87/8421 0.0307    0.054  0.012    657     1


### Subtype pre-T ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_HS_PreT_unique), mart = ensmebl)
Up_HS_PreT_unique$ensembl_gene_id<-row.names(Up_HS_PreT_unique)

merged <- merge(Up_HS_PreT_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ###

up_PreT_hs <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(up_PreT_hs, showCategory = 15)

dotplot(up_PreT_hs, showCategory=15)

heatplot(up_PreT_hs, showCategory = 5)

### We perform KEGG enrichment analysis. HS UP ###

up_PreT_hs_kegg <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(up_PreT_hs_kegg,10)

#Results 
# WP1425 WP1425         Bone morphogenic protein signaling and regulation       1/1 12/8421 0.0014    0.015     NA    657     1
# WP3874 WP3874               Canonical and non canonical TGF B signaling       1/1 17/8421 0.0020    0.015     NA    657     1
# WP1591 WP1591                                         Heart development       1/1 47/8421 0.0056    0.015     NA    657     1
# WP5053 WP5053                 Development of ureteric collection system       1/1 60/8421 0.0071    0.015     NA    657     1
# WP5352 WP5352             10q11 21q11 23 copy number variation syndrome       1/1 61/8421 0.0072    0.015     NA    657     1
# WP474   WP474                                 Endochondral ossification       1/1 63/8421 0.0075    0.015     NA    657     1
# WP4808 WP4808        Endochondral ossification with skeletal dysplasias       1/1 63/8421 0.0075    0.015     NA    657     1
# WP5402 WP5402                            10q22q23 copy number variation       1/1 64/8421 0.0076    0.015     NA    657     1
# WP2840 WP2840 Hair follicle development cytodifferentiation part 3 of 3       1/1 87/8421 0.0103    0.017     NA    657     1
# WP5094 WP5094                                   Orexin receptor pathway       1/1 88/8421 0.0105    0.017     NA    657     1

### Gene Ontology enrichment analysis (biological process) HS Down ###

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Down_HS_PreT_unique), mart = ensmebl)
Down_HS_PreT_unique$ensembl_gene_id<-row.names(Down_HS_PreT_unique)

merged <- merge(Down_HS_PreT_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

down_PreT_hs <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(down_PreT_hs, showCategory = 15)

dotplot(down_PreT_hs, showCategory=15)

heatplot(down_PreT_hs, showCategory = 5)

### We perform KEGG enrichment analysis. HS Down ###

down_PreT_hs_wp <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(down_PreT_hs_wp, 10)

# P5269 WP5269 Genetic causes of porto sinusoidal vascular disease      2/12 37/8421 0.0012    0.069  0.051 653361/2969     2

### Subtype B ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_HS_B_unique), mart = ensmebl)
Up_HS_B_unique$ensembl_gene_id<-row.names(Up_HS_B_unique)

merged <- merge(Up_HS_B_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ###

up_B_hs <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(up_B_hs, showCategory = 15)

dotplot(up_B_hs, showCategory=15)

heatplot(up_B_hs, showCategory = 5)

### We perform KEGG enrichment analysis. HS UP ###

up_B_hs_kegg <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 1, qvalueCutoff = 1)
#Results 
# --> No gene can be mapped....
# --> Expected input gene ID: 4536,2932,56413,4708,3356,9550
# --> return NULL...

head(up_B_hs_kegg,10)
#NULL


#### For Down_HS_B_unique### 

# We didn't compute anything because it's empty!!


### Subtype T ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_HS_T_unique), mart = ensmebl)
Up_HS_T_unique$ensembl_gene_id<-row.names(Up_HS_T_unique)

merged <- merge(Up_HS_T_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ###

up_T_hs <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(up_T_hs, showCategory = 15)

dotplot(up_T_hs, showCategory=15)

heatplot(up_T_hs, showCategory = 5)

### We perform KEGG enrichment analysis. HS UP ###

up_T_hs_kegg <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

# --> No gene can be mapped....
# --> Expected input gene ID: 29933,8390,2878,4722,9551,539
# --> return NULL...

head(up_T_hs_kegg,10)

#Results 
# NULL

#### For Down_HS_T_unique### 

# We didn't compute anything because it's empty!!


#### GSEA subtypes all genes ####

### Subtype pre-B ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_subtype_PreB_unique), mart = ensmebl)
Up_subtype_PreB_unique$ensembl_gene_id<-row.names(Up_subtype_PreB_unique)

merged <- merge(Up_subtype_PreB_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ###

up_PreB <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(up_PreB, showCategory = 15)

dotplot(up_PreB, showCategory=15)

heatplot(up_PreB, showCategory = 5)

###We perform KEGG enrichment analysis. HS UP ###

up_PreB_kegg <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(up_PreB_hs_kegg,10)

#Results 
# WP2328 WP2328                            Allograft rejection      7/68  90/8421 7.186176e-06 0.0009049959 0.0008384876
# WP4217 WP4217                  Ebola virus infection in host      8/68 129/8421 8.340976e-06 0.0009049959 0.0008384876
# WP244   WP244               Alpha 6 beta 4 signaling pathway      4/68  33/8421 1.334009e-04 0.0096493294 0.0089401985
# WP3932 WP3932 Focal adhesion PI3K Akt mTOR signaling pathway      9/68 302/8421 6.669396e-04 0.0361814732 0.0335224904

###Gene Ontology enrichment analysis (biological process) HS Down ###

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Down_subtype_PreB_unique), mart = ensmebl)
Down_subtype_PreB_unique$ensembl_gene_id<-row.names(Down_subtype_PreB_unique)

merged <- merge(Down_subtype_PreB_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

down_PreB <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(down_PreB, showCategory = 15)

dotplot(down_PreB, showCategory=15)

heatplot(down_PreB, showCategory = 5)

### We perform KEGG enrichment analysis. HS Down ###

down_PreB_wp <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(down_PreB_wp, 10)

# Results:
# WP5072 WP5072                          Modulators of TCR signaling and T cell activation    10/125  61/8421 1.719441e-08 4.470546e-06 4.253353e-06
# WP69     WP69                                          T cell receptor signaling pathway    10/125  75/8421 1.326703e-07 1.724714e-05 1.640922e-05
# WP3863 WP3863 T cell antigen receptor TCR pathway during Staphylococcus aureus infection     9/125  62/8421 2.739663e-07 2.374375e-05 2.259020e-05
# WP4884 WP4884                  Pathogenesis of SARS CoV 2 mediated by nsp9 nsp10 complex     6/125  21/8421 4.288529e-07 2.787544e-05 2.652117e-05
# WP5098 WP5098                                               T cell activation SARS CoV 2    10/125  88/8421 6.140328e-07 3.192971e-05 3.037846e-05
# WP4585 WP4585                                      Cancer immunotherapy by PD 1 blockade     5/125  23/8421 1.806370e-05 7.827605e-04 7.447316e-04
# WP5130 WP5130                                          Th17 cell differentiation pathway     7/125  70/8421 7.406990e-05 2.751168e-03 2.617508e-03
# WP5115 WP5115                                Network map of SARS CoV 2 signaling pathway    12/125 218/8421 8.974658e-05 2.916764e-03 2.775059e-03
# WP4494 WP4494     Selective expression of chemokine receptors during T cell polarization     4/125  30/8421 9.404194e-04 2.716767e-02 2.584779e-02
# WP5142 WP5142      Calcium mediated T cell apoptosis involved in inclusion body myositis     3/125  20/8421 3.026474e-03 7.868831e-02 7.486540e-02


### Subtype pre-T ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_subtype_PreT_unique), mart = ensmebl)
Up_subtype_PreT_unique$ensembl_gene_id<-row.names(Up_subtype_PreT_unique)

merged <- merge(Up_subtype_PreT_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ###

up_PreT <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(up_PreT, showCategory = 15)

dotplot(up_PreT, showCategory=15)

heatplot(up_PreT, showCategory = 5)

### We perform KEGG enrichment analysis. HS UP ###

up_PreT_kegg <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(up_PreT_kegg,10)

#Results 
# WP69     WP69                                          T cell receptor signaling pathway    11/147  75/8421 5.926124e-08 9.884704e-06 9.330892e-06
# WP5072 WP5072                          Modulators of TCR signaling and T cell activation    10/147  61/8421 8.143396e-08 9.884704e-06 9.330892e-06
# WP3863 WP3863 T cell antigen receptor TCR pathway during Staphylococcus aureus infection    10/147  62/8421 9.565842e-08 9.884704e-06 9.330892e-06
# WP4884 WP4884                  Pathogenesis of SARS CoV 2 mediated by nsp9 nsp10 complex     6/147  21/8421 1.116928e-06 8.656190e-05 8.171208e-05
# WP5098 WP5098                                               T cell activation SARS CoV 2    10/147  88/8421 2.726557e-06 1.690465e-04 1.595753e-04
# WP5115 WP5115                                Network map of SARS CoV 2 signaling pathway    15/147 218/8421 5.579600e-06 2.882793e-04 2.721279e-04
# WP5130 WP5130                                          Th17 cell differentiation pathway     8/147  70/8421 2.695415e-05 1.193684e-03 1.126805e-03
# WP2583 WP2583                               T cell receptor and co stimulatory signaling     5/147  28/8421 1.076340e-04 4.170819e-03 3.937140e-03
# WP4585 WP4585                                      Cancer immunotherapy by PD 1 blockade     4/147  23/8421 6.095664e-04 2.099617e-02 1.981982e-02
# WP5078 WP5078                                     T cell modulation in pancreatic cancer     5/147  46/8421 1.166540e-03 3.616275e-02 3.413665e-02
# geneID Count
# WP69                     920/5588/9402/8631/925/29851/940/3932/919/916/27040    11
# WP5072                        5588/9402/921/925/53347/940/3932/919/916/27040    10
# WP3863                       920/5588/9402/959/925/29851/940/3932/5133/27040    10
# WP4884                                              920/914/925/3932/919/916     6
# WP5098                      920/5588/9402/64798/29851/940/3932/919/916/27040    10
# WP5115 3075/3675/920/5588/7850/7173/914/8631/814/925/64798/9332/3932/919/916    15
# WP5130                                 920/5588/2625/6097/3932/919/916/27040     8
# WP2583                                                925/5578/940/3932/5133     5
# WP4585                                                     925/3932/5133/916     4
# WP5078                                               959/29851/940/7293/5133     5

### Gene Ontology enrichment analysis (biological process) HS Down ###

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Down_subtype_PreT_unique), mart = ensmebl)
Down_subtype_PreT_unique$ensembl_gene_id<-row.names(Down_subtype_PreT_unique)

merged <- merge(Down_subtype_PreT_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

down_PreT <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(down_PreT, showCategory = 15)

dotplot(down_PreT, showCategory=15)

heatplot(down_PreT, showCategory = 5)

### We perform KEGG enrichment analysis. HS Down ###

down_PreT_wp <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(down_PreT_wp, 10)

# Results:
# WP2328 WP2328                                            Allograft rejection    21/382  90/8421 3.409704e-10 1.711671e-07 1.557696e-07
# WP5218 WP5218 Extrafollicular and follicular B cell activation by SARS CoV 2    17/382  76/8421 3.343284e-08 8.391643e-06 7.636765e-06
# WP23     WP23                              B cell receptor signaling pathway    16/382  92/8421 3.129697e-06 5.237027e-04 4.765925e-04
# WP5033 WP5033  Genes associated with the development of rheumatoid arthritis     7/382  18/8421 7.713134e-06 9.679983e-04 8.809211e-04
# WP4217 WP4217                                  Ebola virus infection in host    18/382 129/8421 1.901616e-05 1.909222e-03 1.737476e-03
# WP3932 WP3932                 Focal adhesion PI3K Akt mTOR signaling pathway    29/382 302/8421 1.012238e-04 8.469060e-03 7.707218e-03
# WP2877 WP2877                                     Vitamin D receptor pathway    20/382 187/8421 3.026050e-04 1.851162e-02 1.684639e-02
# WP5053 WP5053                      Development of ureteric collection system    10/382  60/8421 3.237148e-04 1.851162e-02 1.684639e-02
# WP4172 WP4172                                     PI3K Akt signaling pathway    30/382 341/8421 3.621089e-04 1.851162e-02 1.684639e-02
# WP4658 WP4658                                         Small cell lung cancer    13/382  96/8421 3.687573e-04 1.851162e-02 1.684639e-02

### Subtype B ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_subtype_B_unique), mart = ensmebl)
Up_subtype_B_unique$ensembl_gene_id<-row.names(Up_subtype_B_unique)

merged <- merge(Up_subtype_B_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ###

up_B <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(up_B, showCategory = 15)

dotplot(up_B, showCategory=15)

heatplot(up_B, showCategory = 5)

### We perform KEGG enrichment analysis. HS UP ###

up_B_kegg <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 1, qvalueCutoff = 1)

head(up_B_kegg,10)
#Resutls:
# WP3927 WP3927                 BMP signaling in eyelid development       1/4  20/8421 0.009467949 0.04648352 0.01631001   3625     1
# WP2203 WP2203 Thymic stromal lymphopoietin TSLP signaling pathway       1/4  40/8421 0.018868508 0.04648352 0.01631001  64109     1
# WP2324 WP2324                                    AGE RAGE pathway       1/4  66/8421 0.030989010 0.04648352 0.01631001   3625     1
# WP2369 WP2369                               Histone modifications       1/4  66/8421 0.030989010 0.04648352 0.01631001   9757     1
# WP536   WP536                 Calcium regulation in cardiac cells       1/4 152/8421 0.070281209 0.07207918 0.02529094   6004     1
# WP289   WP289      Myometrial relaxation and contraction pathways       1/4 156/8421 0.072079178 0.07207918 0.02529094   6004     1

#### For Down_HS_B_unique### 

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Down_subtype_B_unique), mart = ensmebl)
Down_subtype_B_unique$ensembl_gene_id<-row.names(Down_subtype_B_unique)

merged <- merge(Down_subtype_B_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ###

down_B <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

# No gene sets have size between 10 and 500 ...
# --> return NULL...
barplot(down_B, showCategory = 15)

dotplot(down_B, showCategory=15)

heatplot(down_B, showCategory = 5)

### We perform KEGG enrichment analysis. HS UP ###

down_B_kegg <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 1, qvalueCutoff = 1)

head(down_B_kegg,10)

# --> No gene can be mapped....
# --> Expected input gene ID: 1268,11245,539,7078,1349,64066

### Subtype T ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_subtype_T_unique), mart = ensmebl)
Up_subtype_T_unique$ensembl_gene_id<-row.names(Up_subtype_T_unique)

merged <- merge(Up_subtype_T_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ###

up_T <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(up_T, showCategory = 15)

dotplot(up_T, showCategory=15)

heatplot(up_T, showCategory = 5)

### We perform KEGG enrichment analysis. HS UP ###

up_T_kegg <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(up_T_kegg,10)

#Results 
# [1] ID          Description GeneRatio   BgRatio     pvalue      p.adjust    qvalue      geneID      Count      
# <0 rows> (or 0-length row.names)

#### For Down_HS_T_unique### 

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Down_subtype_T_unique), mart = ensmebl)
Down_subtype_T_unique$ensembl_gene_id<-row.names(Down_subtype_T_unique)

merged <- merge(Down_subtype_T_unique, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ###

down_T <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(down_T, showCategory = 15)

dotplot(down_T, showCategory=15)

heatplot(down_T, showCategory = 5)

### We perform KEGG enrichment analysis. HS UP ###

down_T_kegg <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 1, qvalueCutoff = 1)

head(down_T_kegg,10)

# --> No gene can be mapped....
# --> Expected input gene ID: 1268,11245,539,7078,1349,64066
# --> return NULL...



############ DEG adults vs pediatric of course just tumor ###### 

# clusterino_pam2 contains in Cell type columns the info on the subtypes
info_age<-clusterino_pam2

edge_c_age <- DGEList(counts = tumor_adjusted1, group=info_age$type, samples=info_age, genes=tumor_adjusted1)
edge_n_age <- calcNormFactors(edge_c_age,method = 'TMM')
# We create the cpm table
cpm_table_age <-as.data.frame(round(cpm(edge_n_age),2)) # the library size is scaled by the normalization factor

# Here we define the experimental design matrix, we build a model with no intercept also we have two varaibles, one for each condition 
# 1 for control and 2 for tumor 
design_age<- model.matrix(~0+group, data = edge_n_age$samples)
colnames(design_age) <- levels(edge_n_age$samples$group)
rownames(design_age) <- edge_n_age$samples$sample

# Calculate dispersion and fit the result with edgeR (necessary for differential expression analysis)
edge_d_age <- estimateDisp(edge_n_age,design_age)

# Fit the data we model the data using a negative binomial distribution
edge_f_age<-glmQLFit(edge_d_age, design_age)

# Definition of the contrast (conditions to be compared)
contro_age <- makeContrasts("pediatric-adult", levels=design_age)

# Fit the model with generalized linear models
edge_t_age<- glmQLFTest(edge_f_age,contrast=contro_age)

# edge_t contains the results of the DE analysis
# -> we can extract the data using the function topTags -> extract the top20, using a cut off and sorting by fold-change
# -> we get the top 20 DE genes
#  We sort for the fold change

# the n are the number of row. this case 21376

DEGs_age <- as.data.frame(topTags(edge_t_age,n=21376,p.value = 0.01,sort.by = "logFC"))

# We add a new column to the DEGs dataframe called class.
# Used to express the values of the fold change of the transcripts.
# The selection is based on the log fold change ratio (>1.5 for up-regulated genes and < (-1.5) for down-regulated genes)
# and a log CPM (>1 for both cases). From the contingency table of our DEGs we can see that the up regulated genes
# correspond to the 3.7% of the total and the down regulated are the 16% of the total.

DEGs_age$class <- '='
DEGs_age$class[which(DEGs_age$logCPM > 1 & DEGs_age$logFC > 1.5)] = '+'
DEGs_age$class[which(DEGs_age$logCPM > 1 & DEGs_age$logFC < (-1.5))] = '-'
DEGs_age <- DEGs_age[order(DEGs_age$logFC, decreasing = T),] # we order based on the fold change

table(DEGs_age$class)
# after adjustment of age, so over 18 adults under and also 18 pediatric  - 94   + 69   = 3754 -> more DEGS?!
     
# Let`s check how many human specific genes we have in the up regulated and down regulated genes in pediatric cancer
DEGs_age_HS <- DEGs_age %>% dplyr::filter(rownames(DEGs_age) %in% Human_genes$`Ensembl ID`)
Up_HS_age <- DEGs_age[DEGs_age$class=='+',] %>% dplyr::filter(rownames(DEGs_age[DEGs_age$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HS_age<- DEGs_age[DEGs_age$class=='-',] %>% dplyr::filter(rownames(DEGs_age[DEGs_age$class=='-',]) %in% Human_genes$`Ensembl ID`) 

table(DEGs_age_HS$class)
# after change  - 9   + 6  = 115
     

#### UMAP DEGs pediatric-adult only HS  #### 

Totalage_HS_unique <- rownames(DEGs_age_HS[DEGs_age_HS$class != '=',])

cpm_table_age_log <- as.data.frame(round(log10(cpm(edge_n_age)+1),2))
Age_cpm_log_only_HS <- cpm_table_age_log[which(rownames(cpm_table_age_log) %in% Totalage_HS_unique),]
Age_cpm_log_only_HS_filtered <- Age_cpm_log_only_HS[,which(apply(Age_cpm_log_only_HS, 2, var) != 0)]
Age_cpm_log_only_HS_filtered <- Age_cpm_log_only_HS_filtered[which(apply(Age_cpm_log_only_HS_filtered, 1, var) != 0),]

# Perform UMAP dimensionality reduction
umap_result6 <- umap(t(Age_cpm_log_only_HS_filtered),
                     n_neighbors = 5,         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)



umap_result6 <- umap(t(Age_cpm_log_only_HS_filtered),
                     n_neighbors = 15,         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)

umap_result6 <- umap(t(Age_cpm_log_only_HS_filtered),
                     n_neighbors = sqrt(dim(t(Age_cpm_log_only_HS_filtered))[1]),         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)

clusterino_umap_age <- clusterino_pam2
umap_df6 <- data.frame(umap_result6$layout,
                       Cell_type = clusterino_umap_age$type, Age =clusterino_umap_age$age)

colnames(umap_df6) <- c("umap_1","umap_2","umap_3", "Age",'Years')
rownames(umap_df6) <-  rownames(t(Age_cpm_log_only_HS_filtered))
# Print UMAP
print(umap_result6)

#umap_df <- cbind(umap_df,disease_state)
# Create a 3D scatter plot using Plotly
fig4U <- plot_ly(umap_df6, 
                 x = ~umap_1, y = ~umap_2, z = ~umap_3,
                 color = umap_df6$Age,
                 colors = c("#003f5c","#bc5090","#ffa600"),   
                 mode = 'markers',
                 text = paste0('Age: ',umap_df6$Years),
                 size=10) %>% layout(title = 'Pediatric-adults stratification HS, metric euclidian, neighbors = square')
# Display the 3D scatter plot
# fig4U <- fig4U %>% add_trace(x = ~umap_1, y = ~umap_2, z = ~umap_3, text = umap_df6$Years, showlegends = TRUE)

#### UMAP DEGs pediatric-adult except HS  #### 

Totalage_whHS_unique <- DEGs_age[DEGs_age$class != '=',]
Totalage_whHS_unique <- Totalage_whHS_unique[Totalage_whHS_unique %nin% Human_genes$`Ensembl ID`,]

Age_cpm_log_wh_HS <- cpm_table_age_log[which(rownames(cpm_table_age_log) %in% rownames(Totalage_whHS_unique)),]
Age_cpm_log_wh_HS_filtered <- Age_cpm_log_wh_HS[,which(apply(Age_cpm_log_wh_HS, 2, var) != 0)]
Age_cpm_log_wh_HS_filtered <- Age_cpm_log_wh_HS_filtered[which(apply(Age_cpm_log_wh_HS_filtered, 1, var) != 0),]



# Perform UMAP dimensionality reduction
umap_result7 <- umap(t(Age_cpm_log_wh_HS_filtered),
                     n_neighbors = 5,         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)



umap_result7 <- umap(t(Age_cpm_log_wh_HS_filtered),
                     n_neighbors = 15,         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)

umap_result7 <- umap(t(Age_cpm_log_wh_HS_filtered),
                     n_neighbors = sqrt(dim(t(Age_cpm_log_wh_HS_filtered))[1]),         #or the square root of the rows 
                     min_dist = 0.1,          #
                     metric = "euclidean",    #you can change it
                     n_components = 3
)

# clusterino_umap_age <- clusterino_pam2$type
umap_df7 <- data.frame(umap_result7$layout,
                       Cell_type = clusterino_umap_age)

colnames(umap_df7) <- c("umap_1","umap_2","umap_3", "Age")
rownames(umap_df7) <-  rownames(t(Age_cpm_log_wh_HS_filtered))
# Print UMAP
print(umap_result7)

#umap_df <- cbind(umap_df,disease_state)
# Create a 2D scatter plot using Plotly
fig5U <- plot_ly(umap_df7, 
                 x = ~umap_1, y = ~umap_2, z = ~umap_3,
                 color = umap_df7$Age,
                 colors =c("#003f5c","#bc5090","#ffa600"),   
                 mode = 'markers',
                 size=10) %>%  layout(title = 'Pediatric-adults stratification without HS, metric euclidian, neighbors = square')

# Display the 2D scatter plot
fig5U






##### GSEA HS Pediatric-Adult ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_HS_age), mart = ensmebl)
Up_HS_age$ensembl_gene_id<-row.names(Up_HS_age)

merged <- merge(Up_HS_age, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) HS UP ####

up_age_hs <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(up_age_hs, showCategory = 15)

dotplot(up_age_hs, showCategory=15)

heatplot(up_age_hs, showCategory = 5)

#### We perform KEGG enrichment analysis. HS UP ####

up_age <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(up_age,10)

#Results 
# WP3664 WP3664 Regulation of Wnt B catenin signaling by small molecule compounds       1/2  17/8421 0.0040    0.046  0.008   8325     1
# WP4150 WP4150                                   Wnt signaling in kidney disease       1/2  37/8421 0.0088    0.046  0.008   8325     1
# WP3680 WP3680        Physico chemical features and toxicity associated pathways       1/2  66/8421 0.0156    0.046  0.008   8325     1
# WP2571 WP2571                                 Polycystic kidney disease pathway       1/2  85/8421 0.0201    0.046  0.008   8325     1
# WP4336 WP4336      ncRNAs involved in Wnt signaling in hepatocellular carcinoma       1/2  89/8421 0.0210    0.046  0.008   8325     1
# WP4258 WP4258           lncRNA in canonical Wnt signaling and colorectal cancer       1/2  98/8421 0.0231    0.046  0.008   8325     1
# WP399   WP399                            Wnt signaling pathway and pluripotency       1/2 102/8421 0.0241    0.046  0.008   8325     1
# WP428   WP428                                                     Wnt signaling       1/2 114/8421 0.0269    0.046  0.008   8325     1
# WP3931 WP3931                         Embryonic stem cell pluripotency pathways       1/2 117/8421 0.0276    0.046  0.008   8325     1
# WP4787 WP4787                   Osteoblast differentiation and related diseases       1/2 119/8421 0.0281    0.046  0.008   8325     1

#### Gene Ontology enrichment analysis (biological process) HS Down ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Down_HS_age), mart = ensmebl)
Down_HS_age$ensembl_gene_id<-row.names(Down_HS_age)

merged <- merge(Down_HS_age, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

down_age_hs <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(down_age_hs, showCategory = 15)

dotplot(down_age_hs, showCategory=15)

heatplot(down_age_hs, showCategory = 5)

#### We perform KEGG enrichment analysis. HS Down ####

down_age_wp <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(down_age_wp, 10)

# WP22     WP22                             IL 9 signaling pathway       1/2  17/8421  0.004    0.016     NA   3581     1
# WP4917 WP4917                          Proximal tubule transport       1/2  57/8421  0.013    0.025     NA    486     1
# WP4538 WP4538 Regulatory circuits of the STAT3 signaling pathway       1/2  78/8421  0.018    0.025     NA   3581     1
# WP536   WP536                Calcium regulation in cardiac cells       1/2 152/8421  0.036    0.036     NA    486     1

##### GSEA HS Pediatric-Adult non HS ####

Up_age <- DEGs_age[DEGs_age$class=='+',]
Down_age<- DEGs_age[DEGs_age$class=='-',]

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_age), mart = ensmebl)
Up_age$ensembl_gene_id<-row.names(Up_age)

merged <- merge(Up_age, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

#### Gene Ontology enrichment analysis (biological process) UP ####

Up_age_BP <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(Up_age_BP, showCategory = 15)

dotplot(Up_age_BP, showCategory=15)

heatplot(Up_age_BP, showCategory = 5)

#### We perform KEGG enrichment analysis. UP ####

up_age <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 1, qvalueCutoff = 0.1 )

head(up_age,10)
#Why? empty

#### Gene Ontology enrichment analysis (biological process) Down ####

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Down_age), mart = ensmebl)
Down_age$ensembl_gene_id<-row.names(Down_age)

merged <- merge(Down_age, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

down_age_BP<- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(down_age_BP, showCategory = 15)

dotplot(down_age_BP, showCategory=15)

heatplot(down_age_BP, showCategory = 5)

#### We perform KEGG enrichment analysis. Down ####

down_age_wp <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(down_age_wp, 10)

#RESULTS:
# WP2877 WP2877 Vitamin D receptor pathway      5/35 187/8421 0.00097    0.068  0.060 1565/7078/9173/6280/820     5
# WP1533 WP1533     Vitamin B12 metabolism      3/35  51/8421 0.00119    0.068  0.060          3039/3043/6352     3
# WP236   WP236               Adipogenesis      4/35 131/8421 0.00202    0.076  0.067    56729/5618/2624/2662     4





############ Validation ####

# VALIDATION OF ONEGENE - STRING
# install.packages('rbioapi')
library(rbioapi)

pre_expansion <- read.csv('Expansion_files/ZNF850.csv')
colnames(pre_expansion)<-c('index', 'HGSymbol')
expansion<-pre_expansion[-1]
protein_mapped<-rba_string_map_ids(ids=expansion$HGSymbol, species=9606) # 9606 equal homo sapiens
# functional enrichment
enriched<-rba_string_enrichment(ids=protein_mapped$preferredName, species=9606)

# protein-protein interaction enrichment

enrich_interaction<- rba_string_enrichment_ppi(ids=protein_mapped$preferredName, species=9606)


annotation<-rba_string_annotations(ids=protein_mapped$preferredName, species=9606)

interactome_string<-rba_string_interactions_network(protein_mapped$preferredName, species = 9606)

# odd ratio and pvalue

expansion <- dplyr::mutate_all(expansion, .funs = toupper)
count<- interactome_string %>% dplyr::filter(interactome_string$preferredName_B %in% expansion$HGSymbol)

library(GeneOverlap)
overlap<-newGeneOverlap(expansion$HGSymbol, interactome_string$preferredName_B)
go.obj<-testGeneOverlap(overlap)
getContbl(go.obj)
print(go.obj)

# Analysis pathway and GO of the single graphs from python network x

graph1 <- read.csv('Graph1_5HS.csv')
graph1<- dplyr::mutate_all(graph1, .funs = toupper)

graph1_BP <- enrichGO(gene = graph1$X0, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.1, qvalueCutoff = 1)

barplot(graph1_BP, showCategory = 15)

dotplot(graph1_BP, showCategory=15)

heatplot(graph1_BP, showCategory = 5)

# kegg

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('external_gene_name'), values = graph1$X0, mart = ensmebl)


merged <- merge(graph1, convert, by.x = 'X0', by.y = 'external_gene_name')


wp_graph1 <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.05, qvalueCutoff = 0.1 )

head(wp_graph1,10)

# graph2

graph2 <- read.csv('Graph2_3HS.csv')
graph2<- dplyr::mutate_all(graph2, .funs = toupper)

graph2_BP <- enrichGO(gene = graph2$X0, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.1, qvalueCutoff = 1)

barplot(graph2_BP, showCategory = 15)

dotplot(graph2_BP, showCategory=15)

heatplot(graph2_BP, showCategory = 5)

# kegg

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('external_gene_name'), values = graph2$X0, mart = ensmebl)


merged <- merge(graph2, convert, by.x = 'X0', by.y = 'external_gene_name')


wp_graph2 <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(wp_graph2,10)

# graph3

graph3 <- read.csv('Graph3_3HS.csv')
graph3<- dplyr::mutate_all(graph3, .funs = toupper)

graph3_BP <- enrichGO(gene = graph3$X0, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(graph3_BP, showCategory = 15)

dotplot(graph3_BP, showCategory=15)

heatplot(graph3_BP, showCategory = 5)

# kegg

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('external_gene_name'), values = graph3$X0, mart = ensmebl)


merged <- merge(graph3, convert, by.x = 'X0', by.y = 'external_gene_name')


wp_graph3 <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(wp_graph3,10)


## TRy again the pipelin on teh genes:
# nbpf14','gtf2i','fam156a''tmem236', 'nrxn3', 'mrc1' grapl', 'lrrc37a', 'lrrc37a3', 'arl17a', 'stag3'

pre_expansion <- read.csv('Expansion_files/GRAPL.csv')
colnames(pre_expansion)<-c('index', 'HGSymbol')
expansion<-pre_expansion[-1]
protein_mapped<-rba_string_map_ids(ids=expansion$HGSymbol, species=9606) # 9606 equal homo sapiens
# functional enrichment
expansion <- dplyr::mutate_all(expansion, .funs = toupper)
count<- interactome_string %>% dplyr::filter(interactome_string$preferredName_B %in% expansion$HGSymbol)
overlap<-newGeneOverlap(expansion$HGSymbol, interactome_string$preferredName_B)
go.obj<-testGeneOverlap(overlap)
getContbl(go.obj)
print(go.obj)

