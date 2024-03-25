# BiocManager::install('sva')
library(sva)
library(tidyverse)
library(dplyr)
library(readxl)
library(ggplot2)
library(stringr)

`%nin%` <- Negate(`%in%`)

setwd("../data/Datasets/Post_manipulation")

# 88 controls and 705 tumor samples 

##### Upload huamn specific genes 

Human_genes <- readxl::read_xlsx('Human-specific.xlsx')

##### Batch effect correction 

#Load datasets 
Tumor <- read.csv('Tumor_dataframe.csv',sep =',',header = T) # needs to be unziped 
Control <- read.csv('Controls_merged.csv',sep = ',',header = T)

# We found out a duplicated ensembl_gene_id due to the fact there isn't a 1 to 1 mapping from ensembl_gene_id and hugo_symbols
# so we are going to eliminate the less informative one 

duplicato <- Tumor$ensembl_gene_id[duplicated(Tumor$ensembl_gene_id)]
#sum <- Tumor %>% dplyr::filter(Tumor$ensembl_gene_id == duplicato) 
# rowSums(sum[2:641]) # the first one is the most informative so we use distinct()
Tumor <- distinct(Tumor,ensembl_gene_id,.keep_all =T )

Tumor_2 <- as.matrix(sapply(Tumor[2:706], as.numeric))
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
batch_tumor <- c(rep(1,38),rep(2,173),rep(3,321),rep(4,108),rep(5,65))
batch_control <- c(rep(1,10),rep(2,20),rep(3,40), rep(4,18))

# application of Combat-Seq and creation of adjusted dataframes 
tumor_adjusted <- as.data.frame(ComBat_seq(Tumor_2,batch = batch_tumor,group = NULL))
control_adjusted <- as.data.frame(ComBat_seq(Control_2, batch = batch_control, group = NULL))

# adding the ensembl_gene_id column 
colnames(control_adjusted)
control_adjusted <- add_column(control_adjusted,'ensembl_gene_id' =Control$ensembl_gene_id, .before = 'TU0049_CD4_HC')
colnames(tumor_adjusted)
tumor_adjusted <- add_column(tumor_adjusted, 'ensembl_gene_id' = Tumor$ensembl_gene_id, .before = 'GSM5491718_16.001')

##### Normalization with edgeR package 

# We use TMM method , which is a normalization method intra and inter-sample and we create CPM matrices 
library(edgeR)
# install.packages('DESeq2')
library(DESeq2)

# Let`s check how many human specific genes we have in our dataset
HSgenes_tumor <- tumor_adjusted %>% dplyr::filter(tumor_adjusted$ensembl_gene_id %in% Human_genes$`Ensembl ID`) 
HSgenes_control <- control_adjusted %>% dplyr::filter(control_adjusted$ensembl_gene_id %in% Human_genes$`Ensembl ID`) 

# Previously the result is that of 873 human specific genes there are present 498 in control 
# But after adding controls ans tumors samples ( 88 control total and 705 tumors total ) in control we have 478 and tumor 603 


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

# in this case there are too many samples so we are gonging to plot 20 samples instead of 640 
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

###### Differential gene expression analysis 

total_adjusted <- merge(control_adjusted,tumor_adjusted,by='ensembl_gene_id')
total_adjusted1 <- total_adjusted %>% column_to_rownames('ensembl_gene_id')

# creating a dataframe containing the info on the samples, this is needed to be able to perform the DGE, we set the conditions of the samples as healty (H) and tumoral (T)
info_sample_1<-data.frame("sample"=colnames(total_adjusted1))
rownames(info_sample_1)<-info_sample_1$sample
info_sample_2<-as.data.frame(str_split(string=info_sample_1$sample, pattern="_", simplify=T)) #? serve? 
colnames(info_sample_2)<-c("condition","replicate")
info_samples<-cbind(info_sample_1, info_sample_2[1:2])
info_samples$condition<-c(rep("H",88),rep("T",705)) # which are 88 healthy and 705 tumors 
info_samples$replicate<-c(rep(1,793))


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

####### 
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
DEGs <- as.data.frame(topTags(edge_t,n=18203,p.value = 0.01,sort.by = "logFC"))

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
#rigth one   -756    + 2693   = 5986

# AFTER ADDING NEW DATA
#  - 778   + 2107   = 7513
  
   
  


# Let`s check how many human specific genes we have in the up regulated and down regulated genes
#  We have 110 down-reg HS genes and 18 up-regulated HS genes
DEGs_Hsgenes <- DEGs %>% dplyr::filter(rownames(DEGs) %in% Human_genes$`Ensembl ID`)
Up_HSgenes <- DEGs[DEGs$class=='+',] %>% dplyr::filter(rownames(DEGs[DEGs$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HSgenes <- DEGs[DEGs$class=='-',] %>% dplyr::filter(rownames(DEGs[DEGs$class=='-',]) %in% Human_genes$`Ensembl ID`) 
table(DEGs_Hsgenes$class)
# after the filter of median more than 5 
# down-reg 19 and 103 up reg
# AFTER ADDING 15 down and 82 up 

# Display the results using a volcano plot (x-axes: log FoldChange, y-axes: inverse function of the p-value).
# We can see the most significant DEGs colored in green, which are genes that surpass a threshold set on both the p-value
# and the Fold Change.
jpeg(filename = '../images/Vulcano_plot_DEGs.jpeg')
input_df<-DEGs
xlabel<- "log2 FC control vs case"
ylabel<-"-log10 p-value"
par(fig=c(0,1, 0,1), mar=c(4,4,1,2), mgp=c(2, 0.75,0))
plot(DEGs$logFC,-log(DEGs$PValue, base=10), xlab=xlabel,ylab = ylabel, col=ifelse(DEGs$class=="=", "grey70", "olivedrab4"), pch=20, frame.plot=TRUE, cex=0.8, main="Volcano plot") %>% 
abline(v = 0, lty = 2, col="grey20")
dev.off()

######### vulcano hs
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

col <- rep('chartreuse4', 670)
col[which(info_samples$condition == 'T')] <- 'burlywood3' 
pal <- c('blue','white','red')
pal <- colorRampPalette(pal)(670)
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

############ heatmap human specific
library(plotly)

col <- rep('chartreuse4', 670)
col[which(info_samples$condition == 'T')] <- 'burlywood3' 
pal <- c('blue','white','red')
pal <- colorRampPalette(pal)(670)
DEGs_selected <- DEGs_Hsgenes %>% dplyr::filter(DEGs_Hsgenes$class != '=')
jpeg(filename = '../images/Heatmap_plot_DEGsHS.jpeg')
heatmap(as.matrix(cpm_table[which(rownames(cpm_table) %in% rownames(DEGs_selected)),]),ColSideColors = col, cexCol = 0.5,margins = c(4,4), col = pal, cexRow = 0.2)
dev.off()

#for improving the clusterization we set the cpm table as logarithmic
cpm_table_log <- as.data.frame(round(log10(cpm(edge_n_total)+1),2))
jpeg(filename = '../images/Heatmap_plot_DEGs_logHS.jpeg')
heatmap(as.matrix(cpm_table_log[which(rownames(cpm_table_log) %in% rownames(DEGs_selected)),]),ColSideColors = col, cexCol = 0.5,margins = c(4,4), col = pal, cexRow = 0.2)
dev.off()



##### PCA analysis 
Diff_expressed <- DEGs[which(DEGs$class != '='),]
PCA_cpm_log_nonHS <- cpm_table_log[which(rownames(cpm_table_log) %in% rownames(Diff_expressed)),]
PCA_cpm_log <- cpm_table_log[which(rownames(cpm_table_log) %in% rownames(DEGs_selected)),]

# # we need to have both for the columns and the row a variance different from zero (because divide for the varaince )
PCA_cpm_log_filtered<-PCA_cpm_log[,which(apply(PCA_cpm_log, 2, var) != 0)]
PCA_cpm_log_filtered<- PCA_cpm_log_filtered[which(apply(PCA_cpm_log_filtered, 1, var) != 0),]
color<- c(rep('darkgreen',30),rep('indianred',640))

PCA_cpm_log_nonHS_filtered <- PCA_cpm_log_nonHS[,which(apply(PCA_cpm_log_nonHS, 2, var) != 0)]
PCA_cpm_log_nonHS_filtered <- PCA_cpm_log_nonHS[which(apply(PCA_cpm_log_nonHS, 1, var) != 0),]

# # PCA plot
data.PC <- prcomp(t(PCA_cpm_log_filtered),scale. = T)
jpeg(filename = '../images/PCA_plot_DEGs_log_HS.jpeg')
plot(data.PC$x[,1:2],col=color,pch = 19) 
dev.off()

data.PC_nonHG <- prcomp(t(PCA_cpm_log_nonHS_filtered),scale. = T)
plot(data.PC_nonHG$x[,1:2],col=color,pch = 19) 

# # PCA plot of tumor only
data.PC.tumor <- prcomp(t(PCA_cpm_log_filtered[31:670]),scale. = T )
jpeg(filename = '../images/PCA_plot_DEGs_log_tumor.jpeg')
plot(data.PC.tumor$x[,1:2],pch = 19)
dev.off()

data.PC_nonHG_tumor <- prcomp(t(PCA_cpm_log_nonHS_filtered[31:670]),scale. = T )
plot(data.PC_nonHG_tumor$x[,1:2],pch = 19)


##### Partitioning around medoids, need to also to install cmake
# install.packages('factoextra')
# install.packages('cluster')
library(factoextra)
library(cluster)

# for control-tumor -> 2 clusters as seen from graphs under
fviz_nbclust(data.PC$x,FUNcluster = cluster::pam,k.max = 10)
fviz_nbclust(data.PC$x,FUNcluster = cluster::pam,k.max = 10,method = 'gap_stat')+ theme_classic()
fviz_nbclust(data.PC$x,FUNcluster = cluster::pam,k.max = 10, method = "wss")

# For subtypes of tumors
fviz_nbclust(data.PC.tumor$x,FUNcluster = cluster::pam,k.max = 15)
fviz_nbclust(data.PC.tumor$x,FUNcluster = cluster::pam,k.max = 15,method = 'gap_stat')+ theme_classic()
fviz_nbclust(data.PC.tumor$x,FUNcluster = cluster::pam,k.max = 15, method = "wss")

# for control-tumor -> 2 clusters as seen from graphs under no HS
fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 10)
fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 10,method = 'gap_stat')+ theme_classic()
fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 10, method = "wss")

# For subtypes of tumors no HS
fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 15)
fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 15,method = 'gap_stat')+ theme_classic()
fviz_nbclust(data.PC_nonHG_tumor$x,FUNcluster = cluster::pam,k.max = 15, method = "wss")

# PAM on control-tumor
pam1<-eclust(data.PC$x, "pam", k=9)
pam1_nonHS <- eclust(data.PC_nonHG$x,'pam',k=9)

# PAM tumors subtypes
pam2<-eclust(data.PC.tumor$x, "pam", k=8)
pam2_nonHS <- eclust(data.PC_nonHG_tumor$x,'pam',k=8)
##### hierarchical clustering

#calculate distances between observations and create a simple dendogram
dm <- dist(data.PC$x)
hc <- hclust(dm,method = 'average')
plot(hc,hang =-1)
rect.hclust(hc,k=2,border = 'red')
clust.vec.2<-cutree(hc,k=2)
fviz_cluster(list(data=data.PC$x, cluster=clust.vec.2))

# same for tumor subtypes
dm2 <- dist(data.PC.tumor$x)
hc2 <- hclust(dm2,method = 'average')
plot(hc2,hang =-1)
rect.hclust(hc2,k=2,border = 'red')
clust.vec<-cutree(hc2,k=8)
fviz_cluster(list(data=data.PC.tumor$x, cluster=clust.vec))

# clusters <- mutate(cpm_table_log[31:670],cluster =clust.vec
library(RColorBrewer)

clusterino_pam2<-as.data.frame((pam2$clustering))
components<-data.PC.tumor[["x"]]
components<-data.frame(components)
components<-cbind(components, clusterino_pam2)
components$PC2<- -components$PC2

components$`(pam2$clustering)` <- as.factor(components$`(pam2$clustering)`)

fig<-plot_ly(components, x=~PC1, y=~PC2, color=clusterino_pam2$`(pam2$clustering)`,colors=c('cadetblue1', 'red', 
'chartreuse3','blueviolet','blue4','darkgoldenrod2','darksalmon','seagreen4') ,type='scatter',mode='markers') #  %>%
# layout(legend = list(title = list(text = 'color')))

fig


fig2<-plot_ly(components, x=~PC1, y=~PC2, z=~PC3,color=components$`(pam2$clustering)`,colors=brewer.pal(n = 8, name = "RdBu") ,mode='markers',marker = list(size = 4)) #  %>%
fig2

#######
clusterino_pam2_nonHS<-as.data.frame((pam2_nonHS$clustering))
components_nonHS<-data.PC_nonHG_tumor$x
# components<-data.frame(components)
components_nonHS<-cbind(components_nonHS, clusterino_pam2_nonHS)
components_nonHS$PC2<- -components_nonHS$PC2
components_nonHS$`(pam2_nonHS$clustering)`<- as.factor(components_nonHS$`(pam2_nonHS$clustering)`)


fig_nonHS<-plot_ly(components_nonHS, x=~PC1, y=~PC2, color=clusterino_pam2_nonHS$`(pam2_nonHS$clustering)`,colors=c('cadetblue1', 'red', 
                                                                                            'chartreuse3','blueviolet','blue4','darkgoldenrod2','darksalmon','seagreen4') ,type='scatter',mode='markers') #  %>%
# layout(legend = list(title = list(text = 'color')))

fig_nonHS

fig2_nonHS<-plot_ly(components_nonHS, x=~PC1, y=~PC2, z=~PC3, color =components_nonHS$`(pam2_nonHS$clustering)`, colors=brewer.pal(n = 8, name = "RdBu") ,mode='markers',marker = list(size = 4)) #  %>%
fig2_nonHS
# layout(legend = list(title = list(text = 'color')))


clusterino_pam2$type <- 'pediatric'
clusterino_pam2$type[533:640] <- 'adult'
#ADJUSTED FOR DIFFERENT AGE CLASSIFICATION BETWEEN THE STUDIES 
clusterino_pam2$type[rownames(clusterino_pam2) %in% c('CMUTALLS4','T59','T91','T89','T87','T82','T81','T74','T59','T112','T102','SIHTALLS32','SIHTALLS25','SIHTALLS12','H301TALLS3','H301TALLS13','H301TALLS11','CMUTALLS9','CMUTALLS13','T67','T77','T103')] <- 'pediatric'
components2 <- as.data.frame(data.PC.tumor$x)
components2<-cbind(components2, clusterino_pam2)
components2$PC2 <- -components2$PC2
fig3<-plot_ly(components2, x=~PC1, y=~PC2, color=clusterino_pam2$type,colors=c('red2', 'blue4') ,type='scatter',mode='markers')
fig3

fig4<-plot_ly(components2, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2$type,colors=c('darkred', 'blue4') ,mode='markers')
fig4

clusterino_pam2_nonHS$type <- 'pediatric'
clusterino_pam2_nonHS$type[533:640] <- 'adult'
#ADJUSTED FOR DIFFERENT AGE CLASSIFICATION BETWEEN THE STUDIES 
clusterino_pam2_nonHS$type[rownames(clusterino_pam2_nonHS) %in% c('CMUTALLS4','T59','T91','T89','T87','T82','T81','T74','T59','T112','T102','SIHTALLS32','SIHTALLS25','SIHTALLS12','H301TALLS3','H301TALLS13','H301TALLS11','CMUTALLS9','CMUTALLS13','T67','T77','T103')] <- 'pediatric'
components2_nonHS <- as.data.frame(data.PC_nonHG_tumor$x)
components2_nonHS<-cbind(components2_nonHS, clusterino_pam2_nonHS)
components2_nonHS$PC2 <- -components2_nonHS$PC2
fig3_nonHS<-plot_ly(components2_nonHS, x=~PC1, y=~PC2, color=clusterino_pam2_nonHS$type,colors=c('red2', 'blue4') ,type='scatter',mode='markers')
fig3_nonHS

fig4_nonHS<-plot_ly(components2_nonHS, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2_nonHS$type,colors=c('darkred', 'blue4') ,mode='markers')
fig4_nonHS


setwd("~/Desktop/magistrale_Qcb/3master_QCB_first_semester_second_year/biological_data_mining_blanzieri/Laboratory_Biological_Data_Mining/Dataset/GSE181157")
metadata<-  readxl::read_xlsx('GSE181157_SampleMetadata.xlsx')
metadata$`DFCI ID` <- rownames(clusterino_pam2)[39:211]
  
setwd("~/Desktop/magistrale_Qcb/3master_QCB_first_semester_second_year/biological_data_mining_blanzieri/Laboratory_Biological_Data_Mining/Datasets_finals")

# metadata$`Final Risk`<- replace(metadata$`Final Risk`, metadata$`Final Risk` == 'Not Available', NA) 

clusterino_pam2$risk <- 'Not Available'
clusterino_pam2$risk[rownames(clusterino_pam2) %in% metadata$`DFCI ID`] <- metadata$`Final Risk`
componet3 <- data.PC.tumor$x
componet3 <- cbind(componet3,clusterino_pam2)
componet3$PC2 <- -componet3$PC2 

fig5<-plot_ly(componet3, x=~PC1, y=~PC2, color=clusterino_pam2$risk,colors=c('red2', 'blue4') ,type='scatter',mode='markers')
fig5

fig6<-plot_ly(componet3, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2$risk,colors=c('darkred', 'blue4') ,mode='markers')
fig6

#####
setwd("~/Desktop/magistrale_Qcb/3master_QCB_first_semester_second_year/biological_data_mining_blanzieri/Laboratory_Biological_Data_Mining/Dataset/GSE181157")
metadata_nonHS<-  readxl::read_xlsx('GSE181157_SampleMetadata.xlsx')
metadata_nonHS$`DFCI ID` <- rownames(clusterino_pam2_nonHS)[39:211]

setwd("~/Desktop/magistrale_Qcb/3master_QCB_first_semester_second_year/biological_data_mining_blanzieri/Laboratory_Biological_Data_Mining/Datasets_finals")

# metadata$`Final Risk`<- replace(metadata$`Final Risk`, metadata$`Final Risk` == 'Not Available', NA) 

clusterino_pam2_nonHS$risk <- 'Not Available'
clusterino_pam2_nonHS$risk[rownames(clusterino_pam2_nonHS) %in% metadata_nonHS$`DFCI ID`] <- metadata_nonHS$`Final Risk`
componet3_nonHS <- data.PC_nonHG_tumor$x
componet3_nonHS <- cbind(componet3_nonHS,clusterino_pam2_nonHS)
componet3_nonHS$PC2 <- -componet3_nonHS$PC2 

fig5_nonHS<-plot_ly(componet3_nonHS, x=~PC1, y=~PC2, color=clusterino_pam2_nonHS$risk,colors=c('red2', 'blue4') ,type='scatter',mode='markers')
fig5_nonHS

fig6_nonHS<-plot_ly(componet3, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2$risk,colors=c('darkred', 'blue4') ,mode='markers')
fig6_nonHS

#########
clusterino_pam2$Cell_type <- 'Unkown'
clusterino_pam2$Cell_type[533:640] <- 'T cell' #by letaruet of only T cells 
clusterino_pam2$Cell_type[rownames(clusterino_pam2) %in% metadata$`DFCI ID`] <- metadata$Diagnosis
componet4 <- data.PC.tumor$x
componet4 <- cbind(componet4,clusterino_pam2)
componet4$PC2 <- -componet4$PC2

fig7<-plot_ly(componet4, x=~PC1, y=~PC2, color=clusterino_pam2$Cell_type,colors=c('red', 'blue','green') ,type='scatter',mode='markers')
fig7

fig8<-plot_ly(componet4, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2$Cell_type,colors=c('darkred', 'blue4','green','orange') ,mode='markers')
fig8

clusterino_pam2_nonHS$Cell_type <- 'Unkown'
clusterino_pam2_nonHS$Cell_type[533:640] <- 'T cell' #by letaruet of only T cells 
clusterino_pam2_nonHS$Cell_type[rownames(clusterino_pam2_nonHS) %in% metadata$`DFCI ID`] <- metadata$Diagnosis
componet4_nonHS <- data.PC_nonHG_tumor$x
componet4_nonHS <- cbind(componet4_nonHS,clusterino_pam2_nonHS)
componet4_nonHS$PC2 <- -componet4_nonHS$PC2

fig7_nonHS<-plot_ly(componet4_nonHS, x=~PC1, y=~PC2, color=clusterino_pam2_nonHS$Cell_type,colors=c('red2', 'blue4') ,type='scatter',mode='markers')
fig7_nonHS

fig8_nonHS<-plot_ly(componet4_nonHS, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam2_nonHS$Cell_typ,colors=c('darkred', 'blue4','green','orange') ,mode='markers')
fig8_nonHS

# ############ 

#HS 
clusterino_pam1<-as.data.frame((pam1$clustering))
clusterino_pam1$C_T <- "Tumor"
clusterino_pam1$C_T[rownames(clusterino_pam1) %in% c("X817_T","X845_B","X845_T","X858_B","X858_T","X867_B","X867_T","X899_B","X899_T","X817_B","TU0049_CD4_HC","TU0049_CD8_HC",
                                                     "TU0051_CD4_HC","TU0051_CD8_HC","TU0054_CD4_HC","TU0054_CD8_HC","XT0130_CD4_HC","XT0130_CD8_HC","XT0133_CD4_HC","XT0133_CD8_HC",
                                                     "XT0108_CD4_HC","XT0108_CD8_HC","XT0115_CD4_HC","XT0115_CD8_HC","XT0127_CD4_HC","XT0127_CD8_HC","XT0131_CD4_HC","XT0131_CD8_HC",
                                                     "XT0141_CD4_HC","XT0141_CD8_HC")] <- 'Control'
clusterino_pam1$type <- 'pediatric'
clusterino_pam1$type[563:670] <- 'adult'
clusterino_pam1$type[rownames(clusterino_pam1) %in% c('CMUTALLS4','T59','T91','T89','T87','T82','T81','T74','T59','T112','T102','SIHTALLS32','SIHTALLS25','SIHTALLS12','H301TALLS3','H301TALLS13','H301TALLS11','CMUTALLS9','CMUTALLS13','T67','T77','T103')] <- 'pediatric'
clusterino_pam1$type[11:30] <- 'adult'
clusterino_pam1$risk <- 'Not Available'
clusterino_pam1$risk[rownames(clusterino_pam1) %in% metadata$`DFCI ID`] <- metadata$`Final Risk`
clusterino_pam1$Cell_type <- 'Unkown'
clusterino_pam1$Cell_type[1:30] <- 'Control'
clusterino_pam1$Cell_type[563:670] <- 'T Cell' #by letaruet of only T cells 
clusterino_pam1$Cell_type[rownames(clusterino_pam1) %in% metadata$`DFCI ID`] <- metadata$Diagnosis
component5 <- data.PC$x
component5<- cbind(component5,clusterino_pam1)
component5$PC2 <- -component5$PC2

### sistemare 
fig9<-plot_ly(component5, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam1$Cell_type,colors=brewer.pal(n = 4, name = "RdBu"),  symbol = clusterino_pam1$type, symbols = c('diamond','circle'), mode='markers',marker = list(size = 4))
fig9

#write.csv(component5,file='ML_HS.csv',row.names = T)
#### Non HS 
clusterino_pam1_nonHS<-as.data.frame((pam1_nonHS$clustering))
clusterino_pam1_nonHS$C_T <- "Tumor"
clusterino_pam1_nonHS$C_T[rownames(clusterino_pam1_nonHS) %in% c("X817_T","X845_B","X845_T","X858_B","X858_T","X867_B","X867_T","X899_B","X899_T","X817_B","TU0049_CD4_HC","TU0049_CD8_HC",
                                                                 "TU0051_CD4_HC","TU0051_CD8_HC","TU0054_CD4_HC","TU0054_CD8_HC","XT0130_CD4_HC","XT0130_CD8_HC","XT0133_CD4_HC","XT0133_CD8_HC",
                                                                 "XT0108_CD4_HC","XT0108_CD8_HC","XT0115_CD4_HC","XT0115_CD8_HC","XT0127_CD4_HC","XT0127_CD8_HC","XT0131_CD4_HC","XT0131_CD8_HC",
                                                                 "XT0141_CD4_HC","XT0141_CD8_HC")] <- 'Control'
clusterino_pam1_nonHS$type <- 'pediatric'
clusterino_pam1_nonHS$type[563:670] <- 'adult'
clusterino_pam1_nonHS$type[rownames(clusterino_pam1_nonHS) %in% c('CMUTALLS4','T59','T91','T89','T87','T82','T81','T74','T59','T112','T102','SIHTALLS32','SIHTALLS25','SIHTALLS12','H301TALLS3','H301TALLS13','H301TALLS11','CMUTALLS9','CMUTALLS13','T67','T77','T103')] <- 'pediatric'
clusterino_pam1_nonHS$type[11:30] <- 'adult'
clusterino_pam1_nonHS$risk <- 'Not Available'
clusterino_pam1_nonHS$risk[rownames(clusterino_pam1_nonHS) %in% metadata$`DFCI ID`] <- metadata$`Final Risk`
clusterino_pam1_nonHS$Cell_type <- 'Unkown'
clusterino_pam1_nonHS$Cell_type[1:30] <- 'Control'
clusterino_pam1_nonHS$Cell_type[563:670] <- 'T Cell' #by letaruet of only T cells 
clusterino_pam1_nonHS$Cell_type[rownames(clusterino_pam1_nonHS) %in% metadata$`DFCI ID`] <- metadata$Diagnosis
component6 <- data.PC_nonHG$x
component6<- cbind(component6,clusterino_pam1_nonHS)
component6$PC2 <- -component6$PC2
#write.csv(component6,file='ML_nonHS.csv',row.names = T)

##### sistemare 
fig9_nonHS<-plot_ly(component6, x=~PC1, y=~PC2,z=~PC3, color=clusterino_pam1_nonHS$Cell_type,colors=brewer.pal(n = 4, name = "RdBu"),  symbol = clusterino_pam1_nonHS$type, symbols = c('diamond','circle'), mode='markers',marker = list(size = 4))
fig9_nonHS

##############

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

##########

library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)


ensmebl <- useMart(biomart = 'ensembl',dataset = 'hsapiens_gene_ensembl')

############ non HS
convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = rownames(DEGs), mart = ensmebl)

DEGs_2 <- rownames_to_column(DEGs, var = 'ensembl_gene_id')

DEGs_merge_convert <- merge(DEGs_2, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

DEGs_merge_convert<-DEGs_merge_convert[which(!is.na(DEGs_merge_convert$entrezgene_id)),] 

DEGs_merge_convert<-DEGs_merge_convert[-which(duplicated(DEGs_merge_convert$entrezgene_id)),]

Up_DEGs_merge_convert <- DEGs_merge_convert %>% dplyr::filter(DEGs_merge_convert$class == '+')

Down_DEGs_merge_convert <- DEGs_merge_convert %>% dplyr::filter(DEGs_merge_convert$class == '-')
# UP
ego_BP_UP <- enrichGO(gene = Up_DEGs_merge_convert$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05)

barplot(ego_BP_UP, showCategory = 15)

dotplot(ego_BP_UP, showCategory=15)

heatplot(ego_BP_UP, showCategory = 2)

eWP_BP_UP <- enrichWP(gene =Up_DEGs_merge_convert$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.05, qvalueCutoff = 0.1 )

head(eWP_BP_UP,10)

#Down
ego_BP_DW <- enrichGO(gene = Down_DEGs_merge_convert$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.05)

barplot(ego_BP_DW, showCategory = 15)

dotplot(ego_BP_DW, showCategory=15)

heatplot(ego_BP_DW, showCategory = 2)

eWP_BP_DW <- enrichWP(gene =Down_DEGs_merge_convert$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.05, qvalueCutoff = 0.1 )

head(eWP_BP_DW@result[["Description"]],10)
################# HS

convert_HS <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = rownames(DEGs_Hsgenes), mart = ensmebl)

DEGs_Hsgenes_2 <- rownames_to_column(DEGs_Hsgenes, var = 'ensembl_gene_id')

DEGs_merge_convert_HS <- merge(DEGs_Hsgenes_2, convert_HS, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

DEGs_merge_convert_HS<-DEGs_merge_convert_HS[which(!is.na(DEGs_merge_convert_HS$entrezgene_id)),] 

DEGs_merge_convert_HS<-DEGs_merge_convert_HS[-which(duplicated(DEGs_merge_convert_HS$entrezgene_id)),]

Up_DEGs_merge_convert_HS<- DEGs_merge_convert_HS %>% dplyr::filter(DEGs_merge_convert_HS$class == '+')

Down_DEGs_merge_convert_HS <- DEGs_merge_convert_HS %>% dplyr::filter(DEGs_merge_convert_HS$class == '-')
# UP
ego_BP_UP_HS <- enrichGO(gene = Up_DEGs_merge_convert_HS$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.2, qvalueCutoff =0.2)

barplot(ego_BP_UP_HS)

dotplot(ego_BP_UP_HS, showCategory=15)

heatplot(ego_BP_UP_HS, showCategory = 2)




# pathway annotation will be done after the gene expansion, due to lack of informations, to low numbers of terms
# eWP_BP_UP_HS <- enrichWP(gene =Up_DEGs_merge_convert_HS$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 1, qvalueCutoff = 0.2)
# 
# head(eWP_BP_UP_HS,10)


#Down
ego_BP_DW_HS <- enrichGO(gene = Down_DEGs_merge_convert_HS$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.2, qvalueCutoff = 0.2)

barplot(ego_BP_DW_HS)

dotplot(ego_BP_DW_HS)

heatplot(ego_BP_DW_HS, showCategory = 2)


############ ############# DEGs

#### DEG subtype vs subtype without controls

# creating a dataframe containing the info on the samples, this is needed to be able to perform the DGE

# clusterino_pam2 contains in Cell type columns the info on the subtypes
info_subtypes<-clusterino_pam2
info_subtypes$sample<-rownames(info_subtypes)
info_subtypes<-info_subtypes[info_subtypes$Cell_type!="Unkown",]
tumors_subtype<-tumor_adjusted1[colnames(tumor_adjusted1) %in% rownames(info_subtypes)]


# Now we can create the DGEList object
# edge_c_total <- DGEList(counts = total_adjusted1, group=info_samples$condition, samples=info_samples, genes=total_adjusted1)
# edge_n_total <- calcNormFactors(edge_c_total,method = 'TMM')

#######
edge_c_subtypes <- DGEList(counts = tumors_subtype, group=info_subtypes$Cell_type, samples=info_subtypes, genes=tumors_subtype)
edge_n_subtypes <- calcNormFactors(edge_c_subtypes,method = 'TMM')
# We create the cpm table
cpm_table_subtypes <-as.data.frame(round(cpm(edge_n_subtypes),2)) # the library size is scaled by the normalization factor

# Here we define the experimental design matrix, we build a model with no intercept also we have two varaibles, one for each condition 
# 1 for control and 2 for tumor 
design_subtype <- model.matrix(~0+group, data = edge_n_subtypes$samples, contrast.arg=list(group='contr.sum'))
colnames(design_subtype) <- levels(edge_n_subtypes$samples$group)
rownames(design_subtype) <- edge_n_subtypes$samples$sample
#design_subtype[,'T cell']<-design_subtype[,'T cell']/2
#design_subtype[,'9837/3 - Pre-T ALL']<-design_subtype[,'9837/3 - Pre-T ALL']/2
colnames(design_subtype)<-c("PreB", "PreT", "T")


# Calculate dispersion and fit the result with edgeR (necessary for differential expression analysis)
edge_d_subtype <- estimateDisp(edge_n_subtypes,design_subtype)

# Fit the data we model the data using a negative binomial distribution
edge_f_subtype<-glmQLFit(edge_d_subtype, design_subtype)

# Definition of the contrast (conditions to be compared)
contro_subtype_B <- makeContrasts("PreB-(PreT+T)/2", levels=design_subtype)
#contro_subtype[,"PreB-PreT-T"]<-c(1, -0.5, -0.5)

# Fit the model with generalized linear models
edge_t_subtype_B <- glmQLFTest(edge_f_subtype,contrast=contro_subtype_B)

# edge_t contains the results of the DE analysis
# -> we can extract the data using the function topTags -> extract the top20, using a cut off and sorting by fold-change
# -> we get the top 20 DE genes
#  We sort for the fold change
DEGs_subtype_B <- as.data.frame(topTags(edge_t_subtype_B,n=21420,p.value = 0.01,sort.by = "logFC"))

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
# +: 462  -:432     =: 7894 

# Let`s check how many human specific genes we have in the up regulated and down regulated genes in Pre B type
#  We have 22 down-reg HS genes and 16 up-regulated HS genes
DEGs_subtype_B_HS <- DEGs_subtype_B %>% dplyr::filter(rownames(DEGs_subtype_B) %in% Human_genes$`Ensembl ID`)
Up_HS_PreB <- DEGs_subtype_B[DEGs_subtype_B$class=='+',] %>% dplyr::filter(rownames(DEGs_subtype_B[DEGs_subtype_B$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HS_PreB<- DEGs_subtype_B[DEGs_subtype_B$class=='-',] %>% dplyr::filter(rownames(DEGs_subtype_B[DEGs_subtype_B$class=='-',]) %in% Human_genes$`Ensembl ID`) 

# Subtype PreT vs all

# Definition of the contrast (conditions to be compared)
contro_subtype_PT <- makeContrasts("PreT-(PreB+T)/2", levels=design_subtype)
#contro_subtype[,"PreB-PreT-T"]<-c(1, -0.5, -0.5)

# Fit the model with generalized linear models
edge_t_subtype_PT<- glmQLFTest(edge_f_subtype,contrast=contro_subtype_PT)

# edge_t contains the results of the DE analysis
# -> we can extract the data using the function topTags -> extract the top20, using a cut off and sorting by fold-change
# -> we get the top 20 DE genes
#  We sort for the fold change
DEGs_subtype_PT <- as.data.frame(topTags(edge_t_subtype_PT,n=21420,p.value = 0.01,sort.by = "logFC"))

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
# +: 482  -:848   =: 6773 

# Let`s check how many human specific genes we have in the up regulated and down regulated genes in Pre B type
#  We have 32 down-reg HS genes and 16 up-regulated HS genes
DEGs_subtype_PT_HS <- DEGs_subtype_PT %>% dplyr::filter(rownames(DEGs_subtype_PT) %in% Human_genes$`Ensembl ID`)
Up_HS_PreT <- DEGs_subtype_PT[DEGs_subtype_PT$class=='+',] %>% dplyr::filter(rownames(DEGs_subtype_PT[DEGs_subtype_PT$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HS_PreT<- DEGs_subtype_PT[DEGs_subtype_PT$class=='-',] %>% dplyr::filter(rownames(DEGs_subtype_PT[DEGs_subtype_PT$class=='-',]) %in% Human_genes$`Ensembl ID`) 

# T subtype vs all
# Subtype PreT vs all

# Definition of the contrast (conditions to be compared)
contro_subtype_T <- makeContrasts("T-(PreB+PreT)/2", levels=design_subtype)
#contro_subtype[,"PreB-PreT-T"]<-c(1, -0.5, -0.5)

# Fit the model with generalized linear models
edge_t_subtype_T<- glmQLFTest(edge_f_subtype,contrast=contro_subtype_T)

# edge_t contains the results of the DE analysis
# -> we can extract the data using the function topTags -> extract the top20, using a cut off and sorting by fold-change
# -> we get the top 20 DE genes
#  We sort for the fold change
DEGs_subtype_T <- as.data.frame(topTags(edge_t_subtype_T,n=21420,p.value = 0.01,sort.by = "logFC"))

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
# +: 381  -:16   =: 3848 

# Let`s check how many human specific genes we have in the up regulated and down regulated genes in Pre B type
#  We have 0 down-reg HS genes and 33 up-regulated HS genes
DEGs_subtype_T_HS <- DEGs_subtype_T %>% dplyr::filter(rownames(DEGs_subtype_T) %in% Human_genes$`Ensembl ID`)
Up_HS_T <- DEGs_subtype_T[DEGs_subtype_T$class=='+',] %>% dplyr::filter(rownames(DEGs_subtype_T[DEGs_subtype_T$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HS_T<- DEGs_subtype_T[DEGs_subtype_T$class=='-',] %>% dplyr::filter(rownames(DEGs_subtype_T[DEGs_subtype_T$class=='-',]) %in% Human_genes$`Ensembl ID`) 

#### COMPARISON BETWEEN THE HS IN THE 3 SUBTYPES -> NOT sure of the code
Common_T_PT <- Up_HS_T[rownames(Up_HS_T) %in% rownames(Up_HS_PreT),]
# 3 in common!
Common_PT_PB <- Up_HS_PreT[rownames(Up_HS_PreT) %in% rownames(Up_HS_PreB),]
# 0 in common
Common_T_PB <- Up_HS_T[rownames(Up_HS_T) %in% rownames(Up_HS_PreB),]
# 4 in common
Common_T_PT_D <- Down_HS_T[rownames(Down_HS_T) %in% rownames(Down_HS_PreT),]
# 0 in common
Common_T_PB_D <- Down_HS_T[rownames(Down_HS_T) %in% rownames(Down_HS_PreB),]
# 0 in common
Common_PT_PB_D <- Down_HS_PreT[rownames(Down_HS_PreT) %in% rownames(Down_HS_PreB),]
# 6 in common -> PreB has 22 down hs while PreT 23

## FIle creation
#write.csv(DEGs_subtype_T,file = 'DEGs_subtype_T.csv',row.names = T, col.names = T)

#write.csv(DEGs_subtype_PT,file = 'DEGs_subtype_PT.csv',row.names = T, col.names = T)

#write.csv(DEGs_subtype_B,file = 'DEGs_subtype_B.csv',row.names = T, col.names = T)

############ DEG adults vs pediatric

# creating a dataframe containing the info on the samples, this is needed to be able to perform the DGE

# clusterino_pam2 contains in Cell type columns the info on the subtypes
info_age<-clusterino_pam2

# Now we can create the DGEList object
# edge_c_total <- DGEList(counts = total_adjusted1, group=info_samples$condition, samples=info_samples, genes=total_adjusted1)
# edge_n_total <- calcNormFactors(edge_c_total,method = 'TMM')

#######
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
DEGs_age <- as.data.frame(topTags(edge_t_age,n=21420,p.value = 0.01,sort.by = "logFC"))

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
# +:55 -: 108  =: 3508  
     

# Let`s check how many human specific genes we have in the up regulated and down regulated genes in pediatric cancer
#  We have 10 down-reg HS genes and 6 up-regulated HS genes
DEGs_age_HS <- DEGs_age %>% dplyr::filter(rownames(DEGs_age) %in% Human_genes$`Ensembl ID`)
Up_HS_age <- DEGs_age[DEGs_age$class=='+',] %>% dplyr::filter(rownames(DEGs_age[DEGs_age$class=='+',]) %in% Human_genes$`Ensembl ID`) 
Down_HS_age<- DEGs_age[DEGs_age$class=='-',] %>% dplyr::filter(rownames(DEGs_age[DEGs_age$class=='-',]) %in% Human_genes$`Ensembl ID`) 



convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_HS_age), mart = ensmebl)
Up_HS_age$ensembl_gene_id<-row.names(Up_HS_age)

merged <- merge(Up_HS_age, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

up_age_hs <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(up_age_hs, showCategory = 15)

dotplot(up_age_hs, showCategory=15)

heatplot(up_age_hs, showCategory = 5)

# kegg

up_age <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(up_age,10)

# down

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Down_HS_age), mart = ensmebl)
Down_HS_age$ensembl_gene_id<-row.names(Down_HS_age)

merged <- merge(Down_HS_age, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

down_age_hs <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(down_age_hs, showCategory = 15)

dotplot(down_age_hs, showCategory=15)

heatplot(down_age_hs, showCategory = 5)

# kegg

down_age_wp <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(down_age_wp, 10)

# non hs

Up_age <- DEGs_age[DEGs_age$class=='+',]
Down_age<- DEGs_age[DEGs_age$class=='-',]

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Up_age), mart = ensmebl)
Up_age$ensembl_gene_id<-row.names(Up_age)

merged <- merge(Up_age, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

Up_age_BP <- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(Up_age_BP, showCategory = 15)

dotplot(Up_age_BP, showCategory=15)

heatplot(Up_age_BP, showCategory = 5)

# kegg

up_age <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 1, qvalueCutoff = 0.1 )

head(up_age,10)

# down

convert <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),filters = c('ensembl_gene_id'), values = row.names(Down_age), mart = ensmebl)
Down_age$ensembl_gene_id<-row.names(Down_age)

merged <- merge(Down_age, convert, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')

down_age_BP<- enrichGO(gene = merged$external_gene_name, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1)

barplot(down_age_BP, showCategory = 15)

dotplot(down_age_BP, showCategory=15)

heatplot(down_age_BP, showCategory = 5)

# kegg

down_age_wp <- enrichWP(gene =merged$entrezgene_id, organism = 'Homo sapiens', pvalueCutoff = 0.1, qvalueCutoff = 0.1 )

head(down_age_wp, 10)

############

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

