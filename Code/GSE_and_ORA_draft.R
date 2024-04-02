setwd("..../Datasets/Post_manipulation")
library(gprofiler2)
library(tidyverse)

#write.table(All_DEG_gene_names$name, file = "draft_gem.txt", sep = "\t", quote = F, row.names = F)
#### 1. Load and extract ensembl id ####
DEG_without_Human_specific <- read.csv("DEGs_nomHS_metrics.csv")
rownames(DEG_without_Human_specific) <- DEG_without_Human_specific$X
# Select the + #
DEG_without_Human_specific <- DEG_without_Human_specific %>%
  select(-X) %>%             # Removing column 'X'
  filter(class == "+")       # Filtering rows where 'class' is equal to '+'

Human_specific_in_DEG <- read.csv("DEGs_HS_metrics.csv")
rownames(Human_specific_in_DEG) <- Human_specific_in_DEG$X
# Select the + #
Human_specific_in_DEG <- Human_specific_in_DEG %>%
  select(-X) %>%             # Removing column 'X'
  filter(class == "+")       # Filtering rows where 'class' is equal to '+'

#### 2. Convert to gene symbol for convenience, gplofiler don't need thi step ####

HS_DEG_gene_names <- gconvert(query = Human_specific_in_DEG$ensembl_id, organism = "hsapiens", 
         target="ENSG", mthreshold = Inf, filter_na = TRUE)

nonHS_DEG_gene_names <- gconvert(query = DEG_without_Human_specific$ensembl_id, organism = "hsapiens", 
                              target="ENSG", mthreshold = Inf, filter_na = TRUE)

#### 2.1 Quick check with biomaRt for consistency ####
# Install and load biomaRt
library(biomaRt)

ensmebl <- useMart(biomart = 'ensembl',dataset = 'hsapiens_gene_ensembl')

gene_symbols <- getBM(attributes =c('ensembl_gene_id','entrezgene_id','external_gene_name'),
                 filters = c('ensembl_gene_id'), 
                 values = Human_specific_in_DEG$Human_specific_in_DEG.X, 
                 mart = ensmebl)
print(gene_symbols)
unique(gene_symbols$external_gene_name) == HS_DEG_gene_names$name

#rm(All_DEG,DEG_without_Human_specific,Human_specific_in_DEG)

#### 3. Use gost to to over-representation analysis (ORA) ####
setwd("..../Images/ORA_and_GSEA")
#Over-Representation Analysis (ORA):
#ORA determines whether genes within a predefined gene set are over-represented or enriched within a list of differentially expressed or 
#selected genes compared to what would be expected by chance.

#--------------------------------------------------------------------------------------------------------------------------------------#
                                                                                                                                       #
#GO:BP: Gene Ontology Biological Process (GO:BP)                                                                                       #
#GO:MF: Gene Ontology Molecular Function (GO:MF)                                                                                       #
#REAC: Reactome pathways are curated pathways representing biological processes.                                                       #
#WP: WikiPathways is a community-curated resource for biological pathways.                                                             #
#HPA: Human Protein Atlas (HPA) provides information on the expression and localization of human proteins in tissues and cells.        #
#CORUM: The CORUM database contains experimentally verified mammalian protein complexes.                                               #
#HP: Human Phenotype Ontology (HP) terms describe phenotypic abnormalities associated with human diseases or genetic variations.       #
                                                                                                                                       #
#--------------------------------------------------------------------------------------------------------------------------------------#

# nonHS ORA

# Perform GO enrichment analysis using gost() function
nonHS_gost <- gost(
  query = nonHS_DEG_gene_names$name,
  organism = "hsapiens",
  significant = TRUE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.01,
  correction_method = "gSCS",
  domain_scope = "annotated",
  numeric_ns = "",
  sources = c("GO:BP", "GO:MF", "REAC", "WP"),
  as_short_link = FALSE,
  exclude_iea = FALSE,
  highlight = TRUE
)

plot_gost <- gostplot(
  nonHS_gost,                             # The GOST results to plot
  capped = TRUE,                          # If TRUE, caps the number of displayed terms
  interactive = F                      # If TRUE, creates an interactive plot
)

plot_gost

highlighted_terms <- nonHS_gost$result$term_id[nonHS_gost$result$highlighted == TRUE] # not very intrsting 

pp <- publish_gostplot(plot_gost, highlight_terms = c("GO:0046629","GO:0045066","GO:0046632","GO:0002335","GO:0002685","GO:0071706"))

pp

# Display the plot
plot_gost

## HS ORA

HS_gost <- gost(query = HS_DEG_gene_names$name, 
                   organism = "hsapiens", 
                   significant = TRUE,
                   measure_underrepresentation = FALSE, 
                   evcodes = FALSE, 
                   user_threshold = 0.01,
                   correction_method = "fdr", 
                   domain_scope = "annotated",
                   numeric_ns = "", 
                   sources = c("GO:BP","GO:MF","REAC","WP"), 
                   as_short_link = FALSE,
                   exclude_iea = FALSE,
                   highlight = TRUE)

plot_gost <- gostplot(
  HS_gost,                             # The GOST results to plot
  capped = TRUE,                          # If TRUE, caps the number of displayed terms
  interactive = F                      # If TRUE, creates an interactive plot
)

plot_gost

#highlighted_terms <- HS_gost$result$term_id[nonHS_gost$result$highlighted == TRUE] # not very interesting 

pp <- publish_gostplot(plot_gost, highlight_terms = c("REAC:R-HSA-198933","WP:WP5092","GO:0001788"))

pp

rm(pp,plot_gost)

#### 4. ORA with EnrichR ####
library("enrichR")
library("ggplot2")

# List available databases from Enrichr
dbs <- listEnrichrDbs()
dbs <- dbs[order(dbs$libraryName),]
Databases <- data.frame(dbs$libraryName)

#### 4.1 Enrichment analysis for GO ####
# Define the databases for enrichment analysis
dbs_dd <- c("GO_Biological_Process_2023","GO_Molecular_Function_2023","WikiPathway_2023_Human")

# Perform enrichment analysis for DEGs_subtype_T
ORA_HS <- enrichr(genes = HS_DEG_gene_names$name, databases = dbs_dd)
BP_HS <- ORA_HS[[1]]
BP_HS <- BP_HS[order(BP_HS$Combined.Score, decreasing = T),]
plotEnrich(ORA_HS[[1]], showTerms = 15, numChar = 60, y = "Count", orderBy = "Combined.Score", title = "ORA on HS for BP")
# Perform enrichment analysis for DEGs_subtype_PT
ORA_nonHS <- enrichr(genes = nonHS_DEG_gene_names$name, databases = dbs_dd)
BP_nonHS <- ORA_nonHS[[1]]
BP_nonHS <- BP_nonHS[order(BP_nonHS$Combined.Score, decreasing = T),]
plotEnrich(ORA_nonHS[[1]], showTerms = 15, numChar = 60, y = "Count", orderBy = "Combined.Score",  title = "ORA on nonHS for BP")

#### 4.2 ORA for DrugMatrix ####
## NOT USEFULL FOR NOW, WAITING FOR SUBTYPES DEGS ##

# Define the databases for enrichment analysis
dbs_dd <- c("DrugMatrix")

# Perform enrichment analysis for DEGs_subtype_T
ORA_HS <- enrichr(genes = HS_DEG_gene_names$name, databases = dbs_dd)
DM_HS <- ORA_HS[[1]]
DM_HS <- DM_HS[order(DM_HS$Combined.Score, decreasing = T),]
plotEnrich(ORA_HS[[1]], showTerms = 15, numChar = 60, y = "Count", orderBy = "Combined.Score", title = "DrugMatrix on whole HS")
# Perform enrichment analysis for DEGs_subtype_PT
ORA_nonHS <- enrichr(genes = nonHS_DEG_gene_names$name, databases = dbs_dd)
DM_nonHS <- ORA_nonHS[[1]]
DM_nonHS <- DM_nonHS[order(DM_nonHS$Combined.Score, decreasing = T),]
plotEnrich(ORA_nonHS[[1]], showTerms = 15, numChar = 60, y = "Count", orderBy = "Combined.Score",  title = "DrugMatrix on whole nonHS")

#### Gene Set Enrichment Analysis (GSEA) with fgsea ####
library(fgsea)
#Gene Set Enrichment Analysis (GSEA):
#GSEA evaluates whether a predefined gene set shows statistically significant differences in expression between two biological states, 
#typically a treatment condition versus a control condition.





