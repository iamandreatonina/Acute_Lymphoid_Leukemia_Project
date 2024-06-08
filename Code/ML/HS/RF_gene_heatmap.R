library("heatmaply")
heatmaply( t(x))
heatmaply( t(x1),column_text_angle = 90,fontsize_col = 5) %>%
  layout(title = "PreB")

#pheatmap(heat_data[,-ncol(heat_data)], cutree_rows = 10, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",legend = T)

# Load the data frame
DEGs_HS_metrics <- read.csv("C:/Users/perso/Desktop/Acute_Lymphoid_Leukemia_Project-main/data/Datasets/Post_manipulation/DEGs_HS_metrics.csv")

# Set row names to the first column and remove the first column
rownames(DEGs_HS_metrics) <- DEGs_HS_metrics$X
DEGs_HS_metrics$X <- NULL

# Vector of gene names to subset
names <- c("ENSG00000164330", "ENSG00000222014", "ENSG00000169994", "ENSG00000176771", "ENSG00000109321")

# Subset the data frame based on the vector of gene names
subset_df <- filter(DEGs_HS_metrics, ensembl_id %in% rownames(RF_var_importance))
subset_df$ID <- subset_df$ensembl_id 
subset_df$ensembl_id <- NULL

RF_var_importance$ID <- rownames(RF_var_importance)

merged_df <- merge(subset_df, RF_var_importance, by = "ID")
merged_df <- merged_df %>%
  arrange(-importance)

merged_df <- subset(merged_df, Category != "Small")
#write.csv(merged_df, "Intrest_HS_genes_info_DEGs.csv", row.names = T)
# Print the subsetted data frame
print(subset_df)

library(readr)
Imp_by_Subtype <- read_csv("Imp_by_Subtype.csv")
Imp_by_Subtype$...1 <- NULL


# Select the first four columns (numerical)
numerical_cols <- names(Imp_by_Subtype)[1:4]

# Scale the numerical columns
Imp_by_Subtype_linera_traformed_1e4 <- Imp_by_Subtype %>%
  mutate(across(all_of(numerical_cols), ~ . * 1e4))

# Define thresholds for categorization
high_threshold <- max(Imp_by_Subtype_linera_traformed_1e4$B_importance) * 0.5
medium_threshold <- max(Imp_by_Subtype_linera_traformed_1e4$B_importance) * 0.25

# Categorize values
Imp_by_Subtype_linera_traformed_1e4$Category <- cut(
  Imp_by_Subtype_linera_traformed_1e4$B_importance,
  breaks = c(-Inf, medium_threshold, high_threshold, Inf),
  labels = c("Small", "Medium", "High")
)

Imp_by_Subtype_linera_traformed_1e4 <- Imp_by_Subtype_linera_traformed_1e4 %>%
  arrange(desc(B_importance))


# Create the plot
ggplot(head(Imp_by_Subtype_linera_traformed_1e4, 25), aes(x = B_importance, y = reorder(ID, B_importance), fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("lightcoral", "lightgreen", "skyblue")) +  # Choose colors for each category
  labs(title = "Top 25 Variable Importance Plot RF", x = "Importance", y = "Variable")+
  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 15),
        legend.position = "right",
        panel.grid.major.y = element_line(color = "grey", linetype = "dashed")) +
  guides(fill = guide_legend(reverse = T, title.position = "top"))

