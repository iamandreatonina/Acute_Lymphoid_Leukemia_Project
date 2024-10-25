#### Loading the data ####
library(data.table)

# List of filenames
filenames <- c(
  "Metrics_HS_AGE_batch.csv",
  "Metrics_HS_AGE.csv",
  "Metrics_HS_subtype_batch.csv",
  "Metrics_HS_Subtype.csv",
  "Metrics_nonHS_AGE_batch.csv",
  "Metrics_nonHS_AGE.csv",
  "Metrics_nonHS_subtype_batch.csv",
  "Metrics_nonHS_subtype.csv"
)

# Function to read and process each file
process_file <- function(file) {
  data <- fread(file)
  if ("V1" %in% names(data)) {
    data[, V1 := NULL]
  }
  data[, Model := gsub("Multiple Layer Perceptron", "MLP", Model)]
  data[, .SD, .SDcols = 1:3]
}

# Apply function to each file and assign to variables
datasets <- lapply(filenames, process_file)


#### Plotting the data ####
library(ggplot2)
library(tidyr)

# Function to create grouped bar plots
create_grouped_barplot <- function(data, title, filename) {
  melted_data <- melt(data, id.vars = "Model", variable.name = "Metric", value.name = "Value")
  
  plot <- ggplot(melted_data, aes(x = Model, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = round(Value, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5, fontface = "bold") +
    labs(title = title, x = "Model", y = "Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    ) +
    scale_fill_brewer(palette = "Set1")
  
  # Save plot as image
  ggsave(filename = filename, plot = plot, width = 8, height = 6, dpi = 300,bg = "white")
}

# List of titles
titles <- c(
  "Metrics HS AGE Batch",
  "Metrics HS AGE",
  "Metrics HS Subtype Batch",
  "Metrics HS Subtype",
  "Metrics nonHS AGE Batch",
  "Metrics nonHS AGE",
  "Metrics nonHS Subtype Batch",
  "Metrics nonHS Subtype"
)
filenames <- paste0(gsub(" ", "_", titles), ".png")

# Create and save grouped bar plots
plots <- lapply(seq_along(datasets), function(i) {
  create_grouped_barplot(datasets[[i]], titles[i], filenames[i])
})

# Print all plots
lapply(plots, print)