# load biocmanager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("NormqPCR")

#library(HTqPCR)
library(NormqPCR)
library(Biobase)
library(ggplot2)

# load data from csv and convert to matrix
data <- read.csv("path/to/pcr_data.csv", row.names = 1, header = TRUE)
expression_matrix <- as.matrix(data) 

samples <- colnames(expression_matrix)
genes <- rownames(expression_matrix)

genes <- as.character(genes)

pdata <- data.frame(Sample = samples, row.names = samples)
fdata <- data.frame(Symbols = genes, row.names = genes)

phenoData <- new("AnnotatedDataFrame", data = pdata)
featureData <- new("AnnotatedDataFrame", data = fdata)

qPCR_data <- new("qPCRBatch", 
                 exprs = expression_matrix, 
                 phenoData = phenoData, 
                 featureData = featureData,
                 experimentData = new("MIAME", name = "GeNorm analysis"))

genorm_results <- selectHKs(qPCR_data, Symbols = genes, method = "geNorm", minNrHK = 2)

stability_data <- data.frame(
  Gene = factor(names(genorm_results$meanM), levels = names(genorm_results$meanM)),
  M_value = genorm_results$meanM
)

ggplot(stability_data, aes(x = Gene, y = M_value)) +
  geom_point(size = 3, color = "red") +
  geom_path(aes(group = 1), color = "blue", size = 1) +
  theme_minimal() +
  labs(
    title = "Gene Stability Ranking (geNorm Analysis)",
    x = "Genes (Least Stable → Most Stable)",
    y = "M Value"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot top 2 stable genes
best_genes <- head(stability_data[order(stability_data$M_value), ], 2)
ggplot(best_genes, aes(x = Gene, y = M_value)) +
  geom_col(fill = "darkgreen", width = 0.5) +
  theme_minimal() +
  labs(
    title = "Top 2 Stable Housekeeping Genes",
    x = "Genes",
    y = "M Value"
  )
