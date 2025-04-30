library(DESeq2)

# Read in your count matrices
hn <- read.csv("highnitrogen_counts_matched333.csv")
ln <- read.csv("lownitrogen_counts_matched333.csv")
hn$condition <- "HN"
ln$condition <- "LN"
combined <- rbind(hn, ln)  # 666 rows (samples)
# Set row names
rownames(combined) <- combined$Samp

# Extract count matrix: all gene columns (excluding Samp and condition)
countdata <- combined[, !(colnames(combined) %in% c("Samp", "condition"))]

# Ensure integers
countdata <- round(as.matrix(countdata))
mode(countdata) <- "integer"

# Create metadata (colData)
coldata <- data.frame(
  row.names = rownames(combined),
  condition = factor(combined$condition, levels = c("HN", "LN"))
)
# Transpose so genes are rows, samples are columns
countdata_t <- t(countdata)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countdata_t,
                              colData = coldata,
                              design = ~ condition)

# Optional: Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition", "HN", "LN"))
# View top DE genes
head(res[order(res$pvalue), ])

# Save to CSV
write.csv(as.data.frame(res), "DESeq2_HN_vs_LN_results.csv")

# Plot
plotMA(res, main = "DESeq2 HN vs LN")
plotMA(res, main = "MA Plot: HN vs LN", ylim = c(-5, 5))
# Create dataframe
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Add significance
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")

# Plot
library(ggplot2)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: HN vs LN", x = "Log2 Fold Change", y = "-log10(FDR)")
vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation

plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA: HN vs LN")
library(pheatmap)

# Calculate distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Plot heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap")
# Select top 30 genes by adjusted p-value
topgenes <- head(order(res$padj), 30)

# Extract expression for these genes
mat <- assay(vsd)[topgenes, ]

# Add condition annotations
annotation <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])

# Heatmap
pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = annotation,
         show_rownames = TRUE, scale = "row",
         main = "Top 30 DE Genes (Scaled Expression)")
up <- sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)
down <- sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)

bar_df <- data.frame(Direction = c("Upregulated", "Downregulated"),
                     Count = c(up, down))

ggplot(bar_df, aes(x = Direction, y = Count, fill = Direction)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  labs(title = "Differentially Expressed Genes", y = "Number of Genes")
##############################################################################


# Convert results to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Remove NA padj values
res_df <- res_df[!is.na(res_df$padj), ]
up_genes <- res_df[res_df$padj < 0.05 & res_df$log2FoldChange > 1, ]

# Save to CSV
write.csv(up_genes, "upregulated_genes_HN_vs_LN.csv", row.names = FALSE)
down_genes <- res_df[res_df$padj < 0.05 & res_df$log2FoldChange < -1, ]

# Save to CSV
write.csv(down_genes, "downregulated_genes_HN_vs_LN.csv", row.names = FALSE)

sig_genes <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

write.csv(sig_genes, "all_significant_DE_genes_HN_vs_LN.csv", row.names = FALSE)






















