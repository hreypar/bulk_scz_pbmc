
library(dplyr)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(pheatmap)


# set variables for paths to object
# parametrizing helps with reproducibility and saves time
counts_path <- "results/dataset/counts.rds"
features_path <- "results/dataset/features.rds"

# read in data
raw_counts <- readRDS(counts_path)
features <- readRDS(features_path)

# explore objects we just read in
head(raw_counts)
dim(raw_counts)

# which features (genes, TEs) did we quantify?
table(features$gene.type)

# this is an object that is helpful to know which samples belong to which phenotype
sample_info <- data.frame(
  sample = c("GBS_T1_GSM1869426", "GBS_T2_GSM1869427", "GBS_T3_GSM1869428",
             "Control_T1_GSM1869429", "Control_T2_GSM1869430", "Control_T3_GSM1869431"),
  condition = factor(c(rep("Condition", 3), rep("Control", 3)),
                     levels = c("Control", "Condition"))
)
# Set row names
rownames(sample_info) <- sample_info$sample

stopifnot(
  identical(colnames(raw_counts), rownames(sample_info))
)


# ============================================================================
# EXPLORING THE RAW MATRIX
# ============================================================================
# how many features have zero counts in all samples?
apply(raw_counts, 1, sum) %>%
  `==`(0) %>%
  table()
  
# remove the rows that are full of zeroes and nothing else
raw_counts <- raw_counts %>%
  .[apply(., 1, sum) != 0, ]


dim(raw_counts)

# Summary statistics per sample
sample_stats <- data.frame(
  Sample = colnames(raw_counts),
  Total_Counts = colSums(raw_counts),
  Detected_Genes = colSums(raw_counts > 0),
  Mean_Counts = colMeans(raw_counts),
  Median_Counts = apply(raw_counts, 2, median),
  Max_Count = apply(raw_counts, 2, max)
)

# median counts are still low (not unexpected for rna-seq data)
sample_stats

# Visualize library sizes
ggplot(sample_stats, aes(x = Sample, y = Total_Counts/1e6, fill = Sample)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Library Sizes (RAW counts)", y = "Total Counts (millions)") +
  geom_hline(yintercept = mean(sample_stats$Total_Counts/1e6), 
             linetype = "dashed", color = "red")

# Number of detected genes per sample
ggplot(sample_stats, aes(x = Sample, y = Detected_Genes, fill = Sample)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of Detected Genes per Sample (RAW counts)", y = "Genes with > 0 counts")

# count distribution

# Log2 transform for visualization (add pseudocount)
log2_counts <- log2(raw_counts + 1)

# Density plot of count distributions
log2_long <- pivot_longer(
  as.data.frame(log2_counts),
  cols = everything(),
  names_to = "Sample",
  values_to = "log2_count"
)

ggplot(log2_long, aes(x = log2_count, color = Sample)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  labs(title = "Density Plot of Log2(RAWcounts + 1)",
       x = "Log2(count + 1)",
       y = "Density") +
  xlim(0, 15)

# Boxplot of distributions
ggplot(log2_long, aes(x = Sample, y = log2_count, fill = Sample)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers for clarity
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Log2(RAWcounts + 1)",
       y = "Log2(count + 1)") +
  coord_cartesian(ylim = c(0, 15))


# Top Expressed Genes

# Calculate mean expression per gene
gene_means <- rowMeans(raw_counts)
names(gene_means) <- rownames(raw_counts)

# Top 20 highest expressed genes
top_genes <- sort(gene_means, decreasing = TRUE)[1:20]
top_genes_df <- data.frame(
  Gene = names(top_genes),
  Mean_Count = top_genes,
  Percent_Total = (top_genes / sum(gene_means)) * 100
)

# this is across all samples, so it is a technical exploration, not a biological one
print(top_genes_df)

# What percentage of reads do top genes account for?
top_n <- c(10, 50, 100, 500, 1000)
top_gene_contribution <- sapply(top_n, function(n) {
  top_n_genes <- sort(gene_means, decreasing = TRUE)[1:n]
  sum(top_n_genes) / sum(gene_means) * 100
})

contribution_df <- data.frame(
  Top_N_Genes = top_n,
  Percent_Reads = top_gene_contribution
)

ggplot(contribution_df, aes(x = Top_N_Genes, y = Percent_Reads)) +
  geom_line() +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Cumulative Read Distribution",
       x = "Top N Genes",
       y = "Percentage of Total Reads")

# ============================================================================
# Preview Filtering Effects by Feature Type
# ============================================================================
# First, make sure features are in the same order as raw_counts
# Assuming features has a gene.name column that matches rownames(raw_counts)
if ("gene.name" %in% colnames(features)) {
  # Reorder features to match raw_counts
  features_ordered <- features[match(rownames(raw_counts), features$gene.name), ]
  
  # Verify the order matches
  if (!all(rownames(raw_counts) == features_ordered$gene.name)) {
    stop("Feature names don't match between counts matrix and features dataframe")
  }
} else {
  # If features is already ordered correctly
  features_ordered <- features
}

# Test different filtering thresholds
filter_thresholds <- c(0, 1, 5, 10, 50, 100, 500, 1000)

# Calculate genes remaining for each feature type at each threshold
filter_results <- expand.grid(
  Threshold = filter_thresholds,
  Feature_Type = unique(features_ordered$gene.type),
  stringsAsFactors = FALSE
)

filter_results$Genes_Remaining <- apply(filter_results, 1, function(row) {
  threshold <- as.numeric(row["Threshold"])
  feature_type <- row["Feature_Type"]
  
  # Get indices for this feature type
  type_indices <- which(features_ordered$gene.type == feature_type)
  
  # Count how many pass the threshold
  sum(rowSums(raw_counts[type_indices, , drop = FALSE]) > threshold)
})

# Calculate total genes per type for percentage calculation
total_per_type <- table(features_ordered$gene.type)
filter_results$Total_Genes <- total_per_type[filter_results$Feature_Type]
filter_results$Percent_Remaining <- (filter_results$Genes_Remaining / filter_results$Total_Genes) * 100

# Print summary table
print("Filtering Effects by Feature Type:")
print(filter_results)

# Create the plot with separate lines for each feature type
p1 <- ggplot(filter_results, aes(x = Threshold, y = Genes_Remaining, 
                                 color = Feature_Type, group = Feature_Type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Effect of Count Filtering by Feature Type",
       x = "Minimum Total Count Threshold",
       y = "Number of Features Remaining",
       color = "Feature Type") +
  scale_x_log10(breaks = filter_thresholds) +
  scale_color_manual(values = c("CG" = "#1f77b4", "LINE" = "#ff7f0e", "LTR" = "#2ca02c"),
                     labels = c("CG" = "Canonical Genes", "LINE" = "LINE Elements", "LTR" = "LTR Elements")) +
  theme(legend.position = "bottom")

# Create percentage version
p2 <- ggplot(filter_results, aes(x = Threshold, y = Percent_Remaining, 
                                 color = Feature_Type, group = Feature_Type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Percentage of Features Remaining After Filtering",
       x = "Minimum Total Count Threshold",
       y = "Percentage of Features Remaining (%)",
       color = "Feature Type") +
  scale_x_log10(breaks = filter_thresholds) +
  scale_color_manual(values = c("CG" = "#1f77b4", "LINE" = "#ff7f0e", "LTR" = "#2ca02c"),
                     labels = c("CG" = "Canonical Genes", "LINE" = "LINE Elements", "LTR" = "LTR Elements")) +
  theme(legend.position = "bottom") +
  ylim(0, 100)

# Combine plots
library(cowplot)
combined_plot <- plot_grid(p1, p2, ncol = 1, labels = c("A", "B"))
print(combined_plot)

# Create a more detailed view at lower thresholds
filter_thresholds_detailed <- c(0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 40, 50)

filter_results_detailed <- expand.grid(
  Threshold = filter_thresholds_detailed,
  Feature_Type = unique(features_ordered$gene.type),
  stringsAsFactors = FALSE
)

filter_results_detailed$Genes_Remaining <- apply(filter_results_detailed, 1, function(row) {
  threshold <- as.numeric(row["Threshold"])
  feature_type <- row["Feature_Type"]
  type_indices <- which(features_ordered$gene.type == feature_type)
  sum(rowSums(raw_counts[type_indices, , drop = FALSE]) > threshold)
})

filter_results_detailed$Total_Genes <- total_per_type[filter_results_detailed$Feature_Type]
filter_results_detailed$Percent_Remaining <- (filter_results_detailed$Genes_Remaining / filter_results_detailed$Total_Genes) * 100

# Detailed plot for lower thresholds
p3 <- ggplot(filter_results_detailed, aes(x = Threshold, y = Percent_Remaining, 
                                          color = Feature_Type, group = Feature_Type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Detailed View: Effect of Low Count Filtering",
       x = "Minimum Total Count Threshold",
       y = "Percentage of Features Remaining (%)",
       color = "Feature Type") +
  scale_color_manual(values = c("CG" = "#1f77b4", "LINE" = "#ff7f0e", "LTR" = "#2ca02c"),
                     labels = c("CG" = "Canonical Genes", "LINE" = "LINE Elements", "LTR" = "LTR Elements")) +
  theme(legend.position = "bottom") +
  xlim(0, 50) +
  ylim(0, 100)

print(p3)

# Summary statistics by feature type
for (ft in unique(features_ordered$gene.type)) {
  cat("\n", ft, ":\n", sep = "")
  type_indices <- which(features_ordered$gene.type == ft)
  type_counts <- raw_counts[type_indices, ]
  
  cat("- Total features:", length(type_indices), "\n")
  cat("- Features with zero counts:", sum(rowSums(type_counts) == 0), 
      "(", round(sum(rowSums(type_counts) == 0) / length(type_indices) * 100, 1), "%)\n")
  cat("- Features with >10 total counts:", sum(rowSums(type_counts) > 10), 
      "(", round(sum(rowSums(type_counts) > 10) / length(type_indices) * 100, 1), "%)\n")
  cat("- Median total counts per feature:", median(rowSums(type_counts)), "\n")
  cat("- Mean total counts per feature:", round(mean(rowSums(type_counts)), 1), "\n")
}

# Create a boxplot showing expression distribution by feature type
expression_by_type <- data.frame(
  Total_Counts = log10(rowSums(raw_counts) + 1),
  Feature_Type = features_ordered$gene.type
)

p4 <- ggplot(expression_by_type, aes(x = Feature_Type, y = Total_Counts, fill = Feature_Type)) +
  geom_boxplot(outlier.alpha = 0.1) +
  theme_bw() +
  labs(title = "Expression Distribution by Feature Type",
       x = "Feature Type",
       y = "Log10(Total Counts + 1)") +
  scale_fill_manual(values = c("CG" = "#1f77b4", "LINE" = "#ff7f0e", "LTR" = "#2ca02c"),
                    labels = c("CG" = "Canonical Genes", "LINE" = "LINE Elements", "LTR" = "LTR Elements")) +
  theme(legend.position = "none")

print(p4)

#################################################################
#### ALL of those plots were BEFORE FILTERING AND NORMALIZING
################################################################

# ============================================================================
# DGE
# ============================================================================
# create deseq dataset
# what kind of object is this?
# 
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = sample_info,
  design = ~ condition
)
# what do you see when you use the command 
counts(dds)

# filtering out low counts
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# Apply variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot using DESeq2's built-in function
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create a PCA plot
#
# add the code to label the samples :)
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "bottom"
  ) +
  ggtitle("PCA of Samples") +
  scale_color_manual(values = c("Control" = "blue", "GBS" = "red"))

# # Save the plot
# ggsave("PCA_plot.pdf", width = 6, height = 5)

# sample distance heatmap 
# Calculate sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, rownames(colData(vsd)), sep = "_")
colnames(sampleDistMatrix) <- NULL


colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample Distance Matrix")


# Run the differential expression analysis
dds <- DESeq(dds)
# Get results
res <- results(dds,
               contrast = c("condition", "Condition", "Control"),
               alpha = 0.05)

# QC plots
plotMA(res, ylim = c(-5, 5))

# Volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)


# Merge the results dataframe with the features dataframe (contains feature type)
res_df <- merge(res_df, 
                features[, c("gene.name", "gene.type")], 
                by.x = "gene", 
                by.y = "gene.name", 
                all.x = TRUE)

# Define significance based on new criteria
res_df$expression <- "NS"
res_df$expression[res_df$log2FoldChange > 2 & res_df$padj < 0.05] <- "Overexpressed"
res_df$expression[res_df$log2FoldChange < -2 & res_df$padj < 0.05] <- "Underexpressed"

# Create a column for labeling LINE/LTR genes according to specific criteria
# even if they dont have adjusted p values (see the table)
#
# we take advantage of logical vectors to achieve this
res_df$to_label <- FALSE
res_df$to_label[(res_df$gene.type %in% c("LINE", "LTR")) & 
                  (res_df$log2FoldChange > 2 | res_df$log2FoldChange < -2) & 
                  res_df$pvalue < 0.01] <- TRUE

# we will label the TE features that have a logfc over 2 and 
# a p value (not adjusted) less than 0.05 because none on the right side are below 0.01
res_df$to_label[(res_df$gene.type %in% c("LINE", "LTR")) & 
                  (res_df$log2FoldChange > 2) & 
                  res_df$pvalue < 0.05] <- TRUE

# we will also label the extreme values (over abs 5 logfc) even if they failed significance
res_df$to_label[(res_df$gene.type %in% c("LINE", "LTR")) & 
                  (res_df$log2FoldChange > 5 | res_df$log2FoldChange < -5)] <- TRUE

# Create a combined labeling condition
res_df$label_gene <- res_df$expression != "NS" | res_df$to_label

# Create a new color column that includes black for labeled but not significant points
res_df$point_color <- res_df$expression
res_df$point_color[res_df$to_label & res_df$expression == "NS"] <- "Labeled_NS"

# Create the plot
#pdf("volcano_plot_gbs.pdf", width = 11, height = 7)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  # Add points with different shapes for gene types
  geom_point(aes(color = point_color, shape = gene.type), 
             alpha = 0.6, size = 1.5) +
  # Add labels for both significant genes and LINE/LTR genes meeting criteria
  geom_text_repel(
    data = subset(res_df, label_gene == TRUE),
    aes(label = gene),
    size = 3,
    max.overlaps = 50,  # adjust this number to show more/fewer labels
    box.padding = 0.5,
    force = 10,
    segment.color = "grey50"
  ) +
  scale_color_manual(values = c("Overexpressed" = "red", 
                                "Underexpressed" = "blue",
                                "NS" = "grey",
                                "Labeled_NS" = "black"),
                     breaks = c("Overexpressed", "Underexpressed")) +
  scale_shape_manual(values = c("CG" = 16, "LINE" = 17, "LTR" = 18)) +
  theme_bw() +
  labs(x = "Log2 Fold Change",
       y = "-Log10 P-value",
       title = "Volcano Plot of Differential Gene Expression: GBS vs Control",
       color = "Expression",
       shape = "Gene Type") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black")

#dev.off()


### save out results
#
# save the full matrix as an RDS object
# Create the directory path
dir_path <- "results/differential_gene_expression"

# Create directory if it doesn't exist
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

# Save the RDS file
saveRDS(res_df, file = file.path(dir_path, "differential_expression_results.rds"))

##### write out a csv matrix with significant results
# Filter for significant results (either padj < 0.05 OR pvalue < 0.05)
sig_results <- subset(res_df, padj < 0.05 | pvalue < 0.05)

# Order by adjusted p-value
sig_results <- sig_results[order(sig_results$padj), ]

# Write to CSV in the same directory
write.csv(sig_results, 
          file = file.path(dir_path, "significant_features_differential_expression.csv"), 
          row.names = FALSE, 
          quote = FALSE)


####### HEATMAP of top DEGs
# Filter for adjusted p-value < 0.05
sig_results_adj <- subset(sig_results, padj < 0.05)

# Get normalized counts from DESeq2 object
normalized_counts <- counts(dds, normalized=TRUE)

# Get counts only for significant genes
sig_genes <- sig_results_adj$gene
sig_counts <- normalized_counts[sig_genes, ]

# Log transform the counts
log_counts <- log2(sig_counts + 1)

# Get sample information from DESeq2 object
sample_info <- colData(dds)

# Create annotation dataframe for samples
sample_annotation <- data.frame(
  Condition = sample_info$condition,
  row.names = colnames(sig_counts)
)

# Create the heatmap
pheatmap(log_counts,
         scale = "row",  # Scale by row (z-score)
         show_rownames = TRUE,  # Show gene names
         show_colnames = TRUE,  # Show sample names
         annotation_col = sample_annotation,  # Add sample annotations
         clustering_distance_rows = "correlation",  # Use correlation for clustering
         clustering_distance_cols = "correlation",
         main = "Expression Heatmap of DE Genes (padj < 0.05)",
         fontsize_row = 8)  # Adjust if too many genes

# Save the heatmap
pdf(file.path(dir_path, "expression_heatmap_padj05.pdf"), height=10, width=8)
pheatmap(log_counts,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = sample_annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Expression Heatmap of DE Genes (padj < 0.05)",
         fontsize_row = 8)
dev.off()


####### HEATMAP of Transposable Elements
# Filter for LINE and LTR genes from significant results
sig_results_TE <- subset(sig_results, gene.type %in% c("LINE", "LTR"))

# Get counts only for LINE and LTR genes
TE_genes <- sig_results_TE$gene
TE_counts <- normalized_counts[TE_genes, ]

# Log transform the counts
log_counts_TE <- log2(TE_counts + 1)

# Create annotation dataframe for genes
gene_annotation_TE <- data.frame(
  Gene_Type = sig_results_TE$gene.type,
  row.names = sig_results_TE$gene
)

# Create the heatmap
pheatmap(log_counts_TE,
         scale = "row",  # Scale by row (z-score)
         show_rownames = TRUE,  # Show gene names
         show_colnames = TRUE,  # Show sample names
         annotation_col = sample_annotation,  # Using same sample annotation from before
         annotation_row = gene_annotation_TE,  # Add gene type annotations
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Expression Heatmap of p<0.05 LINE and LTR Genes",
         fontsize_row = 8)

# Save the heatmap
pdf(file.path(dir_path, "expression_heatmap_TE.pdf"), height=10, width=8)
pheatmap(log_counts_TE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = sample_annotation,
         annotation_row = gene_annotation_TE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Expression Heatmap of p<0.05 LINE and LTR Genes",
         fontsize_row = 8)
dev.off()


