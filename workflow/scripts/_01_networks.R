# =============================================================================
# ARACNe Network Analysis: SCZ vs Control
# =============================================================================


# Goal: To build and compare the initial gene regulatory networks for the SCZ and HC groups
# 
# 1 Data Loading & Filtering
# What it does: It loads your raw count matrix (genes/TEs vs. samples) and filters it. You keep only the features (genes and TEs) that are expressed at a reasonably high level in at least half the samples.
# Why you did this: This is a crucial quality control step. Raw RNA-seq data is noisy. Many features have zero or near-zero counts. Including them would add noise, increase the computational burden massively, and potentially create false connections. You're focusing on the "active players" in the transcriptome.
#
# 2 Normalization (DESeq2 VST)
# What it does: It uses a Variance Stabilizing Transformation (VST) on your filtered counts.
# Why you did this: Raw counts are not directly comparable between samples due to differences in sequencing depth. Normalization corrects for this, putting all samples on a common scale. VST is a good choice because it also deals with another issue: for RNA-seq, genes with higher average expression tend to have higher variance. VST helps to moderate this, which is important for correlation-based methods like network inference.
# 3. Feature Selection (Based on Variance)
# What it does: It calculates the variance for every single gene and TE across all samples. It then selects all of these variable features (selected_features <- c(top_genes, top_tes)) to be the "nodes" in your network. Your variance plots (variance_distributions_...) show that TEs tend to be more variable than canonical genes.
# Why you did this: This is a key decision. A network is built on co-variation. Features that don't vary much across your samples (i.e., are always "on" or "off") provide little information for building a network of relationships. By selecting the most variable features, you are focusing on the most dynamic parts of the transcriptome, which are most likely to be involved in regulation and disease-state changes.
# 4. Building ARACNe Networks
# What it does: This is the core of the script.
# It separates your normalized data into two groups: SCZ and HC.
# For each group, it computes a "mutual information matrix" (MIM). This is like a giant correlation matrix, but it's more powerful because it can capture non-linear relationships, not just linear ones. It asks, "How much information does the expression of Gene A give me about the expression of Gene B?"
# It then runs the aracne algorithm. ARACNe's main job is to remove indirect connections. If A -> B and B -> C, then A and C will appear correlated. ARACNe uses a clever rule (the Data Processing Inequality) to identify and prune the weakest link in every triangle of connections (like the A-C link), leaving behind a "skeleton" of the most direct interactions.
# Why you did this: You are inferring the underlying gene regulatory logic separately for each group. This allows you to directly compare the "wiring diagrams" of healthy vs. diseased cells.
# 5. Your First Set of Graphs (What they mean)
# variance_distributions_...: These plots justify your feature selection. They show that TEs are a highly variable part of the transcriptome, supporting the idea that they are dynamic and potentially regulatory.
# degree_distributions_...: This plot shows the distribution of connections per node. In both SCZ and HC networks, you likely see a "scale-free" pattern: many nodes with few connections and a few nodes (hubs) with many connections. This is a hallmark of real-world biological networks. The key visual insight here is that the x-axis for the HC network probably extends much further to the right, meaning it has higher-degree hubs.
# degree_comparison_scz_vs_hc: This scatter plot directly compares the number of connections (degree) for each gene/TE between the two networks. If most points fall below the y=x line, it means that nearly every feature has fewer connections in the SCZ network than in the HC network. This is the first visual proof of your "network collapse."
# 








library(DESeq2)
library(minet)
library(igraph)

# Create output directory structure
# -----------------------------------------------------------------------------
dir.create("results/network_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("results/network_analysis/plots", showWarnings = FALSE)
dir.create("results/network_analysis/networks", showWarnings = FALSE)
dir.create("results/network_analysis/data", showWarnings = FALSE)
dir.create("results/network_analysis/metrics", showWarnings = FALSE)

save_plot <- function(filename, plot_function, width=10, height=5) {
  # Remove extension if provided
  filename <- sub("\\.pdf$|\\.png$", "", filename)
  
  # Save as PDF
  pdf(paste0(filename, ".pdf"), width=width, height=height)
  plot_function()
  dev.off()
  
  # Save as PNG
  png(paste0(filename, ".png"), width=width*100, height=height*100, res=100)
  plot_function()
  dev.off()
  
  cat("Plot saved:", basename(filename), "(PDF & PNG)\n")
}

# Load Data
# -----------------------------------------------------------------------------

counts <- readRDS("results/dataset/counts.rds")
features <- readRDS("results/dataset/features.rds")
samples <- read.csv("config/scz_samples_metadata.csv")

scz_samples <- samples$sample[samples$phenotype == "schizophrenia"]
hc_samples <- samples$sample[samples$phenotype == "control"]

cat("=============================================================================\n")
cat("DATASET OVERVIEW\n")
cat("=============================================================================\n")
cat("Total features in raw data:", nrow(counts), "\n")
cat("Total samples:", ncol(counts), "\n")
cat("SCZ samples:", length(scz_samples), "\n")
cat("HC samples:", length(hc_samples), "\n\n")

# -----------------------------------------------------------------------------
# EXPRESSION FILTERING
# -----------------------------------------------------------------------------

genes_idx <- features$gene.type == "CG"
tes_idx <- features$gene.type %in% c("LINE", "LTR")

# GENES: ≥10 counts in ≥50% samples, mean expression ≥10
gene_pass <- genes_idx & 
  rowSums(counts >= 10) >= 90 &
  rowMeans(counts) >= 10

# TEs: ≥10 counts in ≥50% samples, mean expression ≥5
te_pass <- tes_idx & 
  rowSums(counts >= 10) >= 90 &
  rowMeans(counts) >= 5

keep <- gene_pass | te_pass

# Create filtered dataset
counts_filtered <- counts[keep, ]
features_filtered <- features[keep, ]

rownames(counts_filtered) <- features_filtered$gene.name
rownames(features_filtered) <- features_filtered$gene.name

cat("=============================================================================\n")
cat("FILTERING RESULTS\n")
cat("=============================================================================\n")
cat("Genes passing filter:", sum(gene_pass), "\n")
cat("TEs passing filter:", sum(te_pass), "\n")
cat("  LINEs:", sum(features_filtered$gene.type == "LINE"), "\n")
cat("  LTRs:", sum(features_filtered$gene.type == "LTR"), "\n")
cat("Total features kept:", sum(keep), "\n")
cat("Filtering rate:", round(100 * sum(keep) / nrow(counts), 1), "%\n\n")


# NORMALIZATION (DESeq2 VST)
# -----------------------------------------------------------------------------

samples_ordered <- samples[match(colnames(counts_filtered), samples$sample), ]

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = samples_ordered,
  design = ~ phenotype
)

vsd <- vst(dds, blind = TRUE)
normalized_counts <- assay(vsd)
rownames(normalized_counts) <- rownames(counts_filtered)
cat("Normalization complete!\n\n")


# VARIANCE ANALYSIS
# -----------------------------------------------------------------------------

feature_variance <- apply(normalized_counts, 1, var)

genes_idx_filtered <- features_filtered$gene.type == "CG"
tes_idx_filtered <- features_filtered$gene.type %in% c("LINE", "LTR")

gene_vars <- feature_variance[genes_idx_filtered]
te_vars <- feature_variance[tes_idx_filtered]

cat("=============================================================================\n")
cat("VARIANCE STATISTICS (All Filtered Features)\n")
cat("=============================================================================\n")
cat("\nGenes (n =", length(gene_vars), "):\n")
print(summary(gene_vars))
cat("\nTEs (n =", length(te_vars), "):\n")
print(summary(te_vars))
cat("\n")

save_plot("results/network_analysis/plots/variance_distributions_separate", 
          function() {
            par(mfrow=c(1,2))
            hist(log10(gene_vars), breaks=50, 
                 main="Gene Variance Distribution", 
                 xlab="log10(variance)", col="lightblue")
            hist(log10(te_vars), breaks=50, 
                 main="TE Variance Distribution", 
                 xlab="log10(variance)", col="salmon")
          }, width=10, height=5)

save_plot("results/network_analysis/plots/variance_distributions_overlay", 
          function() {
            par(mfrow=c(1,1))
            
            breaks_range <- seq(min(log10(c(gene_vars, te_vars))), 
                                max(log10(c(gene_vars, te_vars))), 
                                length.out = 51)
            
            hist(log10(gene_vars), breaks=breaks_range, 
                 main="Variance Distribution: Genes vs TEs", 
                 xlab="log10(variance)", 
                 col=rgb(0.5, 0.7, 1, 0.5), 
                 border="blue",
                 xlim=range(breaks_range),
                 ylim=c(0, max(hist(log10(gene_vars), breaks=breaks_range, plot=FALSE)$counts,
                               hist(log10(te_vars), breaks=breaks_range, plot=FALSE)$counts) * 1.1))
            
            hist(log10(te_vars), breaks=breaks_range, 
                 col=rgb(1, 0.5, 0.5, 0.5),  
                 border="red",
                 add=TRUE)
            
            legend("topright", 
                   legend=c(paste0("Genes (n=", length(gene_vars), 
                                   ", median=", round(median(gene_vars), 2), ")"),
                            paste0("TEs (n=", length(te_vars), 
                                   ", median=", round(median(te_vars), 2), ")")),
                   fill=c(rgb(0.5, 0.7, 1, 0.5), rgb(1, 0.5, 0.5, 0.5)),
                   border=c("blue", "red"),
                   bty="n")
            
            abline(v=log10(median(gene_vars)), col="blue", lwd=2, lty=2)
            abline(v=log10(median(te_vars)), col="red", lwd=2, lty=2)
          }, width=10, height=6)

# Density plot 
save_plot("results/network_analysis/plots/variance_distributions_density", 
          function() {
            par(mfrow=c(1,1))
            
            gene_density <- density(log10(gene_vars))
            te_density <- density(log10(te_vars))
            
            plot(gene_density, 
                 main="Variance Distribution: Genes vs TEs (Density)",
                 xlab="log10(variance)",
                 ylab="Density",
                 col="blue", lwd=2,
                 xlim=range(c(gene_density$x, te_density$x)),
                 ylim=c(0, max(c(gene_density$y, te_density$y)) * 1.1))
            
            lines(te_density, col="red", lwd=2)
            
            polygon(gene_density, col=rgb(0.5, 0.7, 1, 0.3), border=NA)
            polygon(te_density, col=rgb(1, 0.5, 0.5, 0.3), border=NA)
            
            legend("topright", 
                   legend=c(paste0("Genes (n=", length(gene_vars), ")"),
                            paste0("  Median: ", round(median(gene_vars), 2)),
                            paste0("  Mean: ", round(mean(gene_vars), 2)),
                            "",
                            paste0("TEs (n=", length(te_vars), ")"),
                            paste0("  Median: ", round(median(te_vars), 2)),
                            paste0("  Mean: ", round(mean(te_vars), 2))),
                   col=c("blue", "blue", "blue", NA, "red", "red", "red"),
                   lwd=c(2, NA, NA, NA, 2, NA, NA),
                   bty="n",
                   cex=0.9)
            
            abline(v=log10(median(gene_vars)), col="blue", lwd=2, lty=2)
            abline(v=log10(median(te_vars)), col="red", lwd=2, lty=2)
            
            # Add grid
            grid(col="gray80")
          }, width=10, height=6)

# Statistical test
variance_test <- wilcox.test(gene_vars, te_vars)
cat("Wilcoxon rank-sum test p-value:", format(variance_test$p.value, scientific=TRUE), "\n")
cat("Fold difference (median TE / median Gene):", 
    round(median(te_vars) / median(gene_vars), 2), "x\n\n")


# FEATURE SELECTION, use filtered features
# -----------------------------------------------------------------------------

top_genes <- names(sort(gene_vars, decreasing = TRUE))
top_tes <- names(sort(te_vars, decreasing = TRUE))

selected_features <- c(top_genes, top_tes)
final_counts <- normalized_counts[selected_features, ]


cat("Genes:", length(top_genes), "\n")
cat("TEs:", length(top_tes), "\n")
cat("Total features for network:", length(selected_features), "\n")
cat("Matrix dimensions:", paste(dim(final_counts), collapse = " x "), "\n\n")

# Check TE types
te_type_breakdown <- table(features_filtered[top_tes, "gene.type"])
cat("TE type breakdown:\n")
print(te_type_breakdown)
cat("\n")

# -----------------------------------------------------------------------------
# EXPRESSION QUALITY VALIDATION
# -----------------------------------------------------------------------------
gene_expr <- rowMeans(counts_filtered[top_genes, ])
te_expr <- rowMeans(counts_filtered[top_tes, ])

gene_breadth <- rowSums(counts_filtered[top_genes, ] > 0)
te_breadth <- rowSums(counts_filtered[top_tes, ] > 0)

cat("=============================================================================\n")
cat("EXPRESSION QUALITY CHECK\n")
cat("=============================================================================\n")
cat("\nMean expression (raw counts):\n")
cat("Genes:\n")
print(summary(gene_expr))
cat("\nTEs:\n")
print(summary(te_expr))
cat("\nExpression breadth (# samples with >0 counts):\n")
cat("Genes:\n")
print(summary(gene_breadth))
cat("\nTEs:\n")
print(summary(te_breadth))
cat("\n")

# Statistical comparison
expr_test <- wilcox.test(gene_expr, te_expr)
breadth_test <- wilcox.test(gene_breadth, te_breadth)

cat("Statistical Comparisons:\n")
cat("Expression level (Wilcoxon test) p-value:", format(expr_test$p.value, scientific=TRUE), "\n")
cat("  Median gene expression:", round(median(gene_expr), 2), "\n")
cat("  Median TE expression:", round(median(te_expr), 2), "\n")
cat("  Fold difference:", round(median(gene_expr) / median(te_expr), 2), "x\n\n")

cat("Expression breadth (Wilcoxon test) p-value:", format(breadth_test$p.value, scientific=TRUE), "\n")
cat("  Median gene breadth:", round(median(gene_breadth), 2), "samples\n")
cat("  Median TE breadth:", round(median(te_breadth), 2), "samples\n\n")

# Visualize expression quality 
save_plot("results/network_analysis/plots/expression_quality_overview",
          function() {
            par(mfrow=c(2,3), mar=c(4,4,3,1))
            
            boxplot(list(Genes=gene_expr, TEs=te_expr), 
                    main="Mean Expression Comparison", 
                    ylab="Raw counts (log scale)", 
                    log="y", 
                    col=c("lightblue", "salmon"),
                    outline=FALSE)
            text(x=1.5, y=max(c(gene_expr, te_expr))*0.8,
                 labels=paste0("p = ", format(expr_test$p.value, digits=2, scientific=TRUE)),
                 cex=0.9)
            
            hist(log10(gene_expr), breaks=50, 
                 main="Gene Expression Distribution",
                 xlab="log10(mean counts)", 
                 col="lightblue",
                 border="white")
            abline(v=log10(median(gene_expr)), col="darkblue", lwd=2, lty=2)
            legend("topright", 
                   legend=paste0("Median: ", round(median(gene_expr), 1)),
                   bty="n", cex=0.9)
            
            # expression distribution
            hist(log10(te_expr), breaks=50, 
                 main="TE Expression Distribution",
                 xlab="log10(mean counts)", 
                 col="salmon",
                 border="white")
            abline(v=log10(median(te_expr)), col="darkred", lwd=2, lty=2)
            legend("topright", 
                   legend=paste0("Median: ", round(median(te_expr), 1)),
                   bty="n", cex=0.9)
            
            boxplot(list(Genes=gene_breadth, TEs=te_breadth),
                    main="Expression Breadth Comparison", 
                    ylab="# Samples expressing", 
                    col=c("lightblue", "salmon"),
                    outline=FALSE)
            text(x=1.5, y=max(c(gene_breadth, te_breadth))*0.9,
                 labels=paste0("p = ", format(breadth_test$p.value, digits=2, scientific=TRUE)),
                 cex=0.9)
            
            hist(gene_breadth, breaks=30, 
                 main="Gene Expression Breadth",
                 xlab="# Samples", 
                 col="lightblue",
                 border="white",
                 xlim=c(0, ncol(counts_filtered)))
            abline(v=median(gene_breadth), col="darkblue", lwd=2, lty=2)
            abline(v=ncol(counts_filtered)/2, col="gray40", lwd=1, lty=3)
            legend("topleft", 
                   legend=c(paste0("Median: ", round(median(gene_breadth), 1)),
                            "50% samples"),
                   col=c("darkblue", "gray40"),
                   lty=c(2, 3),
                   lwd=c(2, 1),
                   bty="n", cex=0.8)
            
            hist(te_breadth, breaks=30, 
                 main="TE Expression Breadth",
                 xlab="# Samples", 
                 col="salmon",
                 border="white",
                 xlim=c(0, ncol(counts_filtered)))
            abline(v=median(te_breadth), col="darkred", lwd=2, lty=2)
            abline(v=ncol(counts_filtered)/2, col="gray40", lwd=1, lty=3)
            legend("topleft", 
                   legend=c(paste0("Median: ", round(median(te_breadth), 1)),
                            "50% samples"),
                   col=c("darkred", "gray40"),
                   lty=c(2, 3),
                   lwd=c(2, 1),
                   bty="n", cex=0.8)
          }, width=15, height=10)

save_plot("results/network_analysis/plots/expression_quality_overlay",
          function() {
            par(mfrow=c(1,2), mar=c(4,4,3,1))
            
            # Expression overlay
            gene_expr_dens <- density(log10(gene_expr))
            te_expr_dens <- density(log10(te_expr))
            
            plot(gene_expr_dens, 
                 main="Expression Level: Genes vs TEs",
                 xlab="log10(mean counts)",
                 ylab="Density",
                 col="blue", lwd=2,
                 xlim=range(c(gene_expr_dens$x, te_expr_dens$x)),
                 ylim=c(0, max(c(gene_expr_dens$y, te_expr_dens$y)) * 1.1))
            
            lines(te_expr_dens, col="red", lwd=2)
            polygon(gene_expr_dens, col=rgb(0.5, 0.7, 1, 0.3), border=NA)
            polygon(te_expr_dens, col=rgb(1, 0.5, 0.5, 0.3), border=NA)
            
            abline(v=log10(median(gene_expr)), col="blue", lwd=2, lty=2)
            abline(v=log10(median(te_expr)), col="red", lwd=2, lty=2)
            
            legend("topright", 
                   legend=c(paste0("Genes (n=", length(gene_expr), ")"),
                            paste0("  Median: ", round(median(gene_expr), 1)),
                            "",
                            paste0("TEs (n=", length(te_expr), ")"),
                            paste0("  Median: ", round(median(te_expr), 1)),
                            "",
                            paste0("p = ", format(expr_test$p.value, digits=2, scientific=TRUE))),
                   col=c("blue", NA, NA, "red", NA, NA, NA),
                   lwd=c(2, NA, NA, 2, NA, NA, NA),
                   bty="n", cex=0.85)
            grid(col="gray80")
            
            # Breadth overlay
            gene_breadth_dens <- density(gene_breadth)
            te_breadth_dens <- density(te_breadth)
            
            plot(gene_breadth_dens, 
                 main="Expression Breadth: Genes vs TEs",
                 xlab="# Samples expressing",
                 ylab="Density",
                 col="blue", lwd=2,
                 xlim=range(c(gene_breadth_dens$x, te_breadth_dens$x)),
                 ylim=c(0, max(c(gene_breadth_dens$y, te_breadth_dens$y)) * 1.1))
            
            lines(te_breadth_dens, col="red", lwd=2)
            polygon(gene_breadth_dens, col=rgb(0.5, 0.7, 1, 0.3), border=NA)
            polygon(te_breadth_dens, col=rgb(1, 0.5, 0.5, 0.3), border=NA)
            
            abline(v=median(gene_breadth), col="blue", lwd=2, lty=2)
            abline(v=median(te_breadth), col="red", lwd=2, lty=2)
            abline(v=ncol(counts_filtered)/2, col="gray40", lwd=1, lty=3)
            
            legend("topleft", 
                   legend=c(paste0("Genes (n=", length(gene_breadth), ")"),
                            paste0("  Median: ", round(median(gene_breadth), 1)),
                            "",
                            paste0("TEs (n=", length(te_breadth), ")"),
                            paste0("  Median: ", round(median(te_breadth), 1)),
                            "",
                            paste0("p = ", format(breadth_test$p.value, digits=2, scientific=TRUE))),
                   col=c("blue", NA, NA, "red", NA, NA, NA),
                   lwd=c(2, NA, NA, 2, NA, NA, NA),
                   bty="n", cex=0.85)
            grid(col="gray80")
          }, width=14, height=6)

save_plot("results/network_analysis/plots/expression_vs_breadth",
          function() {
            par(mfrow=c(1,2), mar=c(4,4,3,1))
            
            # Genes
            plot(gene_breadth, gene_expr,
                 main="Genes: Expression vs Breadth",
                 xlab="# Samples expressing",
                 ylab="Mean expression (log scale)",
                 log="y",
                 pch=16, cex=0.5, col=rgb(0, 0, 1, 0.3))
            
            gene_cor <- cor(gene_breadth, log10(gene_expr), method="spearman")
            abline(lm(log10(gene_expr) ~ gene_breadth), col="darkblue", lwd=2)
            legend("bottomright",
                   legend=paste0("Spearman ρ = ", round(gene_cor, 3)),
                   bty="n", cex=1.1)
            grid(col="gray80")
            
            # TEs
            plot(te_breadth, te_expr,
                 main="TEs: Expression vs Breadth",
                 xlab="# Samples expressing",
                 ylab="Mean expression (log scale)",
                 log="y",
                 pch=16, cex=0.5, col=rgb(1, 0, 0, 0.3))
            
            te_cor <- cor(te_breadth, log10(te_expr), method="spearman")
            abline(lm(log10(te_expr) ~ te_breadth), col="darkred", lwd=2)
            legend("bottomright",
                   legend=paste0("Spearman ρ = ", round(te_cor, 3)),
                   bty="n", cex=1.1)
            grid(col="gray80")
          }, width=12, height=6)

if(length(top_tes) > 0) {
  te_types <- features_filtered[top_tes, "gene.type"]
  te_type_expr <- tapply(te_expr, te_types, median)
  te_type_breadth <- tapply(te_breadth, te_types, median)
  te_type_counts <- table(te_types)
  
  save_plot("results/network_analysis/plots/expression_by_te_type",
            function() {
              par(mfrow=c(2,2), mar=c(4,4,3,1))
              
              # 1. Count by type
              barplot(te_type_counts,
                      main="TE Count by Type",
                      ylab="Number of TEs",
                      col=c("coral", "lightcoral"),
                      border="white")
              
              # 2. Expression by type
              barplot(te_type_expr,
                      main="Median Expression by TE Type",
                      ylab="Median mean counts",
                      col=c("coral", "lightcoral"),
                      border="white")
              
              # 3. Breadth by type
              barplot(te_type_breadth,
                      main="Median Breadth by TE Type",
                      ylab="Median # samples",
                      col=c("coral", "lightcoral"),
                      border="white")
              abline(h=ncol(counts_filtered)/2, col="gray40", lwd=2, lty=2)
              
              # 4. Expression distribution by type
              te_expr_by_type <- split(log10(te_expr), te_types)
              boxplot(te_expr_by_type,
                      main="Expression Distribution by TE Type",
                      ylab="log10(mean counts)",
                      col=c("coral", "lightcoral"),
                      outline=FALSE)
            }, width=12, height=10)
  
  cat("TE Type-specific statistics:\n")
  for(te_type in names(te_type_counts)) {
    cat("\n", te_type, ":\n", sep="")
    cat("  Count:", te_type_counts[te_type], "\n")
    cat("  Median expression:", round(te_type_expr[te_type], 2), "\n")
    cat("  Median breadth:", round(te_type_breadth[te_type], 1), "samples\n")
  }
  cat("\n")
}


# -----------------------------------------------------------------------------
#ARACNe NETWORKS
# -----------------------------------------------------------------------------

scz_data <- t(final_counts[, scz_samples])
hc_data <- t(final_counts[, hc_samples])

cat("=============================================================================\n")
cat("BUILDING ARACNE NETWORKS\n")
cat("=============================================================================\n")
cat("SCZ data:", paste(dim(scz_data), collapse = " x "), "(samples x features)\n")
cat("HC data:", paste(dim(hc_data), collapse = " x "), "(samples x features)\n\n")

# Build SCZ network
cat("Building SCZ network...\n")
cat("  Step 1/2: Computing mutual information matrix...\n")
scz_time <- system.time({
  scz_mim <- build.mim(scz_data, estimator = "spearman")
  cat("  Step 2/2: Applying ARACNe algorithm (DPI)...\n")
  scz_network <- aracne(scz_mim, eps = 0.05)
})
cat("SCZ network completed in", round(scz_time[3]/60, 2), "minutes\n\n")

# Build HC network
cat("Building HC network...\n")
cat("  Step 1/2: Computing mutual information matrix...\n")
hc_time <- system.time({
  hc_mim <- build.mim(hc_data, estimator = "spearman")
  cat("  Step 2/2: Applying ARACNe algorithm (DPI)...\n")
  hc_network <- aracne(hc_mim, eps = 0.05)
})
cat("HC network completed in", round(hc_time[3]/60, 2), "minutes\n\n")

# CONVERT TO IGRAPH & CALCULATE METRICS
# -----------------------------------------------------------------------------

scz_graph <- graph_from_adjacency_matrix(
  scz_network, mode = "undirected", weighted = TRUE, diag = FALSE
)

hc_graph <- graph_from_adjacency_matrix(
  hc_network, mode = "undirected", weighted = TRUE, diag = FALSE
)

# Annotate nodes
feature_types <- ifelse(V(scz_graph)$name %in% top_genes, "gene", "TE")
V(scz_graph)$type <- feature_types
V(hc_graph)$type <- feature_types

# Calculate degrees
scz_degree <- degree(scz_graph)
hc_degree <- degree(hc_graph)

scz_gene_degree <- scz_degree[V(scz_graph)$type == "gene"]
scz_te_degree <- scz_degree[V(scz_graph)$type == "TE"]
hc_gene_degree <- hc_degree[V(hc_graph)$type == "gene"]
hc_te_degree <- hc_degree[V(hc_graph)$type == "TE"]

cat("=============================================================================\n")
cat("NETWORK SUMMARY\n")
cat("=============================================================================\n")
cat("\nSCZ Network:\n")
cat("  Nodes:", vcount(scz_graph), "\n")
cat("  Edges:", ecount(scz_graph), "\n")
cat("  Density:", round(edge_density(scz_graph), 4), "\n")
cat("  Average degree:", round(mean(scz_degree), 1), "\n")
cat("  Clustering coefficient:", round(transitivity(scz_graph), 3), "\n")

cat("\nHC Network:\n")
cat("  Nodes:", vcount(hc_graph), "\n")
cat("  Edges:", ecount(hc_graph), "\n")
cat("  Density:", round(edge_density(hc_graph), 4), "\n")
cat("  Average degree:", round(mean(hc_degree), 1), "\n")
cat("  Clustering coefficient:", round(transitivity(hc_graph), 3), "\n")

cat("\nDegree by Feature Type:\n")
cat("SCZ - Genes: mean =", round(mean(scz_gene_degree), 1), 
    "| TEs: mean =", round(mean(scz_te_degree), 1), "\n")
cat("HC  - Genes: mean =", round(mean(hc_gene_degree), 1), 
    "| TEs: mean =", round(mean(hc_te_degree), 1), "\n\n")

# -----------------------------------------------------------------------------
# VALIDATION, Degree vs Expression Correlation
# -----------------------------------------------------------------------------

cor_scz_expr <- cor(scz_te_degree, te_expr, method = "spearman")
cor_scz_var <- cor(scz_te_degree, te_vars, method = "spearman")
cor_hc_expr <- cor(hc_te_degree, te_expr, method = "spearman")
cor_hc_var <- cor(hc_te_degree, te_vars, method = "spearman")

cat("=============================================================================\n")
cat("VALIDATION: DEGREE vs EXPRESSION CORRELATION\n")
cat("=============================================================================\n")
cat("TE Degree Correlations (Spearman):\n")
cat("SCZ - Degree vs Mean Expression:", round(cor_scz_expr, 3), "\n")
cat("SCZ - Degree vs Variance:", round(cor_scz_var, 3), "\n")
cat("HC  - Degree vs Mean Expression:", round(cor_hc_expr, 3), "\n")
cat("HC  - Degree vs Variance:", round(cor_hc_var, 3), "\n\n")

# -----------------------------------------------------------------------------
#VISUALIZATION
# -----------------------------------------------------------------------------

save_plot("results/network_analysis/plots/degree_distributions",
          function() {
            par(mfrow=c(2,2))
            hist(scz_gene_degree, breaks=50, main="SCZ: Gene Degree Distribution",
                 xlab="Degree", col="lightblue", xlim=c(0, max(scz_degree)))
            hist(scz_te_degree, breaks=50, main="SCZ: TE Degree Distribution",
                 xlab="Degree", col="salmon", xlim=c(0, max(scz_degree)))
            hist(hc_gene_degree, breaks=50, main="HC: Gene Degree Distribution",
                 xlab="Degree", col="lightblue", xlim=c(0, max(hc_degree)))
            hist(hc_te_degree, breaks=50, main="HC: TE Degree Distribution",
                 xlab="Degree", col="salmon", xlim=c(0, max(hc_degree)))
          }, width=12, height=10)

save_plot("results/network_analysis/plots/degree_vs_expression_variance",
          function() {
            par(mfrow=c(2,3))
            plot(te_expr, scz_te_degree, 
                 main="SCZ: TE Degree vs Expression",
                 xlab="Mean counts (log)", ylab="Degree", log="x", pch=16, cex=0.5)
            abline(lm(scz_te_degree ~ log10(te_expr)), col="red", lwd=2)
            text(x=max(te_expr)/2, y=max(scz_te_degree)*0.9, 
                 labels=paste("ρ =", round(cor_scz_expr, 3)), col="red", cex=1.2)
            
            plot(te_vars, scz_te_degree,
                 main="SCZ: TE Degree vs Variance",
                 xlab="Variance", ylab="Degree", pch=16, cex=0.5)
            abline(lm(scz_te_degree ~ te_vars), col="red", lwd=2)
            text(x=max(te_vars)/2, y=max(scz_te_degree)*0.9, 
                 labels=paste("ρ =", round(cor_scz_var, 3)), col="red", cex=1.2)
            
            plot(te_breadth, scz_te_degree,
                 main="SCZ: TE Degree vs Expression Breadth",
                 xlab="# Samples", ylab="Degree", pch=16, cex=0.5)
            
            plot(te_expr, hc_te_degree, 
                 main="HC: TE Degree vs Expression",
                 xlab="Mean counts (log)", ylab="Degree", log="x", pch=16, cex=0.5)
            abline(lm(hc_te_degree ~ log10(te_expr)), col="blue", lwd=2)
            text(x=max(te_expr)/2, y=max(hc_te_degree)*0.9, 
                 labels=paste("ρ =", round(cor_hc_expr, 3)), col="blue", cex=1.2)
            
            plot(te_vars, hc_te_degree,
                 main="HC: TE Degree vs Variance",
                 xlab="Variance", ylab="Degree", pch=16, cex=0.5)
            abline(lm(hc_te_degree ~ te_vars), col="blue", lwd=2)
            text(x=max(te_vars)/2, y=max(hc_te_degree)*0.9, 
                 labels=paste("ρ =", round(cor_hc_var, 3)), col="blue", cex=1.2)
            
            plot(te_breadth, hc_te_degree,
                 main="HC: TE Degree vs Expression Breadth",
                 xlab="# Samples", ylab="Degree", pch=16, cex=0.5)
          }, width=15, height=10)

degree_change_genes <- hc_gene_degree - scz_gene_degree
degree_change_tes <- hc_te_degree - scz_te_degree

save_plot("results/network_analysis/plots/degree_comparison_scz_vs_hc",
          function() {
            par(mfrow=c(2,2))
            plot(scz_gene_degree, hc_gene_degree,
                 main="Gene Degree: SCZ vs HC",
                 xlab="SCZ Degree", ylab="HC Degree", pch=16, cex=0.5)
            abline(0, 1, col="red", lwd=2, lty=2)
            
            plot(scz_te_degree, hc_te_degree,
                 main="TE Degree: SCZ vs HC",
                 xlab="SCZ Degree", ylab="HC Degree", pch=16, cex=0.5, col="salmon")
            abline(0, 1, col="red", lwd=2, lty=2)
            
            hist(degree_change_genes, breaks=50,
                 main="Gene Degree Change (HC - SCZ)",
                 xlab="Degree Change", col="lightblue")
            abline(v=0, col="red", lwd=2, lty=2)
            
            hist(degree_change_tes, breaks=50,
                 main="TE Degree Change (HC - SCZ)",
                 xlab="Degree Change", col="salmon")
            abline(v=0, col="red", lwd=2, lty=2)
          }, width=12, height=10)


# -----------------------------------------------------------------------------
#SAVE RESULTS
# -----------------------------------------------------------------------------

#  processed data
saveRDS(normalized_counts, "results/network_analysis/data/normalized_counts.rds")
saveRDS(features_filtered, "results/network_analysis/data/features_filtered.rds")
saveRDS(final_counts, "results/network_analysis/data/final_counts.rds")
saveRDS(selected_features, "results/network_analysis/data/selected_features.rds")
saveRDS(list(top_genes = top_genes, top_tes = top_tes), 
        "results/network_analysis/data/feature_lists.rds")

# networks
saveRDS(scz_mim, "results/network_analysis/networks/scz_mim.rds")
saveRDS(hc_mim, "results/network_analysis/networks/hc_mim.rds")
saveRDS(scz_network, "results/network_analysis/networks/scz_network.rds")
saveRDS(hc_network, "results/network_analysis/networks/hc_network.rds")
saveRDS(scz_graph, "results/network_analysis/networks/scz_graph.rds")
saveRDS(hc_graph, "results/network_analysis/networks/hc_graph.rds")

#hc_mim.rds - The "All-Possible-Connections" Matrix. This is not a network yet. It's a fully connected matrix. The values are the mutual information scores. A high score between Gene A and Gene B means they are strongly co-expressed. It contains a massive number of potential connections, many of which are redundant or indirect.
#hc_network.rds - The "Filtered-by-ARACNe" Network. This is a sparse matrix where a value of 1 indicates a direct connection (an edge) exists, and 0 means it doesn't. It has far fewer connections than the MIM. This is the first file that truly represents a "network" in the sense of nodes and edges. 
#hc_graph.rds - The "Initial Usable" Network Object. This is a much more user-friendly data structure than a simple matrix. It contains your nodes, your edges, and any metadata you've added (like V(g)$type <- "gene"). This is the network that includes ALL TEs that passed your initial expression filter (intergenic, intronic, exonic, etc.).
#hc_graph_intergenic_only.rds - The "Scientifically Refined" Network. This graph was created by taking hc_graph.rds and filtering it down to keep only the nodes that are either canonical genes or TEs annotated as "INTERGENIC". All edges connected to the removed nodes are also discarded.


#metrics
degree_info <- list(
  scz_degree = scz_degree,
  hc_degree = hc_degree,
  scz_gene_degree = scz_gene_degree,
  scz_te_degree = scz_te_degree,
  hc_gene_degree = hc_gene_degree,
  hc_te_degree = hc_te_degree,
  degree_change_genes = degree_change_genes,
  degree_change_tes = degree_change_tes
)
saveRDS(degree_info, "results/network_analysis/metrics/degree_info.rds")

#summary statistics
summary_stats <- data.frame(
  network = c("SCZ", "HC"),
  nodes = c(vcount(scz_graph), vcount(hc_graph)),
  edges = c(ecount(scz_graph), ecount(hc_graph)),
  density = c(edge_density(scz_graph), edge_density(hc_graph)),
  avg_degree = c(mean(scz_degree), mean(hc_degree)),
  clustering = c(transitivity(scz_graph), transitivity(hc_graph)),
  gene_avg_degree = c(mean(scz_gene_degree), mean(hc_gene_degree)),
  te_avg_degree = c(mean(scz_te_degree), mean(hc_te_degree))
)
write.csv(summary_stats, "results/network_analysis/metrics/network_summary_stats.csv", 
          row.names=FALSE)
