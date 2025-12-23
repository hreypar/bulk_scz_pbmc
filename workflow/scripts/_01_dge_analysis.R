library(ggplot2)
library(dplyr)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(tidyr)
library(gridExtra)
library(matrixStats)
library(ggrepel)

counts <- readRDS("results/dataset/counts.rds")

# Read metadata
features <- readRDS("results/dataset/features.rds")

samples_metadata <- read.csv("config/scz_samples_metadata.csv")

output_dir <- "results/analysis_output"
figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Helper functions
fig_path <- function(filename) file.path(figures_dir, filename)
table_path <- function(filename) file.path(tables_dir, filename)

# Set consistent theme
theme_set(theme_minimal(base_size = 12))

# Color palette
gene_type_colors <- c("CG" = "steelblue2", "LINE" = "darkorange", "LTR" = "gold3")
phenotype_colors <- c("control" = "#4DAF4A", "schizophrenia" = "#E41A1C")

# 1. DATA PREPROCESSING QC 
###########################################

samples_metadata <- samples_metadata[match(colnames(counts), samples_metadata$sample), ]
samples_metadata$phenotype <- factor(samples_metadata$phenotype, 
                                     levels = c("control", "schizophrenia"))

stopifnot(all(colnames(counts) == samples_metadata$sample))

stopifnot(all(rownames(counts) == features$gene.name))

# Separate counts by gene type (using indices to ensure correct assignment)
idx_cg <- which(features$gene.type == "CG")
idx_line <- which(features$gene.type == "LINE")
idx_ltr <- which(features$gene.type == "LTR")

counts_cg <- counts[idx_cg, ]
counts_line <- counts[idx_line, ]
counts_ltr <- counts[idx_ltr, ]

print("Data summary:")
print(paste("Canonical genes:", nrow(counts_cg)))
print(paste("LINE elements:", nrow(counts_line)))
print(paste("LTR elements:", nrow(counts_ltr)))
print(paste("Total samples:", ncol(counts)))

# EDA
###########################################

# Library size by gene type
lib_size_by_type <- data.frame(
  sample = rep(colnames(counts), 3),
  phenotype = rep(samples_metadata$phenotype, 3),
  gene_type = factor(rep(c("CG", "LINE", "LTR"), each = ncol(counts)), 
                     levels = c("CG", "LINE", "LTR")),
  lib_size = c(colSums(counts_cg), colSums(counts_line), colSums(counts_ltr))
)

p_libsize_type <- ggplot(lib_size_by_type, aes(x = phenotype, y = lib_size, fill = phenotype)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
  facet_wrap(~ gene_type, scales = "free_y", 
             labeller = labeller(gene_type = c("CG" = "Canonical Genes", 
                                               "LINE" = "LINE Elements", 
                                               "LTR" = "LTR Elements"))) +
  scale_y_log10() +
  scale_fill_manual(values = phenotype_colors) +
  labs(title = "Library Size Distribution by Gene Type",
       y = "Total Counts (log10)", x = "Phenotype") +
  theme(legend.position = "bottom")

ggsave(fig_path("library_size_by_gene_type.png"), p_libsize_type, width = 12, height = 6, dpi = 300)

# 2.2 Proportion of reads mapping to each gene type
read_proportions <- data.frame(
  sample = colnames(counts),
  phenotype = samples_metadata$phenotype,
  CG = colSums(counts_cg) / colSums(counts) * 100,
  LINE = colSums(counts_line) / colSums(counts) * 100,
  LTR = colSums(counts_ltr) / colSums(counts) * 100
)

read_prop_long <- read_proportions %>%
  pivot_longer(cols = c(CG, LINE, LTR), names_to = "gene_type", values_to = "percentage") %>%
  mutate(gene_type = factor(gene_type, levels = c("CG", "LINE", "LTR")))

p_read_prop <- ggplot(read_prop_long, aes(x = gene_type, y = percentage, fill = phenotype)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 0.8) +
  scale_fill_manual(values = phenotype_colors) +
  labs(title = "Percentage of Reads by Gene Type",
       y = "Percentage of Total Reads", x = "Gene Type") +
  theme(legend.position = "bottom")

ggsave(fig_path("read_proportion_by_gene_type.png"), p_read_prop, width = 10, height = 6, dpi = 300)

# Statistical test for proportion differences
cat("\n=== Statistical Tests for Read Proportions ===\n")
for (gtype in c("CG", "LINE", "LTR")) {
  control_prop <- read_proportions[read_proportions$phenotype == "control", gtype]
  scz_prop <- read_proportions[read_proportions$phenotype == "schizophrenia", gtype]
  test_result <- wilcox.test(control_prop, scz_prop)
  cat(sprintf("\n%s proportion - Control vs SCZ: p = %.4f\n", gtype, test_result$p.value))
  cat(sprintf("Control mean: %.2f%%, SCZ mean: %.2f%%\n", 
              mean(control_prop), mean(scz_prop)))
}

# 2.4 Gene detection rate by type
detection_by_type <- data.frame(
  gene_type = factor(c(rep("CG", nrow(counts_cg)), 
                       rep("LINE", nrow(counts_line)), 
                       rep("LTR", nrow(counts_ltr))),
                     levels = c("CG", "LINE", "LTR")),
  detection_rate = c(
    rowSums(counts_cg > 0) / ncol(counts_cg),
    rowSums(counts_line > 0) / ncol(counts_line),
    rowSums(counts_ltr > 0) / ncol(counts_ltr)
  )
)

p_detection <- ggplot(detection_by_type, aes(x = detection_rate, fill = gene_type)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  facet_wrap(~ gene_type, scales = "free_y",
             labeller = labeller(gene_type = c("CG" = "Canonical Genes", 
                                               "LINE" = "LINE Elements", 
                                               "LTR" = "LTR Elements"))) +
  scale_fill_manual(values = gene_type_colors) +
  labs(title = "Gene Detection Rate by Type",
       x = "Proportion of Samples with Expression",
       y = "Number of Genes/Elements") +
  theme(legend.position = "none")

ggsave(fig_path("detection_rate_by_gene_type.png"), p_detection, width = 12, height = 6, dpi = 300)

# FILTERING and NORMALIZATION            
###########################################

# Filter function
filter_genes <- function(count_matrix, min_count = 10, min_samples = 5) {
  rowSums(count_matrix >= min_count) >= min_samples
}

# Apply filtering with gene names
keep_cg <- filter_genes(counts_cg, min_count = 10, min_samples = 5)
keep_line <- filter_genes(counts_line, min_count = 5, min_samples = 3)
keep_ltr <- filter_genes(counts_ltr, min_count = 5, min_samples = 3)

# Get gene names to keep
genes_keep_cg <- rownames(counts_cg)[keep_cg]
genes_keep_line <- rownames(counts_line)[keep_line]
genes_keep_ltr <- rownames(counts_ltr)[keep_ltr]

# Combine all genes to keep
all_genes_keep <- c(genes_keep_cg, genes_keep_line, genes_keep_ltr)

# Filter counts and features
counts_filtered <- counts[all_genes_keep, ]
features_filtered <- features[features$gene.name %in% all_genes_keep, ]

# Reorder features to match counts
features_filtered <- features_filtered[match(rownames(counts_filtered), features_filtered$gene.name), ]

# Verify order
stopifnot(all(rownames(counts_filtered) == features_filtered$gene.name))

cat("\n=== Filtering Results ===\n")
cat(sprintf("CG: %d -> %d (%.1f%%)\n", 
            sum(features$gene.type == "CG"), length(genes_keep_cg), 
            length(genes_keep_cg)/sum(features$gene.type == "CG")*100))
cat(sprintf("LINE: %d -> %d (%.1f%%)\n", 
            sum(features$gene.type == "LINE"), length(genes_keep_line), 
            length(genes_keep_line)/sum(features$gene.type == "LINE")*100))
cat(sprintf("LTR: %d -> %d (%.1f%%)\n", 
            sum(features$gene.type == "LTR"), length(genes_keep_ltr), 
            length(genes_keep_ltr)/sum(features$gene.type == "LTR")*100))
cat(sprintf("Total: %d -> %d\n", nrow(counts), nrow(counts_filtered)))

###########################################
# DIFFERENTIAL EXPRESSION ANALYSIS     
###########################################

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = samples_metadata,
  design = ~ phenotype
)

# Add gene type information to dds
mcols(dds)$gene_type <- features_filtered$gene.type

# Run DESeq2
dds <- DESeq(dds)

# Get results
res_all <- results(dds, contrast = c("phenotype", "schizophrenia", "control"))

# Add gene information to results
res_all$gene_name <- rownames(res_all)
res_all$gene_type <- mcols(dds)$gene_type[match(rownames(res_all), rownames(dds))]

# Verify gene types
cat("\n=== Gene type distribution in results ===\n")
print(table(res_all$gene_type))

# Separate results by gene type
res_cg <- res_all[res_all$gene_type == "CG", ]
res_line <- res_all[res_all$gene_type == "LINE", ]
res_ltr <- res_all[res_all$gene_type == "LTR", ]

# Function to summarize DE results
summarize_de <- function(res, name, padj_thresh=0.01, top=5) {
  sig <- sum(res$padj < padj_thresh, na.rm = TRUE)
  sig_up <- sum(res$padj < padj_thresh & res$log2FoldChange > 0, na.rm = TRUE)
  sig_down <- sum(res$padj < padj_thresh & res$log2FoldChange < 0, na.rm = TRUE)
  
  cat(sprintf("\n%s DE Summary:\n", name))
  cat(sprintf("Total tested: %d\n", nrow(res)))
  cat(sprintf("Significant (padj < 0.01): %d (%.1f%%)\n", sig, sig/nrow(res)*100))
  cat(sprintf("  - Upregulated in SCZ: %d\n", sig_up))
  cat(sprintf("  - Downregulated in SCZ: %d\n", sig_down))
  
  # Show top 5 genes
  if (sig > 0) {
    top_genes <- head(res[order(abs(res$log2FoldChange), decreasing = T),
                          c("gene_name", "log2FoldChange", "padj")], top)
    
    cat(sprintf("\nTop %d significant genes:\n", top))
    print(top_genes)
  }
}

cat("\n=== Differential Expression Results ===")
summarize_de(res_cg, "Canonical Genes", top = 10)
summarize_de(res_line, "LINE Elements", top=10)
summarize_de(res_ltr, "LTR Elements", top=10)
summarize_de(res_all, "All Features", top=10)


# Save results tables
write.csv(as.data.frame(res_all), table_path("DE_results_all_genes.csv"), row.names = FALSE)
write.csv(as.data.frame(res_cg), table_path("DE_results_canonical_genes.csv"), row.names = FALSE)
write.csv(as.data.frame(res_line), table_path("DE_results_LINE_elements.csv"), row.names = FALSE)
write.csv(as.data.frame(res_ltr), table_path("DE_results_LTR_elements.csv"), row.names = FALSE)


# save diff expressed features with p value < 0.01
res_all_significant <- res_all[!is.na(res_all$padj) & res_all$padj < 0.01, ]

# Sort by adjusted pval
res_all_significant <- res_all_significant[order(res_all_significant$padj), ]

res_all_significant$direction <- ifelse(res_all_significant$log2FoldChange > 0, "Up", "Down")

# the number of significant results
cat(sprintf("Total significant features (padj < 0.01): %d\n", nrow(res_all_significant)))
cat(sprintf("  - Upregulated: %d\n", sum(res_all_significant$log2FoldChange > 0)))
cat(sprintf("  - Downregulated: %d\n", sum(res_all_significant$log2FoldChange < 0)))

write.csv(as.data.frame(res_all_significant), 
          table_path("DE_results_significant_only_padj001.csv"), 
          row.names = FALSE)


# ANALYZE SIGNIFICANT DE GENES BY TYPE WITH VARYING LOG2FC THRESHOLDS
# =============================================================================

analyze_significant_composition <- function(padj_thresh = 0.01) {
  
  # Get ALL significantly DE genes (padj < 0.01, no logFC filter)
  sig_all <- res_all[res_all$padj < padj_thresh & !is.na(res_all$padj), ]
  
  cat(sprintf("=== SIGNIFICANT DE COMPOSITION (padj < %g) ===\n", padj_thresh))
  cat(sprintf("Total significant genes: %d\n", nrow(sig_all)))
  
  # Base composition (no logFC filter)
  base_composition <- table(sig_all$gene_type)
  cat("\nNo logFC filter:\n")
  for(type in names(base_composition)) {
    cat(sprintf("  %s: %d (%.1f%%)\n", type, base_composition[type], 
                base_composition[type]/sum(base_composition)*100))
  }
  
  # Determine logFC thresholds from the data
  logfc_values <- abs(sig_all$log2FoldChange)
  thresholds <- c(0, 0.25, 0.5, 0.75, 1.0, 1.5)  # 0, 1.19x, 1.41x, 1.68x, 2x, 2.83x fold
  
  cat(sprintf("\nlogFC distribution in significant genes:\n"))
  cat(sprintf("  Min: %.3f, Median: %.3f, Max: %.3f\n", 
              min(logfc_values), median(logfc_values), max(logfc_values)))
  cat(sprintf("  25th percentile: %.3f, 75th percentile: %.3f\n", 
              quantile(logfc_values, 0.25), quantile(logfc_values, 0.75)))
  
  # Test different logFC thresholds
  results_by_threshold <- list()
  
  for(thresh in thresholds) {
    if(thresh == 0) {
      filtered <- sig_all
    } else {
      filtered <- sig_all[abs(sig_all$log2FoldChange) >= thresh, ]
    }
    
    if(nrow(filtered) > 0) {
      composition <- table(filtered$gene_type)
      
      cat(sprintf("\n|log2FC| >= %.2f (%.2fx fold change):\n", thresh, 2^thresh))
      cat(sprintf("  Total: %d genes (%.1f%% of significant)\n", 
                  nrow(filtered), nrow(filtered)/nrow(sig_all)*100))
      
      for(type in c("CG", "LINE", "LTR")) {
        count <- ifelse(type %in% names(composition), composition[type], 0)
        cat(sprintf("    %s: %d (%.1f%%)\n", type, count, count/nrow(filtered)*100))
      }
      
      results_by_threshold[[as.character(thresh)]] <- list(
        threshold = thresh,
        total = nrow(filtered),
        composition = composition,
        percentage_of_significant = nrow(filtered)/nrow(sig_all)*100
      )
    }
  }
  
  return(list(
    all_significant = sig_all,
    by_threshold = results_by_threshold,
    base_composition = base_composition
  ))
}

# SUMMARY TABLE
# =============================================================================
create_composition_table <- function(analysis_results) {
  
  thresholds <- names(analysis_results$by_threshold)
  
  summary_table <- data.frame(
    logFC_threshold = numeric(),
    fold_change = numeric(),
    total_genes = numeric(),
    pct_of_significant = numeric(),
    CG_count = numeric(),
    CG_percent = numeric(),
    LINE_count = numeric(), 
    LINE_percent = numeric(),
    LTR_count = numeric(),
    LTR_percent = numeric()
  )
  
  for(thresh in thresholds) {
    result <- analysis_results$by_threshold[[thresh]]
    comp <- result$composition
    
    cg_count <- ifelse("CG" %in% names(comp), comp[["CG"]], 0)
    line_count <- ifelse("LINE" %in% names(comp), comp[["LINE"]], 0)
    ltr_count <- ifelse("LTR" %in% names(comp), comp[["LTR"]], 0)
    
    summary_table <- rbind(summary_table, data.frame(
      logFC_threshold = as.numeric(thresh),
      fold_change = round(2^as.numeric(thresh), 2),
      total_genes = result$total,
      pct_of_significant = round(result$percentage_of_significant, 1),
      CG_count = cg_count,
      CG_percent = round(cg_count/result$total*100, 1),
      LINE_count = line_count,
      LINE_percent = round(line_count/result$total*100, 1), 
      LTR_count = ltr_count,
      LTR_percent = round(ltr_count/result$total*100, 1)
    ))
  }
  
  return(summary_table)
}

# Run the analysis
composition_analysis <- analyze_significant_composition(padj_thresh = 0.01)
composition_table <- create_composition_table(composition_analysis)

# Display the table
cat("\n=== COMPOSITION TABLE ===\n")
print(composition_table)

# Save results
write.csv(composition_table, table_path("DE_composition_by_logFC_threshold.csv"), row.names = FALSE)


# ANALYZE THE MOST EXTREME CHANGES (logFC >= 1.5)
# =============================================================================

analyze_extreme_changes <- function(logfc_thresh = 1.5, padj_thresh = 0.01) {
  
  # Get the extreme changes
  extreme_changes <- res_all[res_all$padj < padj_thresh & 
                               abs(res_all$log2FoldChange) >= logfc_thresh & 
                               !is.na(res_all$padj), ]
  
  # Order by effect size (most extreme first)
  extreme_changes <- extreme_changes[order(-abs(extreme_changes$log2FoldChange)), ]
  
  cat(sprintf("=== MOST EXTREME CHANGES (|log2FC| >= %.1f, padj < %.2f) ===\n", 
              logfc_thresh, padj_thresh))
  cat(sprintf("Total: %d genes\n", nrow(extreme_changes)))
  
  # Composition
  composition <- table(extreme_changes$gene_type)
  print(composition)
  
  # Separate by direction and type
  analyze_by_type <- function(gene_type_name) {
    subset_data <- extreme_changes[extreme_changes$gene_type == gene_type_name, ]
    
    if(nrow(subset_data) == 0) {
      cat(sprintf("\nNo %s elements with extreme changes\n", gene_type_name))
      return(NULL)
    }
    
    up_reg <- subset_data[subset_data$log2FoldChange > 0, ]
    down_reg <- subset_data[subset_data$log2FoldChange < 0, ]
    
    cat(sprintf("\n=== %s ELEMENTS (n=%d) ===\n", gene_type_name, nrow(subset_data)))
    cat(sprintf("Upregulated in SCZ: %d\n", nrow(up_reg)))
    cat(sprintf("Downregulated in SCZ: %d\n", nrow(down_reg)))
    
    # Show top upregulated
    if(nrow(up_reg) > 0) {
      cat("\nTop 10 UPREGULATED:\n")
      top_up <- head(up_reg[order(-up_reg$log2FoldChange), ], 10)
      print(data.frame(
        gene = top_up$gene_name,
        log2FC = round(top_up$log2FoldChange, 2),
        fold_change = round(2^top_up$log2FoldChange, 1),
        padj = formatC(top_up$padj, format = "e", digits = 2),
        baseMean = round(top_up$baseMean, 1)
      ))
    }
    
    # Show top downregulated  
    if(nrow(down_reg) > 0) {
      cat("\nTop 10 DOWNREGULATED:\n")
      top_down <- head(down_reg[order(down_reg$log2FoldChange), ], 10)
      print(data.frame(
        gene = top_down$gene_name,
        log2FC = round(top_down$log2FoldChange, 2), 
        fold_change = round(2^abs(top_down$log2FoldChange), 1),
        padj = formatC(top_down$padj, format = "e", digits = 2),
        baseMean = round(top_down$baseMean, 1)
      ))
    }
    
    return(list(all = subset_data, up = up_reg, down = down_reg))
  }
  
  # Analyze each gene type
  cg_extreme <- analyze_by_type("CG")
  line_extreme <- analyze_by_type("LINE") 
  ltr_extreme <- analyze_by_type("LTR")
  
  # Overall summary table
  summary_table <- data.frame(
    gene_name = extreme_changes$gene_name,
    gene_type = extreme_changes$gene_type,
    log2FoldChange = round(extreme_changes$log2FoldChange, 3),
    fold_change = round(2^abs(extreme_changes$log2FoldChange), 2),
    padj = extreme_changes$padj,
    baseMean = round(extreme_changes$baseMean, 1),
    direction = ifelse(extreme_changes$log2FoldChange > 0, "Up in SCZ", "Down in SCZ"),
    abs_log2FC = abs(extreme_changes$log2FoldChange)
  )
  
  # Order by absolute effect size
  summary_table <- summary_table[order(-summary_table$abs_log2FC), ]
  
  return(list(
    all_extreme = extreme_changes,
    summary_table = summary_table,
    cg = cg_extreme,
    line = line_extreme,
    ltr = ltr_extreme
  ))
}


extreme_analysis <- analyze_extreme_changes(logfc_thresh = 1.5, padj_thresh = 0.01)

# Show top 20 most extreme changes overall
cat("\n=== TOP 20 MOST EXTREME CHANGES (ALL TYPES) ===\n")
print(head(extreme_analysis$summary_table[, c("gene_name", "gene_type", "log2FoldChange", "fold_change", "padj")], 20))

# Save the results
write.csv(extreme_analysis$summary_table, 
          table_path("Most_Extreme_DE_Changes_logFC1.5.csv"), row.names = FALSE)

# Quick stats
cat("\n=== QUICK STATS ===\n")
cat(sprintf("Most upregulated: %s (%.1f-fold, %s)\n", 
            extreme_analysis$summary_table$gene_name[1],
            extreme_analysis$summary_table$fold_change[1],
            extreme_analysis$summary_table$gene_type[1]))

most_down <- extreme_analysis$summary_table[which.min(extreme_analysis$summary_table$log2FoldChange), ]
cat(sprintf("Most downregulated: %s (%.1f-fold, %s)\n", 
            most_down$gene_name, most_down$fold_change, most_down$gene_type))



# PLOTS FOR EXTREME CHANGES
# =============================================================================

create_extreme_change_plots <- function(extreme_analysis) {
  
  data <- extreme_analysis$summary_table
  
  # Plot 1: Volcano plot highlighting extreme changes
  volcano_data <- res_all
  volcano_data$is_extreme <- abs(volcano_data$log2FoldChange) >= 1.5 & volcano_data$padj < 0.01
  volcano_data$gene_type <- factor(volcano_data$gene_type, levels = c("CG", "LINE", "LTR"))
  
  p1 <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = gene_type, alpha = is_extreme, size = is_extreme)) +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
    scale_size_manual(values = c("FALSE" = 0.5, "TRUE" = 1.5)) +
    scale_color_manual(values = c("CG" = "steelblue2", "LINE" = "darkorange", "LTR" = "gold3")) +
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", alpha = 0.5) +
    labs(title = "Volcano Plot: Extreme Changes Highlighted",
         subtitle = "Points highlighted: |log2FC| ≥ 1.5 & padj < 0.01",
         x = "log2 Fold Change (SCZ vs Control)",
         y = "-log10(adjusted p-value)",
         color = "Gene Type") +
    guides(alpha = "none", size = "none") +
    theme_minimal()
  
  ggsave(fig_path("volcano_extreme_changes.png"), p1, width = 12, height = 8, dpi = 300)
  
  
  # Plot 2: Bar plot of extreme changes by type and direction
  direction_summary <- data %>%
    group_by(gene_type, direction) %>%
    summarise(count = n(), .groups = "drop")
  
  p2 <- ggplot(direction_summary, aes(x = gene_type, y = count, fill = direction)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.3) +
    scale_fill_manual(values = c("Up in SCZ" = "#E31A1C", "Down in SCZ" = "#1F78B4")) +
    labs(title = "Extreme Expression Changes by Gene Type",
         subtitle = "Genes with |log2FC| ≥ 1.5 and padj < 0.01",
         x = "Gene Type", y = "Number of Genes",
         fill = "Direction") +
    theme_minimal()
  
  ggsave(fig_path("extreme_changes_by_type.png"), p2, width = 10, height = 6, dpi = 300)
  
  
  # Plot 3: Effect size distribution for extreme changes
  p3 <- ggplot(data, aes(x = gene_type, y = abs_log2FC, fill = gene_type)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
    scale_fill_manual(values = c("CG" = "steelblue2", "LINE" = "darkorange", "LTR" = "gold3")) +
    labs(title = "Effect Size Distribution Among Extreme Changes",
         subtitle = "All genes with |log2FC| ≥ 1.5",
         x = "Gene Type", y = "Absolute log2 Fold Change",
         fill = "Gene Type") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave(fig_path("extreme_effect_sizes.png"), p3, width = 8, height = 6, dpi = 300)
  
  
  # Plot 4: Top 20 most extreme changes
  top_20 <- head(data, 20)
  top_20$gene_name <- factor(top_20$gene_name, levels = rev(top_20$gene_name))
  
  p4 <- ggplot(top_20, aes(x = log2FoldChange, y = gene_name, fill = gene_type)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(values = c("CG" = "steelblue2", "LINE" = "darkorange", "LTR" = "gold3")) +
    labs(title = "Top 20 Most Extreme Expression Changes",
         subtitle = "Ordered by absolute log2 fold change",
         x = "log2 Fold Change (SCZ vs Control)",
         y = "Gene/Element Name",
         fill = "Gene Type") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(fig_path("top_20_extreme_changes.png"), p4, width = 12, height = 10, dpi = 300)
  
  
  # Plot 5: Scatter plot of effect size vs significance for extreme changes
  p5 <- ggplot(data, aes(x = abs_log2FC, y = -log10(padj), color = gene_type, shape = direction)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = c("CG" = "steelblue2", "LINE" = "darkorange", "LTR" = "gold3")) +
    scale_shape_manual(values = c("Up in SCZ" = 17, "Down in SCZ" = 19)) +
    labs(title = "Effect Size vs Statistical Significance",
         subtitle = "Extreme changes only (|log2FC| ≥ 1.5)",
         x = "Absolute log2 Fold Change",
         y = "-log10(adjusted p-value)",
         color = "Gene Type", shape = "Direction") +
    theme_minimal()
  
  ggsave(fig_path("extreme_effect_vs_significance.png"), p5, width = 10, height = 6, dpi = 300)
  
  return(list(volcano = p1, bars = p2, boxplot = p3, top20 = p4, scatter = p5))
}

# Create all plots
plots <- create_extreme_change_plots(extreme_analysis)


#  VOLCANO PLOT                  #
###########################################

create_enhanced_volcano <- function(padj_thresh = 0.01, lfc_thresh = 0.5, top_n = 20) {
  
  # Prepare data
  plot_data <- res_all[!is.na(res_all$padj), ]
  
  # Create significance categories
  plot_data$significance <- case_when(
    plot_data$padj < padj_thresh & abs(plot_data$log2FoldChange) > lfc_thresh ~ "Significant",
    plot_data$padj < padj_thresh ~ "Significant (low effect)",
    TRUE ~ "Not significant"
  )
  
  # Create labels for top genes
  # Get top upregulated and downregulated from each gene type
  get_top_genes <- function(data, gene_type, direction, n = 5) {
    subset_data <- data[data$gene_type == gene_type & data$padj < padj_thresh, ]
    if(direction == "up") {
      subset_data <- subset_data[subset_data$log2FoldChange > 0, ]
      subset_data <- subset_data[order(-subset_data$log2FoldChange), ]
    } else {
      subset_data <- subset_data[subset_data$log2FoldChange < 0, ]
      subset_data <- subset_data[order(subset_data$log2FoldChange), ]
    }
    return(head(subset_data$gene_name, n))
  }
  
  # Get top genes to label
  top_genes <- c(
    get_top_genes(plot_data, "CG", "up", 3),
    get_top_genes(plot_data, "CG", "down", 3),
    get_top_genes(plot_data, "LINE", "up", 2),
    get_top_genes(plot_data, "LINE", "down", 5),
    get_top_genes(plot_data, "LTR", "up", 2),
    get_top_genes(plot_data, "LTR", "down", 5)
  )
  
  # ADD: Top 5 most significant genes overall
  sig_data <- plot_data[plot_data$padj < padj_thresh, ]
  most_significant <- head(sig_data[order(sig_data$padj), ]$gene_name, 5)
  
  # Combine all labels
  all_labels <- unique(c(top_genes, most_significant))
  
  # Add label column
  plot_data$label <- ifelse(plot_data$gene_name %in% all_labels, plot_data$gene_name, "")
  
  # Define colors
  colors <- c(
    "CG" = "steelblue2",      # Blue for canonical
    "LINE" = "darkorange",    # Red for LINE  
    "LTR" = "gold3"      # Gold for LTR
  )
  
    # Create the plot
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
    
    # Background points (not significant)
    geom_point(data = plot_data[plot_data$significance == "Not significant", ],
               aes(color = gene_type), alpha = 0.3, size = 0.5) +
    
    # Significant low effect points
    geom_point(data = plot_data[plot_data$significance == "Significant (low effect)", ],
               aes(color = gene_type), alpha = 0.6, size = 0.8) +
    
    # Highly significant points
    geom_point(data = plot_data[plot_data$significance == "Significant", ],
               aes(color = gene_type), alpha = 0.9, size = 1.2) +
    
    # Add threshold lines
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", 
               color = "gray30", alpha = 0.7) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", 
               color = "gray30", alpha = 0.7) +
    
    # Add labels for top genes - BIGGER SIZE
    geom_text_repel(aes(label = label, color = gene_type),
                    size = 3.5,  # Increased from 3 to 3.5
                    box.padding = 0.35,
                    point.padding = 0.3,
                    segment.color = 'gray50',
                    segment.size = 0.3,
                    max.overlaps = 40,  # Increased to handle more labels
                    min.segment.length = 0.1) +
    
    # Customize colors
    scale_color_manual(values = colors,
                       name = "Gene Type",
                       labels = c("CG" = "Canonical Genes", 
                                  "LINE" = "L1 Elements", 
                                  "LTR" = "HERV Elements")) +
    
    # Labels and theme - WHITE BACKGROUND
    labs(
      title = "Differential Expression in Schizophrenia PBMCs",
      subtitle = sprintf("Significant: padj < %g, |log2FC| > %g", padj_thresh, lfc_thresh),
      x = "log2 Fold Change (Schizophrenia vs Control)",
      y = "-log10(adjusted p-value)",
      caption = "Top genes from each category + 5 most significant genes are labeled"
    ) +
    
    theme_classic(base_size = 12) +  # CHANGED to theme_classic for white background
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", size = 0.3)  # Light grid lines
    ) +
    
    # Set axis limits for better visualization
    xlim(c(-4, 4)) +
    ylim(c(0, max(-log10(plot_data$padj), na.rm = TRUE) + 5))
  
  return(p)
}

# Create the volcano plot
volcano_plot <- create_enhanced_volcano(padj_thresh = 0.01, lfc_thresh = 0.5)

# Show statistics
volcano_stats <- create_volcano_stats(padj_thresh = 0.01, lfc_thresh = 0.5)

# Save the plot
ggsave(fig_path("volcano_enhanced_all_types.png"), 
       volcano_plot, 
       width = 8, height = 6, dpi = 300)

# Also save as PDF for publications
ggsave(fig_path("volcano_enhanced_all_types.pdf"), 
       volcano_plot, 
       width = 12, height = 10)

# Print the plot
print(volcano_plot)

################################################################################
# str(dds)           # DESeq2 dataset object
# str(vsd)           # Variance stabilized data  
# str(res_all)       # Combined differential expression results
# str(samples_metadata) # Sample information
# 
# # Check dimensions
# dim(dds)           # genes x samples
# dim(assay(vsd))    # genes x samples (VST-transformed counts)
# dim(res_all)       # genes x DE statistics
# 
# # Raw count data
# head(counts(dds)[1:5, 1:5])  # Raw counts per gene per sample
# colData(dds)                 # Sample metadata (phenotype, etc.)
# 
# # VST-transformed counts (log2-like, variance stabilized)
# head(assay(vsd)[1:5, 1:5])   # These are the values we'll plot
# # Range should be roughly 0-15 (like log2 counts)
# 
# head(res_all)
# # Contains: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, gene_type
# # This tells us WHICH genes are significant
# 
# head(samples_metadata)
# # Contains: sample IDs, phenotype (control/schizophrenia), other covariates
# 
# 
# 
# # Find significant genes
# res_sig_cg <- res_all[res_all$gene_type == "CG" & !is.na(res_all$padj) & res_all$padj < 0.01, ]
# # This gets only canonical genes with padj < 0.01
# 
# # Get top genes by significance  
# top_cg <- head(rownames(res_sig_cg[order(res_sig_cg$padj), ]), 20)
# # This gets the 20 most significant canonical gene names
# 
# mat <- assay(vsd)[all_top_genes, ]
# # pull VST counts for selected genes across all samples
# # Rows genes, Cols samples
# 
# # Z-score normalize
# mat_zscore <- t(scale(t(mat)))
# # For each gene (row): subtract gene's mean, divide by gene's SD
# # so that each gene has mean=0, sd=1 across samples
# 
# ####
# # top Variable Genes
# gene_vars <- apply(assay(vsd)[rownames(res_sig_cg), ], 1, var)
# top_cg <- head(names(sort(gene_vars, decreasing = TRUE)), 20)

################################################################################
create_top50_heatmap <- function(res_all, vsd_data, method = "combined_score", n_genes = 50) {
  
  # Get all significant genes
  res_sig <- res_all[!is.na(res_all$padj) & res_all$padj < 0.01, ]
  
  if (nrow(res_sig) < n_genes) {
    cat(sprintf("Only %d significant genes found, using all of them\n", nrow(res_sig)))
    n_genes <- nrow(res_sig)
  }
  
  # Select top genes by chosen method
  if (method == "combined_score") {
    res_sig$combined_score <- -log10(res_sig$padj) * abs(res_sig$log2FoldChange)
    top_genes <- head(rownames(res_sig[order(-res_sig$combined_score), ]), n_genes)
    title_suffix <- "by Combined Score (Significance × Effect Size)"
    filename_suffix <- "combined_score"
  } else if (method == "variable") {
    gene_vars <- apply(assay(vsd_data)[rownames(res_sig), ], 1, var)
    top_genes <- head(names(sort(gene_vars, decreasing = TRUE)), n_genes)
    title_suffix <- "by Variability Across Samples"
    filename_suffix <- "most_variable"
  } else if (method == "pvalue") {
    top_genes <- head(rownames(res_sig[order(res_sig$padj), ]), n_genes)
    title_suffix <- "by Significance (p-value)"
    filename_suffix <- "most_significant"
  }
  
  # Show breakdown by gene type (LTR shown as HERV)
  gene_types_raw <- table(res_all[top_genes, "gene_type"])
  cat(sprintf("Selected top %d genes %s:\n", length(top_genes), title_suffix))
  cat(sprintf("  Canonical: %d, L1: %d, HERV: %d\n", 
              gene_types_raw["CG"], gene_types_raw["LINE"], gene_types_raw["LTR"]))
  
  # Get VST counts and Z-score normalize
  mat <- assay(vsd_data)[top_genes, ]
  mat_zscore <- t(scale(t(mat)))
  
  # Create annotations with HERV labeling
  gene_types_for_heatmap <- res_all[top_genes, "gene_type"]
  # Convert LTR to HERV for display
  gene_types_display <- factor(gene_types_for_heatmap, 
                               levels = c("CG", "LINE", "LTR"),
                               labels = c("CG", "LINE", "HERV"))
  
  row_annotation <- data.frame(
    Gene_Type = gene_types_display,
    row.names = top_genes
  )
  
  col_annotation <- data.frame(
    Phenotype = samples_metadata$phenotype,
    row.names = colnames(mat_zscore)
  )
  
  # Update colors to use HERV instead of LTR
  ann_colors <- list(
    Phenotype = phenotype_colors,
    Gene_Type = c("CG" = gene_type_colors[["CG"]], 
                  "LINE" = gene_type_colors[["LINE"]], 
                  "HERV" = gene_type_colors[["LTR"]])  # Use LTR color for HERV
  )
  p <- pheatmap(mat_zscore,
                annotation_col = col_annotation,
                annotation_row = row_annotation,
                annotation_colors = ann_colors,
                show_colnames = FALSE,
                show_rownames = (length(top_genes) <= 50),
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                main = sprintf("Top %d DE Features %s\n(%d Canonical, %d L1, %d HERV)", 
                               length(top_genes), title_suffix,
                               gene_types_raw["CG"], gene_types_raw["LINE"], gene_types_raw["LTR"]))
  
  # Create and save heatmap
  png(fig_path(paste0("heatmap_top50_", filename_suffix, ".png")), 
      width = 12, height = 10, units = "in", res = 300)
  
  print(p)
 
  dev.off()
  
  cat(sprintf("Heatmap saved: heatmap_top50_%s.png\n", filename_suffix))
  return(p)
}


# Top 50 by combined score (significance × effect size)
p1 <- create_top50_heatmap(res_all, vsd, method = "combined_score", n_genes = 50)






###########################################
# BOXPLOT FOR SPECIFIC GENE/FEATURE      #
###########################################

plot_gene_expression <- function(gene_name, vsd_data, res_all, samples_metadata) {
  
  # Check if gene exists
  if (!gene_name %in% rownames(vsd_data)) {
    cat("Gene", gene_name, "not found in data\n")
    return(NULL)
  }
  
  # Get VST expression data for the gene
  gene_expr <- assay(vsd_data)[gene_name, ]
  
  # Get gene info from results
  gene_info <- res_all[gene_name, ]
  gene_type <- gene_info$gene_type
  log2fc <- round(gene_info$log2FoldChange, 2)
  padj <- gene_info$padj
  
  # Format p-value
  if (padj < 0.001) {
    padj_text <- sprintf("p < 0.001")
  } else {
    padj_text <- sprintf("p = %.3f", padj)
  }
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Expression = gene_expr,
    Phenotype = samples_metadata$phenotype,
    Sample = names(gene_expr)
  )
  
  # Convert gene type for display
  gene_type_display <- ifelse(gene_type == "LTR", "HERV", 
                              ifelse(gene_type == "LINE", "L1", 
                                     ifelse(gene_type == "CG", "Canonical", gene_type)))
  
  # Create boxplot
  library(ggplot2)
  
  p <- ggplot(plot_data, aes(x = Phenotype, y = Expression, fill = Phenotype)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +
    scale_fill_manual(values = phenotype_colors) +
    labs(
      title = sprintf("%s (%s)", gene_name, gene_type_display),
      subtitle = sprintf("log2FC = %s, %s", log2fc, padj_text),
      x = "Phenotype",
      y = "VST Expression",
      caption = "Points represent individual samples"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "none"
    ) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black")
  
  # Add significance stars
  if (padj < 0.001) {
    sig_label <- "***"
  } else if (padj < 0.01) {
    sig_label <- "**"
  } else if (padj < 0.05) {
    sig_label <- "*"
  } else {
    sig_label <- "ns"
  }
  
  # Add significance annotation
  y_max <- max(plot_data$Expression) + 0.1 * diff(range(plot_data$Expression))
  p <- p + annotate("text", x = 1.5, y = y_max, label = sig_label, size = 6, fontface = "bold")
  
  return(p)
}

###########################################
# PLOT MULTIPLE GENES AT ONCE            #
###########################################

plot_multiple_genes <- function(gene_list, vsd_data, res_all, samples_metadata, ncol = 2) {
  
  library(gridExtra)
  
  # Create list of plots
  plot_list <- lapply(gene_list, function(gene) {
    plot_gene_expression(gene, vsd_data, res_all, samples_metadata)
  })
  
  # Remove NULL plots (genes not found)
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  if (length(plot_list) == 0) {
    cat("No valid genes found\n")
    return(NULL)
  }
  
  # Arrange plots
  grid.arrange(grobs = plot_list, ncol = ncol)
}

###########################################

# Plot single gene
#p1 <- plot_gene_expression("CD83", vsd, res_all, samples_metadata)
#print(p1)

# Save single plot
#ggsave(fig_path("boxplot_CD83.png"), p1, width = 6, height = 5, dpi = 300)

# Plot multiple interesting genes
df_res_all_sig <- as.data.frame(res_all_significant)
df_res_all_sig <- df_res_all_sig[order(df_res_all_sig$padj), ]

interesting_genes <- head(df_res_all_sig$gene_name, 6)

# Create and save multi-panel plot
png(fig_path("boxplot_top_genes.png"), width = 12, height = 10, units = "in", res = 300)
plot_multiple_genes(interesting_genes, vsd, res_all, samples_metadata, ncol = 3)
dev.off()

interesting_tes <- head(df_res_all_sig[df_res_all_sig$gene_type %in% c("LINE", "LTR"), "gene_name"], 6)

# Create and save multi-panel plot
png(fig_path("boxplot_top_tes.png"), width = 12, height = 10, units = "in", res = 300)
plot_multiple_genes(interesting_tes, vsd, res_all, samples_metadata, ncol = 3)
dev.off()

# Create and save multi-panel plot
png(fig_path("boxplot_GJB.png"), width = 12, height = 10, units = "in", res = 300)
plot_multiple_genes(c("GJB6", "GJB2"), vsd, res_all, samples_metadata, ncol = 2)
dev.off()

# just specify any gene name
gene_of_interest <- "CD83"  
p <- plot_gene_expression(gene_of_interest, vsd, res_all, samples_metadata)
print(p)
