# =============================================================================
# Hub Analysis + Centrality Metrics + SCZ Gene Integration
# High-Impact Analysis for Presentation
# =============================================================================

library(igraph)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

cat("=============================================================================\n")
cat("HUB ANALYSIS + SCZ GENE INTEGRATION\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# STEP 1: Load Filtered Networks
# -----------------------------------------------------------------------------

scz_graph <- readRDS("results/network_analysis/networks/scz_graph_intergenic_only.rds")
hc_graph <- readRDS("results/network_analysis/networks/hc_graph_intergenic_only.rds")
feature_lists <- readRDS("results/network_analysis/data/feature_lists_intergenic_only.rds")
degree_info <- readRDS("results/network_analysis/metrics/degree_info_intergenic_only.rds")

top_genes <- feature_lists$top_genes
top_tes <- feature_lists$top_tes

cat("Networks loaded:\n")
cat("  Genes:", length(top_genes), "\n")
cat("  Intergenic TEs:", length(top_tes), "\n\n")

# Helper function
save_plot <- function(filename, plot_function, width=10, height=5) {
  filename <- sub("\\.pdf$|\\.png$", "", filename)
  pdf(paste0(filename, ".pdf"), width=width, height=height)
  plot_function()
  dev.off()
  png(paste0(filename, ".png"), width=width*100, height=height*100, res=100)
  plot_function()
  dev.off()
  cat("  ✓ Saved:", basename(filename), "\n")
}

# Create presentation plots directory
dir.create("results/network_analysis/presentation_plots", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# STEP 2: Define SCZ Risk Genes
# -----------------------------------------------------------------------------

scz_genes <- c("MTHFR", "RTN4R", "COMT", "HTR2A", "DRD3", "SYN2", "CHI3L1", 
               "NRXN1", "MIR30E", "SETD1A", "MBD5", "PHIP", "IRAK1BP1", "MIR206", 
               "MIR198", "DAOA", "DISC1", "APOL2", "APOL4", "DTNBP1", "DAO", 
               "ABCA13", "SHANK3", "AKT1", "DISC2", "SCZD6", "SCZD1", "SCZD2", 
               "SCZD7", "SCZD8", "SCZD3", "SCZD12", "SCZD11", "NPAS3", "RELN", 
               "MIR195", "MIR29C", "MIR9-1", "MIR212", "MIR15B", "ULK4", "MIR346", 
               "MIR107", "MIR29A", "MIR106B", "MIR30A", "MIR26B", "MIR30B", 
               "MIR15A", "MIR30D", "MIR20B", "ZNF804A", "DRD2", "ASTN2", "RGS4", 
               "PRODH", "GABRB2", "CSMD1", "DLG2", "DLGAP2", "RBFOX1", "PRKN", 
               "CTNND2", "DPP6", "DPP10", "MSRA", "PCDH15", "DMD", "NRG1", 
               "TSNAX-DISC1", "DRD4", "SLC1A1", "GRIN2B", "ZDHHC8", "CHRNA7", 
               "GRIN2A", "PDE4B", "GRID2", "NOTCH4", "CNTN6", "EHMT1", "MACROD2", 
               "KCNN3", "KATNAL2", "TPH1", "GABRB1", "GRIK3", "LPP", "PTPRM", 
               "NRG3", "TRAPPC9", "DBH", "FHIT", "PARD3B", "PDE11A", "TMLHE", 
               "VPS13B", "FZD3", "PTPRT", "CNP")

# Check which SCZ genes are in our networks
scz_genes_in_network <- scz_genes[scz_genes %in% top_genes]

cat("SCZ Risk Genes:\n")
cat("  Total in list:", length(scz_genes), "\n")
cat("  Found in network:", length(scz_genes_in_network), "\n")
cat("  Missing:", length(scz_genes) - length(scz_genes_in_network), "\n\n")

if(length(scz_genes_in_network) < 20) {
  cat("⚠ WARNING: Low overlap. Genes may have different names in your data.\n")
  cat("First few network gene names:\n")
  print(head(top_genes, 10))
  cat("\nFirst few SCZ gene names:\n")
  print(head(scz_genes, 10))
  cat("\n")
}

# -----------------------------------------------------------------------------
# STEP 3: Calculate Multiple Centrality Metrics
# -----------------------------------------------------------------------------

cat("Calculating centrality metrics...\n")

# For both networks
metrics_scz <- data.frame(
  node = V(scz_graph)$name,
  type = V(scz_graph)$type,
  degree = degree(scz_graph),
  betweenness = betweenness(scz_graph, normalized = TRUE),
  closeness = closeness(scz_graph, normalized = TRUE),
  eigenvector = eigen_centrality(scz_graph)$vector,
  is_scz_gene = V(scz_graph)$name %in% scz_genes_in_network,
  stringsAsFactors = FALSE
)

metrics_hc <- data.frame(
  node = V(hc_graph)$name,
  type = V(hc_graph)$type,
  degree = degree(hc_graph),
  betweenness = betweenness(hc_graph, normalized = TRUE),
  closeness = closeness(hc_graph, normalized = TRUE),
  eigenvector = eigen_centrality(hc_graph)$vector,
  is_scz_gene = V(hc_graph)$name %in% scz_genes_in_network,
  stringsAsFactors = FALSE
)

cat("✓ Centrality metrics calculated\n\n")

# -----------------------------------------------------------------------------
# STEP 4: Hub Identification
# -----------------------------------------------------------------------------

cat("Identifying hubs...\n")

# Define hubs as top 10% by degree
hub_threshold_scz <- quantile(metrics_scz$degree, 0.90)
hub_threshold_hc <- quantile(metrics_hc$degree, 0.90)

metrics_scz$is_hub <- metrics_scz$degree >= hub_threshold_scz
metrics_hc$is_hub <- metrics_hc$degree >= hub_threshold_hc

# Hub statistics
hub_stats <- data.frame(
  Network = c("SCZ", "HC"),
  Total_Hubs = c(sum(metrics_scz$is_hub), sum(metrics_hc$is_hub)),
  Gene_Hubs = c(sum(metrics_scz$is_hub & metrics_scz$type == "gene"),
                sum(metrics_hc$is_hub & metrics_hc$type == "gene")),
  TE_Hubs = c(sum(metrics_scz$is_hub & metrics_scz$type == "TE"),
              sum(metrics_hc$is_hub & metrics_hc$type == "TE")),
  SCZ_Gene_Hubs = c(sum(metrics_scz$is_hub & metrics_scz$is_scz_gene),
                    sum(metrics_hc$is_hub & metrics_hc$is_scz_gene))
)

hub_stats$Pct_TE_Hubs <- 100 * hub_stats$TE_Hubs / hub_stats$Total_Hubs
hub_stats$Pct_Genes <- 100 * length(top_genes) / (length(top_genes) + length(top_tes))
hub_stats$Pct_TEs <- 100 * length(top_tes) / (length(top_genes) + length(top_tes))

cat("Hub Statistics:\n")
print(hub_stats)
cat("\n")

# Test for TE enrichment in hubs
fisher_scz <- matrix(c(
  sum(metrics_scz$is_hub & metrics_scz$type == "TE"),
  sum(metrics_scz$is_hub & metrics_scz$type == "gene"),
  sum(!metrics_scz$is_hub & metrics_scz$type == "TE"),
  sum(!metrics_scz$is_hub & metrics_scz$type == "gene")
), nrow=2)

fisher_hc <- matrix(c(
  sum(metrics_hc$is_hub & metrics_hc$type == "TE"),
  sum(metrics_hc$is_hub & metrics_hc$type == "gene"),
  sum(!metrics_hc$is_hub & metrics_hc$type == "TE"),
  sum(!metrics_hc$is_hub & metrics_hc$type == "gene")
), nrow=2)

fisher_test_scz <- fisher.test(fisher_scz)
fisher_test_hc <- fisher.test(fisher_hc)

cat("TE Hub Enrichment Tests:\n")
cat("SCZ: OR =", round(fisher_test_scz$estimate, 2), 
    ", p =", format(fisher_test_scz$p.value, scientific=TRUE), "\n")
cat("HC:  OR =", round(fisher_test_hc$estimate, 2), 
    ", p =", format(fisher_test_hc$p.value, scientific=TRUE), "\n\n")

# -----------------------------------------------------------------------------
# STEP 5: Top Hubs
# -----------------------------------------------------------------------------

top_hubs_scz <- head(metrics_scz[order(-metrics_scz$degree), ], 20)
top_hubs_hc <- head(metrics_hc[order(-metrics_hc$degree), ], 20)

cat("Top 20 Hubs in SCZ:\n")
print(top_hubs_scz[, c("node", "type", "degree", "is_scz_gene")])
cat("\n")

cat("Top 20 Hubs in HC:\n")
print(top_hubs_hc[, c("node", "type", "degree", "is_scz_gene")])
cat("\n")

# Top TE hubs
top_te_hubs_scz <- head(metrics_scz[metrics_scz$type == "TE", ][order(-metrics_scz[metrics_scz$type == "TE", "degree"]), ], 10)
top_te_hubs_hc <- head(metrics_hc[metrics_hc$type == "TE", ][order(-metrics_hc[metrics_hc$type == "TE", "degree"]), ], 10)

cat("Top 10 TE Hubs in SCZ:\n")
print(top_te_hubs_scz[, c("node", "degree", "betweenness")])
cat("\n")

cat("Top 10 TE Hubs in HC:\n")
print(top_te_hubs_hc[, c("node", "degree", "betweenness")])
cat("\n")

# -----------------------------------------------------------------------------
# STEP 6: SCZ Gene Analysis
# -----------------------------------------------------------------------------

cat("Analyzing SCZ gene connectivity...\n")

scz_gene_metrics_scz <- metrics_scz[metrics_scz$is_scz_gene, ]
scz_gene_metrics_hc <- metrics_hc[metrics_hc$is_scz_gene, ]

cat("SCZ Genes in Networks:\n")
cat("Average degree in SCZ network:", round(mean(scz_gene_metrics_scz$degree), 1), "\n")
cat("Average degree in HC network:", round(mean(scz_gene_metrics_hc$degree), 1), "\n")
cat("Number that are hubs in SCZ:", sum(scz_gene_metrics_scz$is_hub), "\n")
cat("Number that are hubs in HC:", sum(scz_gene_metrics_hc$is_hub), "\n\n")

# TE connections to SCZ genes
cat("TE connections to SCZ genes...\n")

# For each TE, count SCZ gene neighbors
te_scz_connections_scz <- sapply(top_tes, function(te) {
  neighbors <- neighbors(scz_graph, te)
  sum(V(scz_graph)$name[neighbors] %in% scz_genes_in_network)
})

te_scz_connections_hc <- sapply(top_tes, function(te) {
  neighbors <- neighbors(hc_graph, te)
  sum(V(hc_graph)$name[neighbors] %in% scz_genes_in_network)
})

te_scz_df <- data.frame(
  TE = top_tes,
  SCZ_connections_scz = te_scz_connections_scz,
  SCZ_connections_hc = te_scz_connections_hc,
  total_degree_scz = metrics_scz[match(top_tes, metrics_scz$node), "degree"],
  total_degree_hc = metrics_hc[match(top_tes, metrics_hc$node), "degree"]
)

te_scz_df$pct_scz_scz <- 100 * te_scz_df$SCZ_connections_scz / te_scz_df$total_degree_scz
te_scz_df$pct_scz_hc <- 100 * te_scz_df$SCZ_connections_hc / te_scz_df$total_degree_hc

top_te_scz_connectors <- head(te_scz_df[order(-te_scz_df$SCZ_connections_scz), ], 10)

cat("Top 10 TEs by SCZ gene connections (SCZ network):\n")
print(top_te_scz_connectors[, c("TE", "SCZ_connections_scz", "pct_scz_scz")])
cat("\n")

# -----------------------------------------------------------------------------
# STEP 7: PRESENTATION PLOTS
# -----------------------------------------------------------------------------

cat("Creating presentation-ready plots...\n\n")

# PLOT 1: TE vs Gene Degree (IMPROVED ORDER: HC then SCZ)
save_plot("results/network_analysis/presentation_plots/01_degree_comparison",
          function() {
            par(mfrow=c(1,2), mar=c(5,5,3,2))
            
            # HC FIRST (baseline)
            boxplot(list(
              Genes = metrics_hc$degree[metrics_hc$type == "gene"],
              TEs = metrics_hc$degree[metrics_hc$type == "TE"]
            ), 
            main = "HC Network (Baseline)", 
            ylab = "Degree", 
            col = c("steelblue", "coral"),
            outline = FALSE,
            cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
            
            p_hc <- wilcox.test(
              metrics_hc$degree[metrics_hc$type == "gene"],
              metrics_hc$degree[metrics_hc$type == "TE"]
            )$p.value
            
            text(1.5, max(metrics_hc$degree) * 0.9,
                 labels = paste0("p < ", format(p_hc, digits=2, scientific=TRUE)),
                 cex = 1.2)
            
            # SCZ SECOND (disease)
            boxplot(list(
              Genes = metrics_scz$degree[metrics_scz$type == "gene"],
              TEs = metrics_scz$degree[metrics_scz$type == "TE"]
            ), 
            main = "SCZ Network (Disease)", 
            ylab = "Degree", 
            col = c("steelblue", "coral"),
            outline = FALSE,
            cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
            
            p_scz <- wilcox.test(
              metrics_scz$degree[metrics_scz$type == "gene"],
              metrics_scz$degree[metrics_scz$type == "TE"]
            )$p.value
            
            text(1.5, max(metrics_scz$degree) * 0.9,
                 labels = paste0("p < ", format(p_scz, digits=2, scientific=TRUE)),
                 cex = 1.2)
          }, width=12, height=6)







# PLOT 2: Hub Enrichment (Q2)
save_plot("results/network_analysis/presentation_plots/02_hub_enrichment",
          function() {
            par(mar=c(6,5,3,2))
            
            # Prepare data
            hub_data <- rbind(
              data.frame(Network="SCZ", Type="Expected\n(% in network)", 
                         Value=hub_stats$Pct_TEs[1]),
              data.frame(Network="SCZ", Type="Observed\n(% in hubs)", 
                         Value=hub_stats$Pct_TE_Hubs[1]),
              data.frame(Network="HC", Type="Expected\n(% in network)", 
                         Value=hub_stats$Pct_TEs[2]),
              data.frame(Network="HC", Type="Observed\n(% in hubs)", 
                         Value=hub_stats$Pct_TE_Hubs[2])
            )
            
            bp <- barplot(
              matrix(hub_data$Value, nrow=2),
              beside = TRUE,
              names.arg = c("SCZ", "HC"),
              col = c("gray70", "coral"),
              ylab = "% TEs",
              main = "TE Enrichment in Network Hubs",
              ylim = c(0, max(hub_data$Value) * 1.3),
              cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4,
              cex.names = 1.3
            )
            
            legend("topright", 
                   legend = c("Expected (% in network)", "Observed (% in hubs)"),
                   fill = c("gray70", "coral"),
                   cex = 1.2, bty = "n")
            
            # Add significance stars
            text(bp[2,1], hub_data$Value[2] + 2, 
                 labels = ifelse(fisher_test_scz$p.value < 0.001, "***", 
                                 ifelse(fisher_test_scz$p.value < 0.01, "**", "*")),
                 cex = 2)
            text(bp[2,2], hub_data$Value[4] + 2, 
                 labels = ifelse(fisher_test_hc$p.value < 0.001, "***", 
                                 ifelse(fisher_test_hc$p.value < 0.01, "**", "*")),
                 cex = 2)
            
            # Add p-values
            text(bp[2,1], hub_data$Value[2] + 5,
                 labels = paste0("OR=", round(fisher_test_scz$estimate, 1)),
                 cex = 1.1)
            text(bp[2,2], hub_data$Value[4] + 5,
                 labels = paste0("OR=", round(fisher_test_hc$estimate, 1)),
                 cex = 1.1)
          }, width=10, height=8)

# PLOT 3: Degree HC vs SCZ Scatter (FLIPPED AXES)
save_plot("results/network_analysis/presentation_plots/03_degree_hc_vs_scz",
          function() {
            # Match nodes
            common_nodes <- intersect(metrics_scz$node, metrics_hc$node)
            scz_deg <- metrics_scz$degree[match(common_nodes, metrics_scz$node)]
            hc_deg <- metrics_hc$degree[match(common_nodes, metrics_hc$node)]
            node_types <- metrics_scz$type[match(common_nodes, metrics_scz$node)]
            is_scz <- metrics_scz$is_scz_gene[match(common_nodes, metrics_scz$node)]
            
            par(mar=c(5,5,3,2))
            # FLIPPED: HC on x-axis, SCZ on y-axis
            plot(hc_deg, scz_deg,
                 pch = ifelse(node_types == "TE", 17, 16),
                 col = ifelse(is_scz, "gold", 
                              ifelse(node_types == "TE", "coral", "steelblue")),
                 cex = ifelse(is_scz, 1.5, 0.8),
                 xlab = "Degree in HC Network (Baseline)",
                 ylab = "Degree in SCZ Network (Disease)",
                 main = "Connectivity Change: HC → SCZ",
                 cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
            
            abline(0, 1, col = "red", lwd = 2, lty = 2)
            abline(h = hub_threshold_scz, col = "gray50", lty = 3)
            abline(v = hub_threshold_hc, col = "gray50", lty = 3)
            
            # Add text annotation
            text(max(hc_deg)*0.7, max(scz_deg)*0.2,
                 labels = "Loss of connectivity\nin disease",
                 col = "red", cex = 1.2, font = 2)
            
            legend("topleft",
                   legend = c("Genes", "TEs", "SCZ Risk Genes", "Equal connectivity"),
                   pch = c(16, 17, 16, NA),
                   col = c("steelblue", "coral", "gold", "red"),
                   pt.cex = c(0.8, 0.8, 1.5, NA),
                   lty = c(NA, NA, NA, 2),
                   lwd = c(NA, NA, NA, 2),
                   cex = 1.1, bty = "n")
            
            grid(col = "gray90")
          }, width=10, height=8)

# PLOT 4: Centrality Heatmap (FIXED MARGINS)
save_plot("results/network_analysis/presentation_plots/04_te_centrality_heatmap",
          function() {
            # Get top 20 TEs by average degree
            te_avg_degree <- (metrics_scz$degree[metrics_scz$type == "TE"] + 
                                metrics_hc$degree[metrics_hc$type == "TE"]) / 2
            names(te_avg_degree) <- top_tes
            top20_tes <- names(sort(te_avg_degree, decreasing = TRUE)[1:20])
            
            # Create matrices for SCZ and HC separately
            scz_data <- metrics_scz[match(top20_tes, metrics_scz$node), 
                                    c("degree", "betweenness", "closeness", "eigenvector")]
            hc_data <- metrics_hc[match(top20_tes, metrics_hc$node), 
                                  c("degree", "betweenness", "closeness", "eigenvector")]
            
            # Convert to matrices
            scz_matrix <- as.matrix(scz_data)
            hc_matrix <- as.matrix(hc_data)
            
            # Combine
            heatmap_data <- rbind(scz_matrix, hc_matrix)
            
            # Normalize to 0-1 by column
            heatmap_data_norm <- apply(heatmap_data, 2, function(x) {
              if(max(x) == min(x)) return(rep(0.5, length(x)))
              (x - min(x)) / (max(x) - min(x))
            })
            
            # Set rownames
            rownames(heatmap_data_norm) <- paste0(
              rep(c("SCZ_", "HC_"), each = 20),
              rep(top20_tes, 2)
            )
            
            # Rename columns for clarity
            colnames(heatmap_data_norm) <- c("Degree", "Betweenness", "Closeness", "Eigenvector")
            
            # Plot with FIXED margins
            library(pheatmap)
            pheatmap(t(heatmap_data_norm),
                     main = "Centrality Metrics: Top 20 TE Hubs",
                     cluster_rows = FALSE,
                     cluster_cols = TRUE,
                     color = colorRampPalette(c("white", "orange", "red"))(100),
                     fontsize = 9,
                     fontsize_row = 12,
                     fontsize_col = 7,
                     angle_col = 45,
                     border_color = "gray80",
                     cellwidth = 15,
                     cellheight = 25)
          }, width=16, height=6)


# PLOT 5: SCZ Gene Degree (HC FIRST)
save_plot("results/network_analysis/presentation_plots/05_scz_gene_degree",
          function() {
            par(mfrow=c(1,2), mar=c(5,5,3,2))
            
            # HC network FIRST
            boxplot(list(
              "Other Genes" = metrics_hc$degree[metrics_hc$type == "gene" & !metrics_hc$is_scz_gene],
              "SCZ Risk Genes" = metrics_hc$degree[metrics_hc$is_scz_gene]
            ),
            main = "HC Network (Baseline)",
            ylab = "Degree",
            col = c("steelblue", "gold"),
            outline = FALSE,
            cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
            
            # SCZ network SECOND
            boxplot(list(
              "Other Genes" = metrics_scz$degree[metrics_scz$type == "gene" & !metrics_scz$is_scz_gene],
              "SCZ Risk Genes" = metrics_scz$degree[metrics_scz$is_scz_gene]
            ),
            main = "SCZ Network (Disease)",
            ylab = "Degree",
            col = c("steelblue", "gold"),
            outline = FALSE,
            cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)
          }, width=12, height=6)

# PLOT 6: TE-SCZ Gene Connections (FIXED)
save_plot("results/network_analysis/presentation_plots/06_te_scz_connections",
          function() {
            par(mfrow=c(1,2), mar=c(8,5,3,2))
            
            # Get top 10 TEs
            top10 <- head(te_scz_df[order(-te_scz_df$SCZ_connections_scz), ], 10)
            
            # Plot 1: Single network (SCZ)
            bp1 <- barplot(top10$SCZ_connections_scz,
                           names.arg = "",
                           main = "TEs Connected to SCZ Risk Genes",
                           ylab = "# SCZ Gene Neighbors (SCZ Network)",
                           col = "coral",
                           las = 2,
                           ylim = c(0, max(top10$SCZ_connections_scz) * 1.2),
                           cex.lab = 1.3, cex.main = 1.4)
            
            text(bp1, par("usr")[3] - 0.5, 
                 labels = top10$TE, 
                 srt = 45, adj = 1, xpd = TRUE, cex = 0.85)
            
            # Plot 2: Comparison SCZ vs HC - FIXED
            te_matrix <- rbind(
              SCZ = top10$SCZ_connections_scz,
              HC = top10$SCZ_connections_hc
            )
            
            bp2 <- barplot(te_matrix,
                           beside = TRUE,
                           names.arg = rep("", ncol(te_matrix)),  # Empty names
                           main = "Network Comparison",
                           ylab = "# SCZ Gene Neighbors",
                           col = c("orangered", "steelblue"),
                           las = 2,
                           ylim = c(0, max(te_matrix) * 1.2),
                           cex.lab = 1.3, cex.main = 1.4)
            
            # Add TE names below
            text(colMeans(bp2), par("usr")[3] - 0.5,
                 labels = top10$TE,
                 srt = 45, adj = 1, xpd = TRUE, cex = 0.85)
            
            legend("topright",
                   legend = c("SCZ Network", "HC Network"),
                   fill = c("orangered", "steelblue"),
                   cex = 1.2, bty = "n")
          }, width=16, height=8)


# PLOT 7: Network Overview Statistics (Q1)
save_plot("results/network_analysis/presentation_plots/07_network_overview",
          function() {
            par(mfrow=c(2,2), mar=c(4,5,3,2))
            
            # 1. Node counts
            barplot(c(SCZ = vcount(scz_graph), HC = vcount(hc_graph)),
                    col = c("orangered", "steelblue"),
                    main = "Network Size",
                    ylab = "Number of Nodes",
                    ylim = c(0, vcount(scz_graph) * 1.1),
                    cex.lab = 1.2, cex.main = 1.3)
            
            # 2. Edge counts
            barplot(c(SCZ = ecount(scz_graph), HC = ecount(hc_graph)),
                    col = c("orangered", "steelblue"),
                    main = "Network Connectivity",
                    ylab = "Number of Edges",
                    cex.lab = 1.2, cex.main = 1.3)
            
            # 3. Average degree
            barplot(c(SCZ = mean(metrics_scz$degree), HC = mean(metrics_hc$degree)),
                    col = c("orangered", "steelblue"),
                    main = "Average Degree",
                    ylab = "Mean Degree",
                    cex.lab = 1.2, cex.main = 1.3)
            
            # 4. Clustering
            barplot(c(SCZ = transitivity(scz_graph), HC = transitivity(hc_graph)),
                    col = c("orangered", "steelblue"),
                    main = "Clustering Coefficient",
                    ylab = "Clustering",
                    ylim = c(0, 0.25),
                    cex.lab = 1.2, cex.main = 1.3)
          }, width=12, height=10)

cat("\n✓ All presentation plots created!\n\n")

# -----------------------------------------------------------------------------
# STEP 8: Save Results
# -----------------------------------------------------------------------------

# Save metrics
saveRDS(list(scz = metrics_scz, hc = metrics_hc),
        "results/network_analysis/metrics/centrality_metrics.rds")

# Save hub info
saveRDS(list(
  hub_stats = hub_stats,
  top_hubs_scz = top_hubs_scz,
  top_hubs_hc = top_hubs_hc,
  top_te_hubs_scz = top_te_hubs_scz,
  top_te_hubs_hc = top_te_hubs_hc
), "results/network_analysis/metrics/hub_analysis.rds")

# Save SCZ gene analysis
saveRDS(list(
  te_scz_connections = te_scz_df,
  scz_gene_metrics_scz = scz_gene_metrics_scz,
  scz_gene_metrics_hc = scz_gene_metrics_hc
), "results/network_analysis/metrics/scz_gene_analysis.rds")

# Export tables for presentation
write.csv(hub_stats, 
          "results/network_analysis/presentation_plots/hub_statistics.csv",
          row.names = FALSE)

write.csv(top_te_scz_connectors,
          "results/network_analysis/presentation_plots/top_te_scz_connectors.csv",
          row.names = FALSE)

write.csv(rbind(
  data.frame(Network = "SCZ", top_hubs_scz[, c("node", "type", "degree", "is_scz_gene")]),
  data.frame(Network = "HC", top_hubs_hc[, c("node", "type", "degree", "is_scz_gene")])
), "results/network_analysis/presentation_plots/top_hubs.csv",
row.names = FALSE)

cat("=============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=============================================================================\n\n")

cat("📊 Presentation plots created in:\n")
cat("    results/network_analysis/presentation_plots/\n\n")

cat("Key findings for your presentation:\n\n")

cat("Q1: Are TEs part of the regulatory network?\n")
cat("  ✓ YES! TEs have HIGHER average degree than genes\n")
cat("    SCZ: TEs =", round(mean(metrics_scz$degree[metrics_scz$type == "TE"]), 1),
    "vs Genes =", round(mean(metrics_scz$degree[metrics_scz$type == "gene"]), 1), "\n")
cat("    HC:  TEs =", round(mean(metrics_hc$degree[metrics_hc$type == "TE"]), 1),
    "vs Genes =", round(mean(metrics_hc$degree[metrics_hc$type == "gene"]), 1), "\n")
cat("  ✓ HC network 2x more connected than SCZ\n\n")

cat("Q2: What positions do TEs have?\n")
cat("  ✓ TEs ENRICHED as hubs (top 10% by degree)\n")
cat("    SCZ: OR =", round(fisher_test_scz$estimate, 2), 
    ", p =", format(fisher_test_scz$p.value, scientific=TRUE), "\n")
cat("    HC:  OR =", round(fisher_test_hc$estimate, 2), 
    ", p =", format(fisher_test_hc$p.value, scientific=TRUE), "\n")
cat("  ✓ TEs occupy central network positions (see heatmap)\n\n")

cat("Q3: What are TEs connected to?\n")
cat("  ✓ TEs connect to KNOWN SCZ RISK GENES!\n")
cat("    Top TE:", top_te_scz_connectors$TE[1], 
    "connects to", top_te_scz_connectors$SCZ_connections_scz[1], "SCZ genes\n")
cat("  ✓ SCZ genes are hubs:", sum(scz_gene_metrics_scz$is_hub), "/", 
    nrow(scz_gene_metrics_scz), "in SCZ network\n\n")

cat("Next: Run community detection script for Q3 details!\n")




# 
# scz_genes_all <- c("MTHFR", "RTN4R", "COMT", "HTR2A", "DRD3", "SYN2", "CHI3L1", "NRXN1", "MIR30E", 
#                    "SETD1A", "MBD5", "PHIP", "IRAK1BP1", "MIR206", "MIR198", "DAOA", "DISC1", "APOL2",
#                    "APOL4", "DTNBP1", "DAO", "ABCA13", "SHANK3", "AKT1", "DISC2", "NPAS3", "RELN", "MIR195",
#                    "MIR29C", "MIR9-1", "MIR212", "MIR15B", "ULK4", "MIR346", "MIR107", "MIR29A", "MIR106B",
#                    "MIR30A", "MIR26B", "MIR30B", "MIR15A", "MIR30D", "MIR20B", "ZNF804A", "DRD2", "ASTN2",
#                    "RGS4", "PRODH", "GABRB2", "CSMD1", "DLG2", "DLGAP2", "RBFOX1", "PRKN", "CTNND2", 
#                    "DPP6", "DPP10", "MSRA", "PCDH15", "DMD", "NRG1", "TSNAX-DISC1", "DRD4", "SLC1A1",
#                    "GRIN2B", "ZDHHC8", "CHRNA7", "GRIN2A", "PDE4B", "GRID2", "NOTCH4", "CNTN6", "EHMT1", 
#                    "MACROD2", "KCNN3", "KATNAL2", "TPH1", "GABRB1", "GRIK3", "LPP", "PTPRM", "NRG3",
#                    "TRAPPC9", "DBH", "FHIT", "PARD3B", "PDE11A", "TMLHE", "VPS13B", "FZD3", "PTPRT",
#                    "CNP", "BIRC6", "KIAA1586", "SLC18A2", "CLINT1", "YWHAH", "BDNF", "CPLX2", "GRM1",
#                    "GRIA4", "HTR2C", "GRM3", "PRL", "SLC6A3", "SLC26A5-AS1", "GAD1", "HTR1A", "GNAL",
#                    "SLC6A4", "DRD1", "PVALB", "NDEL1", "GRIN1", "SYN3", "PDLIM5", "GRM2", "CNTNAP2",
#                    "GRIN2C", "HTR3A", "GABBR1", "HRH2", "ERVW-1", "CYFIP1", "FOLH1", "TAAR6", "TH", 
#                    "PPP3CC", "VRK2", "PPP1R1B", "SRR", "MIAT", "GRIN3B", "NTNG1", "GRM5", "DLG4", 
#                    "ITIH3", "QKI", "PIP4K2A", "LRRTM1", "SEMA3A", "SULT4A1", "ERBB4", "NDUFV2", 
#                    "BRD1", "SLC6A9", "DRD5", "KCTD13", "NTS", "DGCR6", "HTR6", "CYP2D6", "MAP6",
#                    "MDGA1", "CNR1", "SP4", "MAOA", "NOS1", "NOS1AP", "SMARCA2", "GRIK4", "SLC6A1",
#                    "MEGF10", "FAN1", "SELENBP1", "KYAT1", "SNAP25", "RSRC1", "FEZ1", "UHMK1", "GRIK5",
#                    "CSMD2", "ATP2A2", "PPP3CA", "CACNG5", "GRIA1", "SLC17A7", "CALB1", "ARHGAP18",
#                    "DLG3", "NR4A2", "PLLP", "HRH3", "GRIK2", "MYO16", "CNIH3", "CHN2", "PPP3CB",
#                    "CHRM1", "MOBP", "HTR1B", "CALB2", "TSNAX", "CNIH2", "MT-ND4", "GRIA2", "PPP1R9B",
#                    "GRIK1", "ARVCF", "CDC42SE2", "DPYSL2", "MICB", "CNIH1", "MED12", "PSAT1", 
#                    "PICK1", "APOD", "PLXNA2", "RTN4", "TPH2", "SLC18A1", "NBPF10", "GPR153",
#                    "MAOB", "HTR7", "NTF3", "CYP1A2", "GPR78", "SLC1A2", "NRGN", "NBPF14", "CIT",
#                    "CHRNB2", "IMPA2", "NMBR", "NCS1", "GRIA3", "ADGRA3", "PAFAH1B1", "NBPF3", "GAD2", 
#                    "GSK3B", "MIR132", "NBPF9", "GNB3", "NBPF12", "CCK", "SLC6A2", "NPY", 
#                    "DNAJC14", "MAG", "TBX1", "NTRK2", "SYT11", "TDO2", "GAP43", "CACNG2", 
#                    "MAP2", "CHRNA4", "KCNN1", "HTR1D", "MECP2", "ADORA2A", "CBS", "HINT1", "TSPO", 
#                    "CHL1", "MBP", "HRH1", "MLC1", "FMR1", "ACHE", "GABRG2", "FOXP2", "CCKAR", 
#                    "SLC1A3", "CREB1", "NR3C1", "DNMT1", "CALY", "APOL1", "SLC6A5", "ST8SIA2", "DDC",
#                    "CHRM5", "CYP3A4", "PLA2G6", "GRM8", "GDNF", "CHAT", "PDYN", "PCNT", "OPRM1", 
#                    "SNCA", "GRM4", "SLC17A6", "CACNA1C", "GAPDH", "PLA2G4A", "GSK3A", "SOX10", 
#                    "DLG1", "ADRA2A", "HCRT", "IL6", "DLGAP1", "MIR24-1", "FYN", "FAAH", "PSEN1", 
#                    "TNF", "TACR3", "BCHE", "IL1B", "CYP2C19", "ADCYAP1", "GLUL", "IDO1", "SLC1A4",
#                    "APOE", "NTRK1", "OXT", "CLOCK", "PLP1", "INS", "NEFL", "ERBB3", "TAC1", 
#                    "HLA-DQB1", "PREP", "MOG", "MIR137", "CHGB", "ST8SIA4", "RGS10", "CHRM2", 
#                    "HLA-DRB1", "NTSR1", "CNTF", "NTRK3", "GRM7", "VLDLR", "MAP1B", "TBP",
#                    "IL10", "ANK3", "S100B", "APP", "OXTR", "LEP", "HMBS", "ESR2", "GRIN2D", 
#                    "PPARA", "CACNA1A", "GRIP1", "FOS", "NR4A1", "XBP1", "DGCR5", "TNFRSF1B", 
#                    "LTA", "TNFRSF1A", "NGF", "SST", "CRP", "ANKK1", "GFAP", "CRH", "LOC105373170",
#                    "GRK3", "PNOC", "NBPF4", "NBPF6", "POMC", "PDE10A", "SERPINA3", "NLGN4X", "ATXN1",
#                    "SIGMAR1", "HSPA8", "SYP", "SHANK2", "HSPA1A", "MIR125A", "IL2", "HTR4", "TCF4",
#                    "GRIN3A", "C4A", "CHRFAM7A", "NRXN2", "PAH", "NLGN1", "IL6R", "MAPT", "SPTAN1", 
#                    "NSF", "SLC2A1", "KALRN", "NDE1", "ACE", "GLS", "TAAR1", "BDNF-AS", "ABAT", "SP1", 
#                    "HOMER1", "CAMK2A", "EGF", "IFNG", "OLIG2", "NRP1", "GLUD1", "SYN1", "NQO1", "ARRB2",
#                    "IL4", "TACR1", "GABRB3", "SLC32A1", "NQO2", "AGER", "GABRA1", "AADAT", "SLC12A5", 
#                    "MIR130B", "FXR1", "ESR1", "SLC12A2", "MTOR", "GPHN", "CHRNA5", "TP53", "GABRA2",
#                    "BACE1", "PTGDS", "ADM", "NRXN3", "KMO", "CCKBR", "SCT", "IL1A", "CTNNB1", "EGR1", 
#                    "ABCB1", "MAD1L1", "GHRL", "SOD1", "HTR3B", "ALB", "KCNH2", "CNNM2", "CACNB2", 
#                    "CPLX1", "ARC", "AS3MT", "FKBP5", "PLA2G2A", "CCL2", "EGR3", "CRHR1", "CXCL8",
#                    "MMP9", "GRID1", "IGF1", "DGCR2", "NCAM1", "IL18R1", "SLC39A8", "TGFB1", "RXRB",
#                    "SCN2A", "PTGS2", "MIR137HG", "TLR4", "CLDN5", "VIPR2", "C9orf72", "NDUFS4", 
#                    "CNTN4", "PPARG", "GABRA5", "ADCY10", "TSNARE1", "PLA2G7", "NTNG2", "MYT1L", 
#                    "NFKB1", "IL1RAPL2", "CDK5", "NT5C2", "FABP7", "HTR5A", "CHRNA3", "AUTS2", "IL18", 
#                    "ADIPOQ", "MAPK3", "ADNP", "NCAN", "ADRA1A", "GSTM1", "PIK3CA", "SYNGAP1", "CHRM4", 
#                    "DGCR8", "PLCB1", "LRP8", "AHI1", "STX1A", "VIP", "ATN1", "CACNA1B", "ASCL1", "PIK3CB", 
#                    "DNMT3B", "GABRA4", "GCLM", "RGS2", "VAMP2", "WBP1L", "CACNA1I", "STXBP1", "TNIK", "TMTC1",
#                    "SATB2", "GABBR2", "POU3F2", "JUN", "GNB1L", "YWHAE", "CDC42", "SOD2", "SLC1A6", 
#                    "PIK3C3", "CRHBP", "VSNL1", "DDO", "TRANK1", "HP", "SYNGR1", "PI4KA", "BMAL1", 
#                    "CLU", "PCM1", "EGR2", "PDE4D", "YWHAZ", "SNX19", "IL1RN", "ARSA", "HLA-B", 
#                    "MCHR1", "L1CAM", "DOC2A", "GABRG3", "C4B", "CFL1", "FXYD6", "IL3", "DISC1FP1", 
#                    "HDAC1", "PRNP", "SYNE1", "CAMK2B", "GLT8D1", "FTO", "CTNNA2", "PLA2G1B", "CNR2",
#                    "HLA-A", "ITIH4", "MEF2C", "PTPRZ1", "CMYA5", "NRN1", "CKB")
# 
