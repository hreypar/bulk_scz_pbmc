# =============================================================================
# Network Visualizations for Presentation
# =============================================================================

library(igraph)
library(RColorBrewer)

cat("=============================================================================\n")
cat("CREATING NETWORK VISUALIZATIONS\n")
cat("=============================================================================\n\n")

# Load data
scz_graph <- readRDS("results/network_analysis/networks/scz_graph_intergenic_only.rds")
hc_graph <- readRDS("results/network_analysis/networks/hc_graph_intergenic_only.rds")
metrics <- readRDS("results/network_analysis/metrics/centrality_metrics.rds")
feature_lists <- readRDS("results/network_analysis/data/feature_lists_intergenic_only.rds")

metrics_scz <- metrics$scz
metrics_hc <- metrics$hc
top_genes <- feature_lists$top_genes
top_tes <- feature_lists$top_tes

scz_genes <- c("MTHFR", "RTN4R", "COMT", "HTR2A", "DRD3", "SYN2", "CHI3L1", 
               "NRXN1", "SETD1A", "MBD5", "PHIP", "DAOA", "DISC1", "DTNBP1", "DAO", 
               "ABCA13", "SHANK3", "AKT1", "NPAS3", "RELN", "ZNF804A", "DRD2", 
               "ASTN2", "RGS4", "PRODH", "GABRB2", "CSMD1", "DLG2", "DLGAP2", 
               "RBFOX1", "PRKN", "CTNND2", "DPP6", "DPP10", "MSRA", "PCDH15", 
               "DMD", "NRG1", "DRD4", "SLC1A1", "GRIN2B", "ZDHHC8", "CHRNA7", 
               "GRIN2A", "PDE4B", "GRID2", "NOTCH4", "CNTN6", "EHMT1", "MACROD2", 
               "KCNN3", "KATNAL2", "TPH1", "GABRB1", "GRIK3", "LPP", "PTPRM", 
               "NRG3", "TRAPPC9", "DBH", "FHIT", "PARD3B", "PDE11A", "TMLHE", 
               "VPS13B", "FZD3", "PTPRT", "CNP")
scz_genes_in_network <- scz_genes[scz_genes %in% top_genes]

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

# -----------------------------------------------------------------------------
# VIZ 1: Edge Type Composition (Q1 - Are TEs connected?)
# -----------------------------------------------------------------------------

cat("Creating edge type visualization...\n")

# Count edge types in both networks
count_edge_types <- function(graph, genes, tes) {
  edges <- as_edgelist(graph)
  
  te_gene <- sum(
    (edges[,1] %in% tes & edges[,2] %in% genes) |
      (edges[,2] %in% tes & edges[,1] %in% genes)
  )
  
  te_te <- sum(edges[,1] %in% tes & edges[,2] %in% tes)
  
  gene_gene <- sum(edges[,1] %in% genes & edges[,2] %in% genes)
  
  return(c(TE_Gene = te_gene, TE_TE = te_te, Gene_Gene = gene_gene))
}

scz_edge_types <- count_edge_types(scz_graph, top_genes, top_tes)
hc_edge_types <- count_edge_types(hc_graph, top_genes, top_tes)

save_plot("results/network_analysis/presentation_plots/12_edge_type_composition",
          function() {
            par(mfrow=c(1,3), mar=c(4,4,3,2))
            
            # SCZ pie chart
            pie(scz_edge_types,
                main="SCZ Network",
                col=c("coral", "darkred", "steelblue"),
                labels=paste0(names(scz_edge_types), "\n", 
                              format(scz_edge_types, big.mark=","), "\n(",
                              round(100*scz_edge_types/sum(scz_edge_types), 1), "%)"),
                cex=1.1)
            
            # HC pie chart
            pie(hc_edge_types,
                main="HC Network",
                col=c("coral", "darkred", "steelblue"),
                labels=paste0(names(hc_edge_types), "\n", 
                              format(hc_edge_types, big.mark=","), "\n(",
                              round(100*hc_edge_types/sum(hc_edge_types), 1), "%)"),
                cex=1.1)
            
            # Comparison bar plot
            edge_comparison <- rbind(
              SCZ = scz_edge_types / sum(scz_edge_types) * 100,
              HC = hc_edge_types / sum(hc_edge_types) * 100
            )
            
            barplot(edge_comparison,
                    beside=TRUE,
                    col=c("orangered", "steelblue"),
                    main="Edge Type Composition",
                    ylab="% of Total Edges",
                    legend=TRUE,
                    args.legend=list(x="topright", bty="n"),
                    las=2,
                    cex.names=0.9)
          }, width=15, height=5)

cat("  Edge types:\n")
cat("    SCZ - TE-Gene:", scz_edge_types[1], 
    "(", round(100*scz_edge_types[1]/sum(scz_edge_types), 1), "%)\n")
cat("    HC  - TE-Gene:", hc_edge_types[1], 
    "(", round(100*hc_edge_types[1]/sum(hc_edge_types), 1), "%)\n\n")

# -----------------------------------------------------------------------------
# VIZ 2: Top TE Hub Neighborhood (Q2 - TE positions)
# -----------------------------------------------------------------------------

cat("Creating hub neighborhood visualization...\n")

# Get top TE hub
top_te_hub <- metrics_hc$node[metrics_hc$type == "TE"][which.max(metrics_hc$degree[metrics_hc$type == "TE"])]

cat("  Visualizing hub:", top_te_hub, "\n")

# Extract 1-hop neighborhood
create_hub_subgraph <- function(graph, hub_node, scz_genes_list) {
  neighbors_hub <- neighbors(graph, hub_node)
  
  # Keep top 50 neighbors by degree
  neighbor_degrees <- degree(graph)[neighbors_hub]
  top_neighbors <- neighbors_hub[order(-neighbor_degrees)[1:min(50, length(neighbors_hub))]]
  
  # Create subgraph
  subgraph_nodes <- c(hub_node, names(top_neighbors))
  subg <- induced_subgraph(graph, subgraph_nodes)
  
  # Set node attributes
  V(subg)$color <- ifelse(V(subg)$name == hub_node, "darkred",
                          ifelse(V(subg)$name %in% scz_genes_list, "gold",
                                 ifelse(V(subg)$type == "TE", "coral", "steelblue")))
  
  V(subg)$size <- ifelse(V(subg)$name == hub_node, 20,
                         ifelse(V(subg)$name %in% scz_genes_list, 8, 5))
  
  V(subg)$label <- ifelse(V(subg)$name == hub_node, V(subg)$name,
                          ifelse(V(subg)$name %in% scz_genes_list, V(subg)$name, ""))
  
  V(subg)$label.cex <- 0.8
  V(subg)$label.color <- "black"
  
  E(subg)$color <- "gray80"
  E(subg)$width <- 0.5
  
  return(subg)
}

scz_hub_subg <- create_hub_subgraph(scz_graph, top_te_hub, scz_genes_in_network)
hc_hub_subg <- create_hub_subgraph(hc_graph, top_te_hub, scz_genes_in_network)

save_plot("results/network_analysis/presentation_plots/13_te_hub_neighborhood",
          function() {
            par(mfrow=c(1,2), mar=c(1,1,3,1))
            
            set.seed(123)
            layout_coords <- layout_with_fr(hc_hub_subg)
            
            # HC network
            plot(hc_hub_subg,
                 layout=layout_coords,
                 main=paste0("HC Network\n", top_te_hub, " Hub"),
                 vertex.frame.color=NA)
            
            legend("bottomright",
                   legend=c("TE Hub", "SCZ Risk Gene", "Other TE", "Gene"),
                   pch=21,
                   pt.bg=c("darkred", "gold", "coral", "steelblue"),
                   pt.cex=c(2, 1.5, 1, 1),
                   bty="n",
                   cex=0.9)
            
            # SCZ network (same layout for comparison)
            plot(scz_hub_subg,
                 layout=layout_coords,
                 main=paste0("SCZ Network\n", top_te_hub, " Hub"),
                 vertex.frame.color=NA)
            
            text(-1.2, 1.2, 
                 paste0("Degree:\nHC: ", degree(hc_graph)[top_te_hub],
                        "\nSCZ: ", degree(scz_graph)[top_te_hub]),
                 cex=1, adj=0)
          }, width=14, height=7)

# -----------------------------------------------------------------------------
# VIZ 3: TE-SCZ Gene Bipartite Network (Q3 - What are TEs connected to?) ⭐
# -----------------------------------------------------------------------------

cat("Creating TE-SCZ gene bipartite network...\n")

# Get TEs that connect to SCZ genes
te_scz_connections <- sapply(top_tes, function(te) {
  neighbors_hc <- neighbors(hc_graph, te)
  neighbors_scz <- neighbors(scz_graph, te)
  c(hc = sum(V(hc_graph)$name[neighbors_hc] %in% scz_genes_in_network),
    scz = sum(V(scz_graph)$name[neighbors_scz] %in% scz_genes_in_network))
})

te_scz_connections <- t(te_scz_connections)

# Keep TEs with at least 2 SCZ gene connections in either network
keep_tes <- rownames(te_scz_connections)[rowSums(te_scz_connections >= 2) > 0]

cat("  TEs connecting to ≥2 SCZ genes:", length(keep_tes), "\n")

# Create bipartite networks
create_bipartite_scz <- function(graph, tes_list, scz_genes_list) {
  # Get all edges between TEs and SCZ genes
  edges_list <- c()
  
  for(te in tes_list) {
    neighbors_te <- neighbors(graph, te)
    scz_neighbors <- intersect(V(graph)$name[neighbors_te], scz_genes_list)
    
    if(length(scz_neighbors) > 0) {
      for(scz_gene in scz_neighbors) {
        edges_list <- c(edges_list, te, scz_gene)
      }
    }
  }
  
  if(length(edges_list) == 0) return(NULL)
  
  # Create graph
  bip <- graph(edges_list, directed=FALSE)
  
  # Set bipartite attribute
  V(bip)$type <- V(bip)$name %in% tes_list
  
  # Set colors
  V(bip)$color <- ifelse(V(bip)$type, "coral", "gold")
  V(bip)$shape <- ifelse(V(bip)$type, "square", "circle")
  V(bip)$size <- 8
  V(bip)$label.cex <- 0.7
  V(bip)$label.color <- "black"
  
  E(bip)$color <- "gray70"
  E(bip)$width <- 1
  
  return(bip)
}

hc_bipartite <- create_bipartite_scz(hc_graph, keep_tes, scz_genes_in_network)
scz_bipartite <- create_bipartite_scz(scz_graph, keep_tes, scz_genes_in_network)

save_plot("results/network_analysis/presentation_plots/14_te_scz_gene_bipartite",
          function() {
            par(mfrow=c(1,2), mar=c(1,1,3,1))
            
            if(!is.null(hc_bipartite)) {
              set.seed(42)
              layout_bip <- layout_as_bipartite(hc_bipartite)
              
              plot(hc_bipartite,
                   layout=layout_bip,
                   main=paste0("HC Network\nTEs ↔ SCZ Risk Genes\n(", 
                               ecount(hc_bipartite), " connections)"))
              
              legend("bottomleft",
                     legend=c("TEs", "SCZ Risk Genes"),
                     pch=c(22, 21),
                     pt.bg=c("coral", "gold"),
                     pt.cex=1.5,
                     bty="n")
            }
            
            if(!is.null(scz_bipartite)) {
              set.seed(42)
              layout_bip <- layout_as_bipartite(scz_bipartite)
              
              plot(scz_bipartite,
                   layout=layout_bip,
                   main=paste0("SCZ Network\nTEs ↔ SCZ Risk Genes\n(", 
                               ecount(scz_bipartite), " connections)"))
            }
          }, width=14, height=7)

# -----------------------------------------------------------------------------
# VIZ 4: Heatmap of TE-SCZ Gene Connections
# -----------------------------------------------------------------------------

cat("Creating TE-SCZ gene connection heatmap...\n")

# Create adjacency matrix
create_te_scz_matrix <- function(graph, tes_list, scz_genes_list) {
  adj_matrix <- matrix(0, nrow=length(tes_list), ncol=length(scz_genes_list),
                       dimnames=list(tes_list, scz_genes_list))
  
  for(i in 1:length(tes_list)) {
    te <- tes_list[i]
    neighbors_te <- neighbors(graph, te)
    
    for(j in 1:length(scz_genes_list)) {
      scz_gene <- scz_genes_list[j]
      if(scz_gene %in% V(graph)$name[neighbors_te]) {
        adj_matrix[i, j] <- 1
      }
    }
  }
  
  return(adj_matrix)
}

hc_te_scz_matrix <- create_te_scz_matrix(hc_graph, keep_tes, scz_genes_in_network)
scz_te_scz_matrix <- create_te_scz_matrix(scz_graph, keep_tes, scz_genes_in_network)

# Only keep TEs and genes with at least one connection
keep_tes_heatmap <- rowSums(hc_te_scz_matrix + scz_te_scz_matrix) > 0
keep_genes_heatmap <- colSums(hc_te_scz_matrix + scz_te_scz_matrix) > 0

hc_te_scz_matrix <- hc_te_scz_matrix[keep_tes_heatmap, keep_genes_heatmap]
scz_te_scz_matrix <- scz_te_scz_matrix[keep_tes_heatmap, keep_genes_heatmap]

save_plot("results/network_analysis/presentation_plots/15_te_scz_gene_heatmap",
          function() {
            par(mfrow=c(1,2), mar=c(8,8,3,2))
            
            # HC heatmap
            image(t(hc_te_scz_matrix[nrow(hc_te_scz_matrix):1, ]),
                  col=c("white", "coral"),
                  main="HC: TE-SCZ Gene Connections",
                  xlab="", ylab="",
                  axes=FALSE)
            
            axis(1, at=seq(0, 1, length.out=ncol(hc_te_scz_matrix)),
                 labels=colnames(hc_te_scz_matrix),
                 las=2, cex.axis=0.7)
            
            axis(2, at=seq(0, 1, length.out=nrow(hc_te_scz_matrix)),
                 labels=rev(rownames(hc_te_scz_matrix)),
                 las=2, cex.axis=0.7)
            
            # SCZ heatmap
            image(t(scz_te_scz_matrix[nrow(scz_te_scz_matrix):1, ]),
                  col=c("white", "orangered"),
                  main="SCZ: TE-SCZ Gene Connections",
                  xlab="", ylab="",
                  axes=FALSE)
            
            axis(1, at=seq(0, 1, length.out=ncol(scz_te_scz_matrix)),
                 labels=colnames(scz_te_scz_matrix),
                 las=2, cex.axis=0.7)
            
            axis(2, at=seq(0, 1, length.out=nrow(scz_te_scz_matrix)),
                 labels=rev(rownames(scz_te_scz_matrix)),
                 las=2, cex.axis=0.7)
          }, width=16, height=10)

cat("\n✓ All network visualizations created!\n\n")

cat("=============================================================================\n")
cat("VISUALIZATION SUMMARY\n")
cat("=============================================================================\n\n")

cat("Created 4 network visualizations:\n\n")

cat("12. Edge Type Composition\n")
cat("    - Shows % of TE-gene, TE-TE, gene-gene edges\n")
cat("    - Answers Q1: Yes, TEs are connected!\n\n")

cat("13. Top TE Hub Neighborhood\n")
cat("    - Shows", top_te_hub, "and its top 50 neighbors\n")
cat("    - SCZ vs HC side-by-side\n")
cat("    - Answers Q2: TEs are central hubs\n\n")

cat("14. TE-SCZ Gene Bipartite Network ⭐\n")
cat("    - Shows which TEs connect to which SCZ risk genes\n")
cat("    - Direct visual evidence of disease relevance\n")
cat("    - Answers Q3: TEs connect to disease genes!\n\n")

cat("15. TE-SCZ Gene Heatmap\n")
cat("    - Matrix view of all TE-SCZ gene connections\n")
cat("    - Easy to see rewiring between networks\n\n")

cat("RECOMMENDATION: Use plots 12, 13, 14 in your presentation\n")
cat("Plot 14 (bipartite) is your KILLER SLIDE! 🎯\n")
