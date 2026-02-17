# ============================================================
# PANEL C: TE Rewiring to SCZ Risk Genes
# ============================================================

library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(ggplot2)

tes_oi <- c("L1FLnI_2p13.1f", "L1FLnI_8q13.3e", "HERVL_4q32.1a")
scz_in_network <- intersect(scz_risk_genes, V(hc_graph)$name)

# ---- Build edge list ----
edges_list <- list()

for(te in tes_oi) {
  hc_neighbors <- intersect(neighbors(hc_graph, te)$name, scz_in_network)
  scz_neighbors <- intersect(neighbors(scz_graph, te)$name, scz_in_network)
  
  all_neighbors <- union(hc_neighbors, scz_neighbors)
  
  for(gene in all_neighbors) {
    in_hc <- gene %in% hc_neighbors
    in_scz <- gene %in% scz_neighbors
    
    edge_type <- ifelse(in_hc & in_scz, "Conserved",
                        ifelse(in_hc & !in_scz, "Lost in SCZ",
                               "Gained in SCZ"))
    
    edges_list[[length(edges_list) + 1]] <- data.frame(
      from = te, to = gene, edge_type = edge_type, stringsAsFactors = FALSE
    )
  }
}

edges_df <- do.call(rbind, edges_list)
cat("Edge types:\n")
print(table(edges_df$edge_type))

# ---- Build graph ----
all_nodes <- unique(c(edges_df$from, edges_df$to))

nodes_df <- data.frame(name = all_nodes, stringsAsFactors = FALSE) %>%
  mutate(
    node_type = ifelse(name %in% tes_oi, "TE", "SCZ Risk Gene"),
    color_group = case_when(
      grepl("^L1", name) & node_type == "TE" ~ "LINE",
      node_type == "TE" ~ "LTR",
      TRUE ~ "SCZ Risk Gene"
    )
  )

g_rewire <- graph_from_data_frame(edges_df, directed = FALSE, vertices = nodes_df)

# ---- Custom bipartite layout ----
te_names <- nodes_df$name[nodes_df$node_type == "TE"]
gene_names <- sort(nodes_df$name[nodes_df$node_type == "SCZ Risk Gene"])

n_tes <- length(te_names)
n_genes <- length(gene_names)

layout_mat <- matrix(0, nrow = nrow(nodes_df), ncol = 2)
rownames(layout_mat) <- nodes_df$name

# TEs on the left
te_y_positions <- seq(from = n_genes * 0.35, to = -n_genes * 0.35, length.out = n_tes)
for(i in seq_along(te_names)) {
  layout_mat[te_names[i], ] <- c(-1.5, te_y_positions[i])
}

# Genes on the right
gene_y_positions <- seq(from = n_genes * 0.5, to = -n_genes * 0.5, length.out = n_genes)
for(i in seq_along(gene_names)) {
  layout_mat[gene_names[i], ] <- c(1.5, gene_y_positions[i])
}

# Reorder to match vertex order
layout_ordered <- layout_mat[V(g_rewire)$name, ]

# ---- Edge colors ----
edge_colors <- c(
  "Conserved" = "#777777",
  "Lost in SCZ" = "#4DAF4A",
  "Gained in SCZ" = "#E41A1C"
)

# ---- Node aesthetics ----
node_colors <- c("LINE" = "#FF8C00", "LTR" = "#CD950C", "SCZ Risk Gene" = "#63B8FF")
node_shapes <- c("LINE" = 23, "LTR" = 23, "SCZ Risk Gene" = 21)
node_sizes <- c("LINE" = 7, "LTR" = 7, "SCZ Risk Gene" = 5)

# ---- Determine label hjust based on position ----
label_hjust <- ifelse(layout_ordered[,1] < 0, 1.2, -0.2)

# ---- Build plot ----
tg_rewire <- as_tbl_graph(g_rewire)

p_rewire <- ggraph(tg_rewire, layout = "manual", 
                   x = layout_ordered[,1], y = layout_ordered[,2]) +
  # Edges with curvature for visual interest
  geom_edge_link(
    aes(color = edge_type),
    alpha = 0.6,
    width = 0.9,
    show.legend = TRUE
  ) +
  scale_edge_color_manual(
    values = edge_colors,
    name = "Connection Status",
    labels = c("Conserved" = "Present in both",
               "Gained in SCZ" = "Gained in schizophrenia", 
               "Lost in SCZ" = "Lost in schizophrenia")
  ) +
  # Nodes
  geom_node_point(
    aes(fill = color_group, shape = color_group, size = color_group),
    color = "gray30", stroke = 0.8
  ) +
  scale_fill_manual(values = node_colors, name = "Node Type") +
  scale_shape_manual(values = node_shapes, name = "Node Type") +
  scale_size_manual(values = node_sizes, name = "Node Type") +
  # Labels
  geom_node_text(
    aes(label = name),
    hjust = label_hjust,
    size = 3.5,
    fontface = ifelse(V(g_rewire)$node_type == "TE", "bold", "plain")
  ) +
  # Annotations
  annotate("text", x = -1.5, y = max(te_y_positions) + 1.5, 
           label = "Transposable\nElements", fontface = "bold", size = 4.5, color = "gray30") +
  annotate("text", x = 1.5, y = max(gene_y_positions) + 1.5, 
           label = "SCZ Risk\nGenes", fontface = "bold", size = 4.5, color = "gray30") +
  # Theme
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 10),
    plot.margin = margin(10, 50, 10, 50)
  ) +
  labs(
    title = "TE Regulatory Rewiring to Schizophrenia Risk Genes",
    subtitle = "Connections absent in healthy controls emerge in schizophrenia network"
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 5)),
    shape = guide_legend(order = 1),
    size = guide_legend(order = 1),
    edge_color = guide_legend(order = 2, override.aes = list(edge_width = 2, edge_alpha = 1))
  )

print(p_rewire)

# ---- Save ----
ggsave(paste0(output_dir, "panel_C_te_rewiring.pdf"),
       p_rewire, width = 10, height = 8, dpi = 300)
ggsave(paste0(output_dir, "panel_C_te_rewiring.png"),
       p_rewire, width = 10, height = 8, dpi = 300)

cat("Panel C saved!\n")
