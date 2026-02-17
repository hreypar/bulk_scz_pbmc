# ============================================================
# PANELS A & B: Side-by-side network visualization in R
# ============================================================

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(patchwork)

# ---- Redefine node selection ----
tes_oi <- c("L1FLnI_2p13.1f", "L1FLnI_8q13.3e", "HERVL_4q32.1a")
scz_in_network <- intersect(scz_risk_genes, V(hc_graph)$name)

te_nodes <- V(hc_graph)$name[V(hc_graph)$type == "TE"]
te_deg_hc <- sort(degree(hc_graph, v = te_nodes), decreasing = TRUE)
top_tes <- setdiff(names(te_deg_hc)[1:20], tes_oi)

gene_nodes <- V(hc_graph)$name[V(hc_graph)$type == "gene"]
gene_deg_hc <- sort(degree(hc_graph, v = gene_nodes), decreasing = TRUE)
top_genes <- setdiff(names(gene_deg_hc)[1:20], scz_in_network)

selected_nodes <- unique(c(tes_oi, scz_in_network, top_tes[1:15], top_genes[1:15]))

# ---- Extract subgraphs ----
hc_sub <- induced_subgraph(hc_graph, vids = selected_nodes)
scz_sub <- induced_subgraph(scz_graph, vids = selected_nodes)

# ---- Assign display attributes ----
assign_display_attrs <- function(g, full_graph) {
  node_names <- V(g)$name
  node_type <- V(g)$type
  
  # Category
  category <- rep("Gene Hub", length(node_names))
  category[node_names %in% scz_in_network] <- "SCZ Risk Gene"
  category[node_type == "TE"] <- "TE Hub"
  category[node_names %in% tes_oi] <- "TE of Interest"
  V(g)$category <- category
  
  # Color group (for node coloring by molecular type)
  color_group <- rep("Gene", length(node_names))
  for(i in seq_along(node_names)) {
    if(node_type[i] == "TE") {
      if(grepl("^L1", node_names[i])) {
        color_group[i] <- "LINE"
      } else {
        color_group[i] <- "LTR"
      }
    }
  }
  V(g)$color_group <- color_group
  
  # Size by degree in full network
  full_deg <- degree(full_graph, v = node_names)
  V(g)$degree_full <- full_deg
  
  # Labels: show all but vary size
  V(g)$label <- node_names
  
  return(g)
}

hc_sub <- assign_display_attrs(hc_sub, hc_graph)
scz_sub <- assign_display_attrs(scz_sub, scz_graph)

# ---- Compute shared layout on HC ----
set.seed(42)
shared_layout <- create_layout(as_tbl_graph(hc_sub), layout = "fr", niter = 5000)

# Extract x, y coordinates
layout_coords <- data.frame(
  name = shared_layout$name,
  x = shared_layout$x,
  y = shared_layout$y
)

# ---- Build HC plot ----
hc_tbl <- as_tbl_graph(hc_sub) %>%
  activate(nodes) %>%
  left_join(layout_coords, by = "name")

# Create manual layout for HC
hc_layout <- create_layout(hc_tbl, layout = "manual", x = layout_coords$x, y = layout_coords$y)

# Node size scaling
size_range <- range(hc_layout$degree_full)
hc_layout$node_size <- 2 + (hc_layout$degree_full - size_range[1]) / 
  (size_range[2] - size_range[1]) * 8

# Node colors
node_fill_colors <- c("Gene" = "#63B8FF", "LINE" = "#FF8C00", "LTR" = "#CD950C")
node_shape_values <- c("Gene" = 21, "LINE" = 23, "LTR" = 23)  # 21=circle, 23=diamond

# Label size: bigger for TEs and SCZ genes
hc_layout$label_size <- ifelse(hc_layout$category %in% c("TE of Interest", "TE Hub"), 2.8,
                               ifelse(hc_layout$category == "SCZ Risk Gene", 2.2, 1.8))

# Border color for emphasis
hc_layout$border_col <- ifelse(hc_layout$category == "TE of Interest", "#E41A1C",
                               ifelse(hc_layout$category == "SCZ Risk Gene", "#333333",
                                      ifelse(hc_layout$category == "TE Hub", "#FF8C00", "#63B8FF")))

p_hc <- ggraph(hc_layout) +
  geom_edge_link(alpha = 0.08, color = "gray50", width = 0.2) +
  geom_node_point(aes(fill = color_group, shape = color_group, size = node_size),
                  color = hc_layout$border_col, stroke = 0.6) +
  geom_node_text(aes(label = label), size = hc_layout$label_size, 
                 repel = TRUE, max.overlaps = 20, 
                 segment.color = "gray70", segment.size = 0.2,
                 fontface = ifelse(hc_layout$color_group != "Gene", "bold", "plain")) +
  scale_fill_manual(values = node_fill_colors, name = "Node Type") +
  scale_shape_manual(values = node_shape_values, name = "Node Type") +
  scale_size_identity() +
  labs(title = "Healthy Control Network",
       subtitle = paste0(vcount(hc_sub), " nodes, ", ecount(hc_sub), " edges")) +
  theme_void() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "#4DAF4A"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

# ---- Build SCZ plot with SAME layout ----
scz_tbl <- as_tbl_graph(scz_sub) %>%
  activate(nodes) %>%
  left_join(layout_coords, by = "name")

scz_layout <- create_layout(scz_tbl, layout = "manual", x = layout_coords$x, y = layout_coords$y)

# Use same size scaling as HC (important for comparison)
scz_layout$node_size <- 2 + (scz_layout$degree_full - size_range[1]) / 
  (size_range[2] - size_range[1]) * 8

scz_layout$label_size <- ifelse(scz_layout$category %in% c("TE of Interest", "TE Hub"), 2.8,
                                ifelse(scz_layout$category == "SCZ Risk Gene", 2.2, 1.8))

scz_layout$border_col <- ifelse(scz_layout$category == "TE of Interest", "#E41A1C",
                                ifelse(scz_layout$category == "SCZ Risk Gene", "#333333",
                                       ifelse(scz_layout$category == "TE Hub", "#FF8C00", "#63B8FF")))

p_scz <- ggraph(scz_layout) +
  geom_edge_link(alpha = 0.08, color = "gray50", width = 0.2) +
  geom_node_point(aes(fill = color_group, shape = color_group, size = node_size),
                  color = scz_layout$border_col, stroke = 0.6) +
  geom_node_text(aes(label = label), size = scz_layout$label_size, 
                 repel = TRUE, max.overlaps = 20, 
                 segment.color = "gray70", segment.size = 0.2,
                 fontface = ifelse(scz_layout$color_group != "Gene", "bold", "plain")) +
  scale_fill_manual(values = node_fill_colors, name = "Node Type") +
  scale_shape_manual(values = node_shape_values, name = "Node Type") +
  scale_size_identity() +
  labs(title = "Schizophrenia Network",
       subtitle = paste0(vcount(scz_sub), " nodes, ", ecount(scz_sub), " edges (−75.7%)")) +
  theme_void() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "#E41A1C"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

# ---- Combine side by side ----
p_combined <- p_hc + p_scz + 
  plot_layout(ncol = 2) +
  plot_annotation(
    caption = "Node color: Gene (blue) | LINE TE (orange) | LTR TE (gold)\nNode shape: Circle (gene) | Diamond (TE) | Size proportional to full network degree",
    theme = theme(
      plot.caption = element_text(size = 9, hjust = 0.5, color = "gray40")
    )
  )

# ---- Save ----
output_dir <- "results/network_analysis/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(paste0(output_dir, "panels_AB_network_comparison.pdf"),
       p_combined, width = 20, height = 10, dpi = 300)
ggsave(paste0(output_dir, "panels_AB_network_comparison.png"),
       p_combined, width = 20, height = 10, dpi = 300)

cat("Panels A & B saved!\n")

# ---- Also save individually for flexible poster layout ----
ggsave(paste0(output_dir, "panel_A_hc_network.pdf"),
       p_hc, width = 10, height = 10, dpi = 300)
ggsave(paste0(output_dir, "panel_B_scz_network.pdf"),
       p_scz, width = 10, height = 10, dpi = 300)
