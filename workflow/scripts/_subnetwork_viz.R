# ============================================================
# PANELS A & B: BIG CANVAS, ALL LABELS READABLE
# ============================================================

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(patchwork)

hc_graph <- readRDS("results/network_analysis/networks/hc_graph_intergenic_only.rds")
scz_graph <- readRDS("results/network_analysis/networks/scz_graph_intergenic_only.rds")

scz_risk_genes <- c("MTHFR", "RTN4R", "COMT", "HTR2A", "DRD3", "SYN2", "CHI3L1", "NRXN1", "SETD1A", "MBD5", "PHIP", "DAOA", "DISC1", "DTNBP1", "DAO", "ABCA13", "SHANK3", "AKT1", "NPAS3", "RELN", "ZNF804A", "DRD2", "ASTN2", "RGS4", "PRODH", "GABRB2", "CSMD1", "DLG2", "DLGAP2", "RBFOX1", "PRKN", "CTNND2", "DPP6", "DPP10", "MSRA", "PCDH15", "DMD", "NRG1", "DRD4", "SLC1A1", "GRIN2B", "ZDHHC8", "CHRNA7", "GRIN2A", "PDE4B", "GRID2", "NOTCH4", "CNTN6", "EHMT1", "MACROD2", "KCNN3", "KATNAL2", "TPH1", "GABRB1", "GRIK3", "LPP", "PTPRM", "NRG3", "TRAPPC9", "DBH", "FHIT", "PARD3B", "PDE11A", "TMLHE", "VPS13B", "FZD3", "PTPRT", "CNP")
tes_oi <- c("L1FLnI_2p13.1f", "L1FLnI_8q13.3e", "HERVL_4q32.1a")

# ---- Node selection ----
scz_in_network <- intersect(scz_risk_genes, V(hc_graph)$name)

te_nodes <- V(hc_graph)$name[V(hc_graph)$type == "TE"]
te_deg_hc <- sort(degree(hc_graph, v = te_nodes), decreasing = TRUE)
top_tes <- setdiff(names(te_deg_hc)[1:20], tes_oi)

gene_nodes <- V(hc_graph)$name[V(hc_graph)$type == "gene"]
gene_deg_hc <- sort(degree(hc_graph, v = gene_nodes), decreasing = TRUE)
top_genes <- setdiff(names(gene_deg_hc)[1:20], scz_in_network)

selected_nodes <- unique(c(tes_oi, scz_in_network, top_tes[1:15], top_genes[1:15]))

hc_sub <- induced_subgraph(hc_graph, vids = selected_nodes)
scz_sub <- induced_subgraph(scz_graph, vids = selected_nodes)

# ---- Assign attributes ----
assign_display_attrs <- function(g, full_graph) {
  node_names <- V(g)$name
  node_type <- V(g)$type
  
  category <- rep("Gene Hub", length(node_names))
  category[node_names %in% scz_in_network] <- "SCZ Gene"
  category[node_type == "TE"] <- "TE Hub"
  category[node_names %in% tes_oi] <- "TE of Interest"
  V(g)$category <- category
  
  color_group <- rep("Gene", length(node_names))
  for(i in seq_along(node_names)) {
    if(node_type[i] == "TE") {
      color_group[i] <- ifelse(grepl("^L1", node_names[i]), "LINE", "LTR")
    }
  }
  V(g)$color_group <- color_group
  V(g)$degree_full <- degree(full_graph, v = node_names)
  V(g)$label <- node_names
  return(g)
}

hc_sub <- assign_display_attrs(hc_sub, hc_graph)
scz_sub <- assign_display_attrs(scz_sub, scz_graph)

# ---- Shared layout ----
set.seed(42)
shared_layout <- create_layout(as_tbl_graph(hc_sub), layout = "fr", niter = 5000)
layout_coords <- data.frame(
  name = shared_layout$name,
  x = shared_layout$x,
  y = shared_layout$y
)

size_range <- range(V(hc_sub)$degree_full)

# ---- Panel builder ----
build_panel <- function(g, layout_coords, size_range,
                        title, title_color, subtitle) {
  
  tbl <- as_tbl_graph(g)
  lay <- create_layout(tbl, layout = "manual",
                       x = layout_coords$x, y = layout_coords$y)
  
  lay$node_size <- 6 + (lay$degree_full - size_range[1]) /
    (size_range[2] - size_range[1]) * 16
  
  lay$label_size <- ifelse(lay$category %in% c("TE of Interest", "TE Hub"), 5,
                           ifelse(lay$category == "SCZ Gene", 4.5, 4))
  
  lay$font_face <- ifelse(lay$color_group != "Gene", "bold", "plain")
  
  node_fill_colors <- c("Gene" = "#63B8FF", "LINE" = "#FF8C00", "LTR" = "#CD950C")
  node_shape_values <- c("Gene" = 21, "LINE" = 23, "LTR" = 23)
  
  p <- ggraph(lay) +
    geom_edge_link(alpha = 0.15, color = "gray35", width = 0.5) +
    geom_node_point(aes(fill = color_group, shape = color_group, size = node_size),
                    color = "gray30", stroke = 0.4) +
    geom_node_text(aes(label = label),
                   size = lay$label_size,
                   repel = TRUE,
                   max.overlaps = 100,
                   box.padding = 0.8,
                   point.padding = 0.5,
                   min.segment.length = 0,
                   force = 2,
                   force_pull = 0.5,
                   segment.color = "gray60",
                   segment.size = 0.3,
                   fontface = lay$font_face,
                   seed = 42) +
    scale_fill_manual(values = node_fill_colors, name = "Node Type") +
    scale_shape_manual(values = node_shape_values, name = "Node Type") +
    scale_size_identity() +
    labs(title = title, subtitle = subtitle) +
    theme_void() +
    theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5, color = title_color),
      plot.subtitle = element_text(size = 18, hjust = 0.5, color = "gray40"),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  return(p)
}

# ---- Build panels ----
p_hc <- build_panel(hc_sub, layout_coords, size_range,
                    "Healthy Control Network", "#4DAF4A",
                    paste0(vcount(hc_sub), " nodes | ", format(ecount(hc_sub), big.mark = ","), " edges"))

p_scz <- build_panel(scz_sub, layout_coords, size_range,
                     "Schizophrenia Network", "#E41A1C",
                     paste0(vcount(scz_sub), " nodes | ", format(ecount(scz_sub), big.mark = ","), " edges (\u221275.7%)"))

p_combined <- p_hc + p_scz +
  plot_layout(ncol = 2) +
  plot_annotation(
    caption = "Node color: Gene (blue) | LINE TE (orange) | LTR TE (gold) | Node shape: Circle (gene) | Diamond (TE) | Size \u221d full network degree",
    theme = theme(
      plot.caption = element_text(size = 14, hjust = 0.5, color = "gray40")
    )
  )

# ---- Save BIG ----
output_dir <- "results/network_analysis/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(paste0(output_dir, "panels_AB_poster.pdf"),
       p_combined, width = 30, height = 16, dpi = 300)
ggsave(paste0(output_dir, "panels_AB_poster.png"),
       p_combined, width = 30, height = 16, dpi = 300)

# Individual panels — same settings, 15x16 each
ggsave(paste0(output_dir, "panel_A_hc_poster.pdf"),
       p_hc, width = 15, height = 16, dpi = 300)
ggsave(paste0(output_dir, "panel_A_hc_poster.png"),
       p_hc, width = 15, height = 16, dpi = 300)

ggsave(paste0(output_dir, "panel_B_scz_poster.pdf"),
       p_scz, width = 15, height = 16, dpi = 300)
ggsave(paste0(output_dir, "panel_B_scz_poster.png"),
       p_scz, width = 15, height = 16, dpi = 300)


