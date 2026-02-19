# ============================================================  
# PANEL C: TE Rewiring - BIGGER TEXT VERSION  
# ============================================================

library(igraph)  
library(ggplot2)  
library(dplyr)  
library(patchwork)

tes_oi_expanded <- c("L1FLnI_2p13.1f", "HERV30_4q13.2c", "HERVL_4q32.1a",  
                     "L1FLnI_8q13.3e", "HERVFH21_16p11.2")  
scz_in_network <- intersect(scz_risk_genes, V(hc_graph)$name)

# ---- Collect all edges ----  
all_edges <- list()

for(te in tes_oi_expanded) {  
  hc_neighbors <- intersect(neighbors(hc_graph, te)$name, scz_in_network)  
  scz_neighbors <- intersect(neighbors(scz_graph, te)$name, scz_in_network)  
  all_neighbors <- union(hc_neighbors, scz_neighbors)  
  
  for(gene in all_neighbors) {  
    in_hc <- gene %in% hc_neighbors  
    in_scz <- gene %in% scz_neighbors  
    
    edge_type <- case_when(  
      in_hc & in_scz ~ "Conserved",  
      in_hc & !in_scz ~ "Lost in SCZ",  
      !in_hc & in_scz ~ "Gained in SCZ"  
    )  
    
    mi_val <- NA_real_  
    if(in_scz) {  
      eid <- get_edge_ids(scz_graph, c(te, gene))  
      if(eid > 0) mi_val <- E(scz_graph)$weight[eid]  
    } else if(in_hc) {  
      eid <- get_edge_ids(hc_graph, c(te, gene))  
      if(eid > 0) mi_val <- E(hc_graph)$weight[eid]  
    }  
    
    all_edges[[length(all_edges) + 1]] <- data.frame(  
      te = te, gene = gene, edge_type = edge_type, mi = mi_val,  
      stringsAsFactors = FALSE  
    )  
  }  
}

edges_df <- do.call(rbind, all_edges)

# ---- Edge colors ----  
edge_colors <- c(  
  "Conserved" = "#777777",  
  "Lost in SCZ" = "#4DAF4A",  
  "Gained in SCZ" = "#E41A1C"  
)

te_fill <- ifelse(grepl("^L1", tes_oi_expanded), "#FF8C00", "#CD950C")  
names(te_fill) <- tes_oi_expanded

# ---- Build one radial plot per TE ----  
build_star_plot <- function(te_name, edges_df, te_fill_color) {  
  
  te_edges <- edges_df %>% filter(te == te_name)  
  n <- nrow(te_edges)  
  if(n == 0) return(NULL)  
  
  te_edges <- te_edges %>%  
    mutate(edge_order = case_when(  
      edge_type == "Gained in SCZ" ~ 1,  
      edge_type == "Conserved" ~ 2,  
      edge_type == "Lost in SCZ" ~ 3  
    )) %>%  
    arrange(edge_order, gene)  
  
  angles <- seq(0, 2 * pi, length.out = n + 1)[1:n]  
  radius <- 1  
  
  gene_x <- cos(angles) * radius  
  gene_y <- sin(angles) * radius  
  
  nodes <- data.frame(  
    x = c(0, gene_x),  
    y = c(0, gene_y),  
    label = c(te_name, te_edges$gene),  
    type = c("TE", rep("Gene", n)),  
    stringsAsFactors = FALSE  
  )  
  
  mi_range <- range(edges_df$mi, na.rm = TRUE)  
  
  segments <- data.frame(  
    x = 0, y = 0,  
    xend = gene_x, yend = gene_y,  
    edge_type = te_edges$edge_type,  
    mi = te_edges$mi,  
    stringsAsFactors = FALSE  
  )  
  segments$width <- 0.8 + (segments$mi - mi_range[1]) / (mi_range[2] - mi_range[1]) * 2.5  
  segments$width[is.na(segments$width)] <- 1  
  
  n_gained <- sum(te_edges$edge_type == "Gained in SCZ")  
  n_lost <- sum(te_edges$edge_type == "Lost in SCZ")  
  n_conserved <- sum(te_edges$edge_type == "Conserved")  
  
  subtitle_parts <- c()  
  if(n_gained > 0) subtitle_parts <- c(subtitle_parts, paste0("+", n_gained, " gained"))  
  if(n_lost > 0) subtitle_parts <- c(subtitle_parts, paste0("-", n_lost, " lost"))  
  if(n_conserved > 0) subtitle_parts <- c(subtitle_parts, paste0(n_conserved, " conserved"))  
  subtitle <- paste(subtitle_parts, collapse = " | ")  
  
  # FIXED: push labels further out  
  label_x <- c(0, cos(angles) * (radius + 0.45))  
  label_y <- c(0, sin(angles) * (radius + 0.45))  
  
  label_h <- ifelse(cos(c(0, angles)) > 0.1, 0,  
                    ifelse(cos(c(0, angles)) < -0.1, 1, 0.5))  
  label_h[1] <- 0.5  
  
  p <- ggplot() +  
    geom_segment(data = segments,  
                 aes(x = x, y = y, xend = xend, yend = yend, color = edge_type,  
                     linewidth = width),  
                 alpha = 0.7, show.legend = FALSE) +  
    scale_color_manual(values = edge_colors) +  
    scale_linewidth_identity() +  
    geom_point(data = nodes %>% filter(type == "Gene"),  
               aes(x = x, y = y),  
               shape = 21, size = 10, fill = "#63B8FF", color = "gray30", stroke = 0.5) +  
    geom_point(data = nodes %>% filter(type == "TE"),  
               aes(x = x, y = y),  
               shape = 23, size = 14, fill = te_fill_color, color = "gray30", stroke = 0.5) +  
    geom_text(aes(x = label_x[-1], y = label_y[-1], label = nodes$label[-1],  
                  hjust = label_h[-1]),  
              size = 5.5, fontface = "plain") +  
    geom_text(aes(x = 0, y = -0.18, label = te_name),  
              size = 5, fontface = "bold", vjust = 1) +  
    # FIXED: wider limits  
    coord_equal(xlim = c(-2.2, 2.2), ylim = c(-2.2, 2.2)) +  
    labs(subtitle = subtitle) +  
    theme_void() +  
    theme(  
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray30",  
                                   margin = margin(b = 5)),  
      plot.margin = margin(5, 5, 5, 5)  
    )  
  
  return(p)  
}  


# ---- Build all 5 panels ----  
star_plots <- list()  
for(i in seq_along(tes_oi_expanded)) {  
  te <- tes_oi_expanded[i]  
  star_plots[[i]] <- build_star_plot(te, edges_df, te_fill[te])  
}

# ---- Combine: 1 row ----  
p_stars <- wrap_plots(star_plots, nrow = 1) +  
  plot_annotation(  
    title = "TE Regulatory Rewiring to Schizophrenia Genes",  
    subtitle = "Each panel shows one TE's connections to SCZ genes across conditions",  
    caption = "Edge color:  
Red = gained in SCZ | Green = lost in SCZ | Gray = conserved  
Edge width proportional to mutual information | Diamond = TE | Circle = SCZ gene",  
    theme = theme(  
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5),  
      plot.subtitle = element_text(size = 18, hjust = 0.5, color = "gray40"),  
      plot.caption = element_text(size = 14, hjust = 0.5, color = "gray40",  
                                  lineheight = 1.3)  
    )  
  )

# ---- Combine: 2 rows ----  
p_stars_2row <- wrap_plots(star_plots, nrow = 2) +  
  plot_annotation(  
    title = "TE Regulatory Rewiring to Schizophrenia Genes",  
    subtitle = "Each panel shows one TE's connections to SCZ genes across conditions",  
    caption = "Edge color:  
Red = gained in SCZ | Green = lost in SCZ | Gray = conserved  
Edge width proportional to mutual information | Diamond = TE | Circle = SCZ gene",  
    theme = theme(  
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5),  
      plot.subtitle = element_text(size = 18, hjust = 0.5, color = "gray40"),  
      plot.caption = element_text(size = 14, hjust = 0.5, color = "gray40",  
                                  lineheight = 1.3)  
    )  
  )

# ---- Save ----  
output_dir <- "results/network_analysis/figures/"

ggsave(paste0(output_dir, "panel_C_star_rewiring.pdf"),  
       p_stars, width = 30, height = 8, dpi = 300)  
ggsave(paste0(output_dir, "panel_C_star_rewiring.png"),  
       p_stars, width = 30, height = 8, dpi = 300)

ggsave(paste0(output_dir, "panel_C_star_rewiring_2row.pdf"),  
       p_stars_2row, width = 18, height = 14, dpi = 300)  
ggsave(paste0(output_dir, "panel_C_star_rewiring_2row.png"),  
       p_stars_2row, width = 18, height = 14, dpi = 300)

cat("Panel C saved!\n")  
