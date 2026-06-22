# =============================================================================        
# SCRIPT 1: GENERATE PANELS for JCI LETTER FIGURE 1        
# JOURNAL SPECS:    
#   - 17.88 cm × 23.88 cm final dimensions    
#   - All text 6–8 pt, Arial/Helvetica, NO bold anywhere    
#   - Export as TIFF at 600 ppi (primary), plus SVG for Inkscape    
#   - Panels: A and B (networks), B (barplots), D (boxplots),   
#             E (volcano), F (heatmap)  
# =============================================================================

message("====== SCRIPT 1: GENERATE PANELS ======")        
suppressPackageStartupMessages({        
  library(ggplot2)        
  library(patchwork)        
  library(igraph)        
  library(ggraph)        
  library(tidygraph)        
  library(pheatmap)        
  library(RColorBrewer)        
  library(ggrepel)        
  library(dplyr)        
  library(tidyr)        
  library(DESeq2)        
  library(scales)        
  library(grid)    
  library(svglite)    
})

# ---- 1.1 Paths ----        
output_dir <- "results/final_figure"        
panel_dir  <- file.path(output_dir, "panels")        
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1.2 Font ----        
FONT_FAMILY <- "Arial"

# ---- 1.3 FIGURE DIMENSIONS (cm → inches) ----        
FIG_W_CM <- 17.88    
FIG_H_CM <- 23.88    
FIG_W <- FIG_W_CM / 2.54   # 7.039 in    
FIG_H <- FIG_H_CM / 2.54   # 9.402 in

# ---- 1.4 TEXT SIZE CONSTANTS (pts) ----    
TEXT_AXIS_TITLE <- 7    
TEXT_AXIS_LABEL <- 6    
TEXT_LEGEND_TITLE <- 6.5    
TEXT_LEGEND_TEXT <- 6    
TEXT_STRIP <- 7    
TEXT_PANEL_TAG <- 8    
TEXT_ANNOTATION <- 5.5    
TEXT_REPEL_LABEL <- 5    
TEXT_BAR_LABEL <- 4.5    
TEXT_NETWORK_LABEL_TE <- 5    
TEXT_NETWORK_LABEL_GENE <- 4.5    
TEXT_SUBTITLE <- 6.5

# ---- 1.5 MASTER PALETTES ----        
color_phenotype <- c("HC" = "#4DAF4A", "SCZ" = "#984EA3")        
color_biotype   <- c("Gene" = "steelblue2", "L1" = "darkorange", "HERV" = "gold3")        
shape_biotype   <- c("Gene" = 21, "L1" = 23, "HERV" = 23)

# ---- 1.6 JCI THEME (no bold anywhere) ----        
theme_jci <- function(base_size = 7) {        
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +        
    theme(        
      plot.title       = element_text(face = "plain", size = TEXT_AXIS_TITLE),      
      plot.subtitle    = element_blank(),      
      axis.title       = element_text(face = "plain", size = TEXT_AXIS_TITLE),        
      axis.text        = element_text(color = "black", size = TEXT_AXIS_LABEL),        
      axis.line        = element_line(linewidth = 0.3),    
      axis.ticks       = element_line(linewidth = 0.2),    
      axis.ticks.length = unit(0.03, "cm"),    
      legend.title     = element_text(face = "plain", size = TEXT_LEGEND_TITLE),        
      legend.text      = element_text(size = TEXT_LEGEND_TEXT),        
      legend.key.size  = unit(0.2, "cm"),        
      legend.margin    = margin(0, 0, 0, 0),        
      legend.spacing   = unit(0.05, "cm"),        
      legend.box.spacing = unit(0.05, "cm"),    
      strip.text       = element_text(face = "plain", size = TEXT_STRIP),        
      strip.background = element_rect(fill = "gray92", linetype = "blank"),        
      plot.tag          = element_text(face = "plain", size = TEXT_PANEL_TAG,     
                                       family = FONT_FAMILY),        
      plot.margin      = margin(2, 2, 2, 2)        
    )        
}

# ---- 1.7 Utilities ----        
classify_te <- function(te_names) {        
  ifelse(grepl("^L1", te_names), "L1", "HERV")        
}

remap_biotype <- function(raw_type) {        
  dplyr::case_when(        
    raw_type == "CG"   ~ "Gene",        
    raw_type == "LINE" ~ "L1",        
    raw_type == "LTR"  ~ "HERV",        
    TRUE               ~ "Gene"        
  )        
}

# ---- 2. LOAD DATA ----        
message("Loading data...")        
de_results    <- read.csv("data/DE_results_all_genes.csv")        
vsd           <- readRDS("data/vsd_object.rds")        
vst_matrix    <- assay(vsd)        
hc_graph      <- readRDS("data/hc_graph_intergenic_only.rds")        
scz_graph     <- readRDS("data/scz_graph_intergenic_only.rds")        
degree_info   <- readRDS("data/degree_info_intergenic_only.rds")        
feature_lists <- readRDS("data/feature_lists_intergenic_only.rds")        
meta          <- read.csv("data/scz_samples_metadata.csv")        
meta$phenotype <- factor(meta$phenotype,        
                         levels = c("control", "schizophrenia"),        
                         labels = c("HC", "SCZ"))        
message("Data loaded.\n")

# ============================================================================        
# ---- 3. PANEL A: NETWORKS ----        
# ============================================================================        
message("--- Panel A: Networks ---")

scz_risk_genes <- c(        
  "MTHFR","RTN4R","COMT","HTR2A","DRD3","SYN2","CHI3L1","NRXN1",        
  "SETD1A","MBD5","PHIP","DAOA","DISC1","DTNBP1","DAO","ABCA13",        
  "SHANK3","AKT1","NPAS3","RELN","ZNF804A","DRD2","ASTN2","RGS4",        
  "PRODH","GABRB2","CSMD1","DLG2","DLGAP2","RBFOX1","PRKN","CTNND2",        
  "DPP6","DPP10","MSRA","PCDH15","DMD","NRG1","DRD4","SLC1A1",        
  "GRIN2B","ZDHHC8","CHRNA7","GRIN2A","PDE4B","GRID2","NOTCH4",        
  "CNTN6","EHMT1","MACROD2","KCNN3","KATNAL2","TPH1","GABRB1",        
  "GRIK3","LPP","PTPRM","NRG3","TRAPPC9","DBH","FHIT","PARD3B",        
  "PDE11A","TMLHE","VPS13B","FZD3","PTPRT","CNP"        
)        
tes_of_interest <- c("L1FLnI_2p13.1f", "L1FLnI_8q13.3e", "HERVL_4q32.1a")

scz_in_network <- intersect(scz_risk_genes, V(hc_graph)$name)

te_node_names <- V(hc_graph)$name[V(hc_graph)$type == "TE"]        
te_deg_hc     <- sort(degree(hc_graph, v = te_node_names), decreasing = TRUE)        
top_te_hubs   <- setdiff(names(te_deg_hc), tes_of_interest)[1:15]

gene_node_names <- V(hc_graph)$name[V(hc_graph)$type == "gene"]        
gene_deg_hc     <- sort(degree(hc_graph, v = gene_node_names), decreasing = TRUE)        
top_gene_hubs   <- setdiff(names(gene_deg_hc), scz_in_network)[1:15]

selected_nodes <- unique(c(tes_of_interest, scz_in_network, top_te_hubs, top_gene_hubs))

hc_sub  <- induced_subgraph(hc_graph,  vids = selected_nodes)        
scz_sub <- induced_subgraph(scz_graph, vids = selected_nodes)

assign_display_attrs <- function(g, full_graph) {        
  nms <- V(g)$name        
  tp  <- V(g)$type        
  V(g)$biotype <- case_when(        
    tp == "gene"       ~ "Gene",        
    grepl("^L1", nms)  ~ "L1",        
    TRUE               ~ "HERV"        
  )        
  V(g)$degree_full <- degree(full_graph, v = nms)        
  V(g)$label <- nms        
  return(g)        
}

hc_sub  <- assign_display_attrs(hc_sub,  hc_graph)        
scz_sub <- assign_display_attrs(scz_sub, scz_graph)

set.seed(42)        
shared_layout <- create_layout(as_tbl_graph(hc_sub), layout = "fr", niter = 5000,        
                               area = vcount(hc_sub)^2 * 8)        
layout_coords <- data.frame(        
  name = shared_layout$name,        
  x    = shared_layout$x,        
  y    = shared_layout$y        
)

build_network_panel <- function(g, layout_coords, full_degree_range,      
                                subtitle_text, title_color) {        
  nms <- V(g)$name        
  idx <- match(nms, layout_coords$name)
  
  tg  <- as_tbl_graph(g)        
  lay <- create_layout(tg, layout = "manual",        
                       x = layout_coords$x[idx],        
                       y = layout_coords$y[idx])
  
  dr <- full_degree_range        
  lay$node_size <- 1.0 + (lay$degree_full - dr[1]) / max(1, dr[2] - dr[1]) * 2.5
  
  label_df <- data.frame(        
    x     = lay$x,        
    y     = lay$y,        
    label = lay$label,        
    is_te = lay$biotype != "Gene",        
    stringsAsFactors = FALSE        
  )        
  label_df$lbl_size <- ifelse(label_df$is_te,     
                              TEXT_NETWORK_LABEL_TE / ggplot2::.pt,     
                              TEXT_NETWORK_LABEL_GENE / ggplot2::.pt)
  
  p <- ggraph(lay) +        
    geom_edge_link(alpha = 0.12, color = "gray50", width = 0.15) +        
    geom_node_point(aes(fill = biotype, shape = biotype, size = node_size),        
                    color = "gray30", stroke = 0.2) +        
    geom_text_repel(        
      data            = label_df,        
      aes(x = x, y = y, label = label),        
      size            = label_df$lbl_size,        
      fontface        = "plain",        
      family          = FONT_FAMILY,        
      box.padding     = 0.15,        
      point.padding   = 0.1,        
      max.overlaps    = 40,        
      force           = 1.0,        
      force_pull      = 0.5,        
      min.segment.length = 0.15,        
      segment.size    = 0.1,        
      segment.color   = "gray60",        
      seed            = 42        
    ) +        
    scale_fill_manual(values = color_biotype, name = "Feature type") +        
    scale_shape_manual(values = shape_biotype, name = "Feature type") +        
    scale_size_identity() +        
    labs(subtitle = subtitle_text) +      
    theme_jci() +        
    theme(        
      axis.text       = element_blank(),        
      axis.title      = element_blank(),        
      axis.ticks      = element_blank(),        
      axis.line       = element_blank(),        
      legend.position = "none",        
      plot.subtitle   = element_text(size = TEXT_SUBTITLE, hjust = 0.5,      
                                     color = title_color, face = "plain"),      
      plot.margin     = margin(1, 1, 1, 1)        
    )        
  return(p)        
}

full_degree_range <- range(c(V(hc_sub)$degree_full, V(scz_sub)$degree_full))

pA_hc <- build_network_panel(        
  hc_sub, layout_coords, full_degree_range,        
  paste0("HC: ", vcount(hc_sub), " nodes, ", format(ecount(hc_sub), big.mark = ","), " edges"),      
  color_phenotype["HC"]        
)        
pA_scz <- build_network_panel(        
  scz_sub, layout_coords, full_degree_range,        
  paste0("SCZ: ", vcount(scz_sub), " nodes, ", format(ecount(scz_sub), big.mark = ","), " edges"),      
  color_phenotype["SCZ"]        
)

saveRDS(pA_hc,  file.path(panel_dir, "panel_A_hc.rds"))        
saveRDS(pA_scz, file.path(panel_dir, "panel_A_scz.rds"))

ggsave(file.path(panel_dir, "panel_A_hc.svg"), pA_hc,    
       device = svglite, width = FIG_W * 0.5, height = FIG_H * 0.35)    
ggsave(file.path(panel_dir, "panel_A_scz.svg"), pA_scz,    
       device = svglite, width = FIG_W * 0.5, height = FIG_H * 0.35)

message("Panel A saved.\n")


# ============================================================================        
# ---- 4. PANEL C: EDGE-TYPE BARPLOTS ----        
# ============================================================================        
message("--- Panel C: Edge-type barplots ---")

count_edges_from_graph <- function(g, condition_label) {        
  el    <- ends(g, E(g), names = TRUE)        
  ntype <- setNames(V(g)$type, V(g)$name)        
  t1 <- ntype[el[, 1]]        
  t2 <- ntype[el[, 2]]
  
  broad_cat <- case_when(        
    t1 == "gene" & t2 == "gene" ~ "Gene-Gene",        
    t1 == "TE"   & t2 == "TE"   ~ "TE-TE",        
    TRUE                         ~ "TE-Gene"        
  )
  
  te_class1 <- ifelse(t1 == "TE", classify_te(el[, 1]), "Gene")        
  te_class2 <- ifelse(t2 == "TE", classify_te(el[, 2]), "Gene")        
  pair_a <- pmin(te_class1, te_class2)        
  pair_b <- pmax(te_class1, te_class2)        
  fine_cat <- paste0(pair_a, "-", pair_b)
  
  gg_count <- sum(broad_cat == "Gene-Gene")        
  gg_df <- data.frame(EdgeType = "Gene-Gene", Count = gg_count,        
                      Condition = condition_label, stringsAsFactors = FALSE)
  
  te_mask <- broad_cat != "Gene-Gene"        
  fine_tbl <- as.data.frame(table(fine_cat[te_mask]), stringsAsFactors = FALSE)        
  names(fine_tbl) <- c("EdgeType", "Count")        
  fine_tbl$Condition <- condition_label
  
  return(list(gg = gg_df, fine = fine_tbl))        
}

message("  Counting HC edges...")        
hc_ec  <- count_edges_from_graph(hc_graph, "HC")        
message("  Counting SCZ edges...")        
scz_ec <- count_edges_from_graph(scz_graph, "SCZ")

gg_df   <- rbind(hc_ec$gg, scz_ec$gg)        
fine_df <- rbind(hc_ec$fine, scz_ec$fine)

gg_df$Condition   <- factor(gg_df$Condition, levels = c("HC", "SCZ"))        
fine_df$Condition <- factor(fine_df$Condition, levels = c("HC", "SCZ"))

te_edge_order <- c("Gene-L1", "Gene-HERV", "HERV-L1", "L1-L1", "HERV-HERV")        
fine_df$EdgeType <- factor(fine_df$EdgeType, levels = te_edge_order)        
fine_df <- fine_df %>% tidyr::complete(EdgeType, Condition, fill = list(Count = 0))

bar_label_size <- TEXT_BAR_LABEL / ggplot2::.pt

pC_left <- ggplot(gg_df, aes(x = Condition, y = Count, fill = Condition)) +        
  geom_col(width = 0.55, color = "gray30", linewidth = 0.2) +        
  geom_text(aes(label = paste0(round(Count / 1e6, 2), "M")),        
            vjust = -0.4, size = bar_label_size, fontface = "plain",     
            family = FONT_FAMILY) +      
  scale_fill_manual(values = color_phenotype) +        
  scale_y_continuous(labels = function(x) paste0(x / 1e6, "M"),        
                     expand = expansion(mult = c(0, 0.18))) +        
  labs(y = "Gene-Gene edges", x = NULL) +      
  theme_jci() +        
  theme(legend.position = "none")

pC_right <- ggplot(fine_df, aes(x = EdgeType, y = Count, fill = Condition)) +        
  geom_col(position = position_dodge(width = 0.7), width = 0.6,        
           color = "gray30", linewidth = 0.15) +        
  geom_text(        
    aes(label = format(Count, big.mark = ",")),        
    position = position_dodge(width = 0.7),        
    vjust = -0.3, size = bar_label_size, family = FONT_FAMILY, fontface = "plain"    
  ) +        
  scale_fill_manual(values = color_phenotype, name = NULL) +        
  scale_y_continuous(labels = label_comma(),        
                     expand = expansion(mult = c(0, 0.20))) +        
  labs(y = "TE-involving edges", x = NULL) +      
  theme_jci() +        
  theme(  
    legend.position = "top",        
    legend.key.size = unit(0.15, "cm"),        
    legend.margin = margin(0, 0, 0, 0),  
    axis.text.x = element_text(  
      size = TEXT_AXIS_LABEL,   
      angle = 45,   
      hjust = 1,   
      vjust = 1  
    )  
  )

saveRDS(pC_left,  file.path(panel_dir, "panel_C_left.rds"))        
saveRDS(pC_right, file.path(panel_dir, "panel_C_right.rds"))  
ggsave(file.path(panel_dir, "panel_C_left.svg"), pC_left,    
       device = svglite, width = 1.8, height = 2.2)    
ggsave(file.path(panel_dir, "panel_C_right.svg"), pC_right,    
       device = svglite, width = 3.0, height = 2.2)    
message("Panel C saved.\n")


# ============================================================================        
# ---- 5. PANEL D: NODE DEGREE BOXPLOT WITH BOOTSTRAP CI ----        
# ============================================================================        
message("--- Panel D: Node degree boxplot with bootstrap test ---")

build_degree_df <- function(g, condition_label) {        
  nms  <- V(g)$name        
  tp   <- V(g)$type        
  degs <- degree(g)        
  biotype <- case_when(        
    tp == "gene"       ~ "Gene",        
    grepl("^L1", nms)  ~ "L1",        
    TRUE               ~ "HERV"        
  )        
  data.frame(Degree = degs, Category = biotype, Phenotype = condition_label,        
             stringsAsFactors = FALSE)        
}

degree_df <- rbind(build_degree_df(hc_graph, "HC"),        
                   build_degree_df(scz_graph, "SCZ"))        
degree_df$Phenotype <- factor(degree_df$Phenotype, levels = c("HC", "SCZ"))        
degree_df$Category  <- factor(degree_df$Category,  levels = c("Gene", "L1", "HERV"))

# ---- Bootstrap test ----      
set.seed(42)      
n_boot <- 100000

boot_results <- lapply(c("HC", "SCZ"), function(pheno) {      
  gene_deg <- degree_df$Degree[degree_df$Category == "Gene" & degree_df$Phenotype == pheno]      
  l1_deg   <- degree_df$Degree[degree_df$Category == "L1"   & degree_df$Phenotype == pheno]      
  herv_deg <- degree_df$Degree[degree_df$Category == "HERV" & degree_df$Phenotype == pheno]  
  
  n_gene <- length(gene_deg)      
  n_l1   <- length(l1_deg)      
  n_herv <- length(herv_deg)  
  
  obs_l1   <- median(l1_deg)   - median(gene_deg)      
  obs_herv <- median(herv_deg) - median(gene_deg)  
  
  boot_l1   <- replicate(n_boot, {      
    bg <- sample(gene_deg, n_gene, replace = TRUE)      
    bl <- sample(l1_deg,   n_l1,   replace = TRUE)      
    median(bl) - median(bg)      
  })  
  
  boot_herv <- replicate(n_boot, {      
    bg <- sample(gene_deg, n_gene, replace = TRUE)      
    bh <- sample(herv_deg, n_herv, replace = TRUE)      
    median(bh) - median(bg)      
  })  
  
  p_l1   <- mean(boot_l1   <= 0)      
  p_herv <- mean(boot_herv <= 0)      
  if (p_l1   == 0) p_l1   <- 1 / n_boot      
  if (p_herv == 0) p_herv <- 1 / n_boot  
  
  rbind(      
    data.frame(      
      Phenotype = pheno, Comparison = "L1 vs Gene",      
      median_Gene = median(gene_deg), median_TE = median(l1_deg),      
      observed_diff = obs_l1,      
      ci_lower = quantile(boot_l1, 0.025), ci_upper = quantile(boot_l1, 0.975),      
      p_value = p_l1, n_Gene = n_gene, n_TE = n_l1,      
      stringsAsFactors = FALSE      
    ),      
    data.frame(      
      Phenotype = pheno, Comparison = "HERV vs Gene",      
      median_Gene = median(gene_deg), median_TE = median(herv_deg),      
      observed_diff = obs_herv,      
      ci_lower = quantile(boot_herv, 0.025), ci_upper = quantile(boot_herv, 0.975),      
      p_value = p_herv, n_Gene = n_gene, n_TE = n_herv,      
      stringsAsFactors = FALSE      
    )      
  )      
})

boot_df <- do.call(rbind, boot_results)    
boot_df$p_label <- ifelse(      
  boot_df$p_value <= 1/n_boot,      
  paste0("p<", format(1/n_boot, scientific = TRUE)),      
  ifelse(boot_df$p_value < 0.001,      
         formatC(boot_df$p_value, format = "e", digits = 1),      
         formatC(boot_df$p_value, format = "f", digits = 3))      
)

message("  Bootstrap results:")      
print(boot_df)    
write.csv(boot_df, file.path(panel_dir, "panel_D_bootstrap_stats.csv"), row.names = FALSE)

# ---- Build bracket annotations with annotate() ----  
# position_dodge(width=0.75) with 3 groups (Gene, L1, HERV):  
#   offset between adjacent groups = 0.75 / 3 = 0.25  
#   Gene = x_base - 0.25,  L1 = x_base,  HERV = x_base + 0.25  
# x_base: HC = 1, SCZ = 2  
dodge_width  <- 0.75  
group_offset <- dodge_width / 3  # 0.25

# Y positions in DATA space (before log10 transform)  
# We need actual degree+1 values, then they get log-transformed by scale_y_log10  
max_val <- max(degree_df$Degree + 1)  
# Use log-space to position brackets evenly in visual space  
y_l1_pos    <- 10^(log10(max_val) * 0.85)  
y_herv_pos  <- 10^(log10(max_val) * 0.95)  
# Tick length in log space  
tick_factor <- 10^(log10(max_val) * 0.02)

anno_text_size <- TEXT_ANNOTATION / ggplot2::.pt

# Start building the plot  
pD <- ggplot(degree_df, aes(x = Phenotype, y = Degree + 1, fill = Category)) +        
  geom_boxplot(outlier.shape = NA, linewidth = 0.2,        
               position = position_dodge(width = dodge_width), width = 0.65)

# Add brackets for each phenotype  
for (pheno in c("HC", "SCZ")) {  
  x_base <- ifelse(pheno == "HC", 1, 2)  
  x_gene <- x_base - group_offset  
  x_l1   <- x_base  
  x_herv <- x_base + group_offset  
  
  p_lab_l1   <- boot_df$p_label[boot_df$Phenotype == pheno &   
                                  boot_df$Comparison == "L1 vs Gene"]  
  p_lab_herv <- boot_df$p_label[boot_df$Phenotype == pheno &   
                                  boot_df$Comparison == "HERV vs Gene"]  
  
  # L1 vs Gene bracket  
  pD <- pD +  
    annotate("segment", x = x_gene, xend = x_l1,   
             y = y_l1_pos, yend = y_l1_pos, linewidth = 0.2) +  
    annotate("segment", x = x_gene, xend = x_gene,  
             y = y_l1_pos / tick_factor, yend = y_l1_pos, linewidth = 0.2) +  
    annotate("segment", x = x_l1, xend = x_l1,  
             y = y_l1_pos / tick_factor, yend = y_l1_pos, linewidth = 0.2) +  
    annotate("text", x = (x_gene + x_l1) / 2, y = y_l1_pos * 1.2,  
             label = p_lab_l1, size = anno_text_size,   
             family = FONT_FAMILY, fontface = "plain")  
  
  # HERV vs Gene bracket  
  pD <- pD +  
    annotate("segment", x = x_gene, xend = x_herv,  
             y = y_herv_pos, yend = y_herv_pos, linewidth = 0.2) +  
    annotate("segment", x = x_gene, xend = x_gene,  
             y = y_herv_pos / tick_factor, yend = y_herv_pos, linewidth = 0.2) +  
    annotate("segment", x = x_herv, xend = x_herv,  
             y = y_herv_pos / tick_factor, yend = y_herv_pos, linewidth = 0.2) +  
    annotate("text", x = (x_gene + x_herv) / 2, y = y_herv_pos * 1.2,  
             label = p_lab_herv, size = anno_text_size,  
             family = FONT_FAMILY, fontface = "plain")  
}

pD <- pD +  
  scale_fill_manual(values = color_biotype, name = NULL,        
                    guide = guide_legend(override.aes = list(shape = NA))) +        
  scale_y_log10(labels = label_comma(),  
                expand = expansion(mult = c(0.05, 0.15))) +        
  labs(y = "Degree (log)", x = NULL) +      
  theme_jci() +        
  theme(        
    legend.position  = "top",        
    legend.direction = "horizontal",        
    legend.key.size  = unit(0.15, "cm"),        
    legend.margin    = margin(0, 0, 0, 0)        
  )

saveRDS(pD, file.path(panel_dir, "panel_D.rds"))    
ggsave(file.path(panel_dir, "panel_D.svg"), pD,    
       device = svglite, width = 2.5, height = 2.8)    
message("panel D saved.\n")  


# ============================================================================        
# ---- 6. PANEL E: VOLCANO PLOT ----        
# ============================================================================        
message("--- Panel E: Volcano plot ---")

de <- de_results %>%        
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%        
  mutate(        
    Biotype = factor(remap_biotype(gene_type), levels = c("Gene", "L1", "HERV")),        
    neg_log10_padj = -log10(pmax(padj, 1e-300)),        
    Significance = case_when(        
      padj < 0.01 & abs(log2FoldChange) >= 0.5 ~ "Significant",        
      TRUE ~ "NS"        
    ),        
    combined_score = neg_log10_padj * abs(log2FoldChange)        
  )

get_top <- function(data, bio, dir, n) {        
  sub <- data %>% filter(Biotype == bio, Significance == "Significant")        
  if (dir == "up")   sub <- sub %>% filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange))        
  if (dir == "down") sub <- sub %>% filter(log2FoldChange < 0) %>% arrange(log2FoldChange)        
  head(sub$gene_name, n)        
}

label_names <- unique(c(        
  get_top(de, "Gene", "up", 3), get_top(de, "Gene", "down", 2),        
  get_top(de, "L1",   "up", 1), get_top(de, "L1",   "down", 3),        
  get_top(de, "HERV", "up", 1), get_top(de, "HERV", "down", 2),        
  head(de %>% filter(Significance == "Significant") %>%        
         arrange(padj) %>% pull(gene_name), 3)        
))        
de$show_label <- de$gene_name %in% label_names

repel_label_size <- TEXT_REPEL_LABEL / ggplot2::.pt

pE <- ggplot(de, aes(x = log2FoldChange, y = neg_log10_padj)) +        
  geom_point(data = de %>% filter(Significance == "NS"),        
             color = "gray82", size = 0.15, alpha = 0.3) +        
  geom_point(data = de %>% filter(Significance == "Significant"),        
             aes(color = Biotype), size = 0.35, alpha = 0.65) +        
  geom_hline(yintercept = -log10(0.01), linetype = "dashed",        
             color = "gray40", linewidth = 0.2) +        
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",        
             color = "gray40", linewidth = 0.2) +        
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dotted",        
             color = "gray55", linewidth = 0.15) +        
  geom_text_repel(        
    data          = de %>% filter(show_label),        
    aes(label = gene_name, color = Biotype),        
    size          = repel_label_size,      
    fontface      = "italic",        
    family        = FONT_FAMILY,        
    box.padding   = 0.15,        
    point.padding = 0.1,        
    max.overlaps  = 25,        
    force         = 1.5,        
    force_pull    = 0.5,        
    min.segment.length = 0,        
    segment.size  = 0.1,        
    segment.color = "gray50",        
    seed          = 42,        
    show.legend   = FALSE        
  ) +        
  scale_color_manual(values = color_biotype, name = NULL) +        
  scale_x_continuous(breaks = seq(-4, 4, by = 1)) +        
  labs(        
    x = expression(log[2] ~ "fold change"),        
    y = expression(-log[10] ~ "(adj. p-value)")        
  ) +        
  coord_cartesian(xlim = c(-4.5, 4.5)) +        
  theme_jci() +        
  theme(        
    legend.position  = "top",        
    legend.direction = "horizontal",        
    legend.key.size  = unit(0.15, "cm"),        
    legend.margin    = margin(0, 0, 0, 0)        
  )

saveRDS(pE, file.path(panel_dir, "panel_E.rds"))    
ggsave(file.path(panel_dir, "panel_E.svg"), pE,    
       device = svglite, width = 3.5, height = 2.8)    
message("Panel E saved.\n")


# ============================================================================        
# ---- 7. PANEL F: HEATMAP ----        
# ============================================================================        
message("--- Panel F: Heatmap ---")

de_for_heatmap <- de_results %>%        
  filter(!is.na(padj), padj < 0.01) %>%        
  mutate(combined_score = -log10(pmax(padj, 1e-300)) * abs(log2FoldChange)) %>%        
  arrange(desc(combined_score))

top50_names <- head(de_for_heatmap$gene_name, 50)        
available   <- top50_names[top50_names %in% rownames(vst_matrix)]        
message(sprintf("  %d of 50 features in VST matrix", length(available)))

heatmap_mat   <- vst_matrix[available, ]        
heatmap_mat_z <- t(scale(t(heatmap_mat)))        
heatmap_mat_z[heatmap_mat_z >  3] <-  3        
heatmap_mat_z[heatmap_mat_z < -3] <- -3

sample_pheno <- meta$phenotype[match(colnames(heatmap_mat_z), meta$sample)]        
names(sample_pheno) <- colnames(heatmap_mat_z)

# Defensive check  
stopifnot(  
  "Some samples in heatmap matrix not found in metadata" =   
    !any(is.na(sample_pheno))  
)

col_anno <- data.frame(        
  Phenotype = factor(sample_pheno, levels = c("SCZ", "HC")),        
  row.names = colnames(heatmap_mat_z)        
)

row_types_raw     <- de_results$gene_type[match(available, de_results$gene_name)]        
row_types_display <- factor(remap_biotype(row_types_raw), levels = c("Gene", "L1", "HERV"))        
row_anno <- data.frame(Type = row_types_display, row.names = available)

anno_colors <- list(        
  Phenotype = color_phenotype,        
  Type      = color_biotype        
)

heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# ---- Cell dimensions: fit within figure width ----  
n_rows <- length(available)  
n_cols <- ncol(heatmap_mat_z)

overhead_width_pt  <- 15 + 15 + 80  # tree + row_anno + legend + row labels  
available_width_pt <- (FIG_W * 72) - overhead_width_pt  
cellwidth_val  <- floor(available_width_pt / n_cols)  
cellwidth_val  <- max(cellwidth_val, 2)  # minimum 2pt  
cellwidth_val  <- min(cellwidth_val, 5)  # maximum 5pt

overhead_height_pt  <- 15 + 15 + 25  # col tree + col anno + title  
available_height_pt <- (FIG_H * 0.38 * 72) - overhead_height_pt  
cellheight_val <- floor(available_height_pt / n_rows)  
cellheight_val <- max(cellheight_val, 4)  # minimum 4pt so row labels fit  
cellheight_val <- min(cellheight_val, 7)  # maximum 7pt

message(sprintf("  Heatmap cell dimensions: h=%.0f x w=%.0f pt (%d rows x %d cols)",  
                cellheight_val, cellwidth_val, n_rows, n_cols))

heatmap_obj <- pheatmap(        
  heatmap_mat_z,        
  color                    = heatmap_colors,        
  clustering_distance_rows = "euclidean",        
  clustering_distance_cols = "euclidean",        
  clustering_method        = "complete",        
  cluster_rows             = TRUE,        
  cluster_cols             = TRUE,        
  annotation_col           = col_anno,        
  annotation_row           = row_anno,        
  annotation_colors        = anno_colors,        
  border_color             = NA,        
  show_colnames            = FALSE,        
  show_rownames            = TRUE,        
  fontsize_row             = 5,          
  fontsize                 = 6,          
  fontsize_col             = 6,        
  fontfamily               = FONT_FAMILY,    
  annotation_names_row     = TRUE,        
  annotation_names_col     = TRUE,        
  annotation_legend        = TRUE,        
  treeheight_row           = 12,         
  treeheight_col           = 12,        
  legend                   = TRUE,        
  legend_breaks            = c(-2, -1, 0, 1, 2),        
  legend_labels            = c("-2", "-1", "0", "1", "2"),        
  cellwidth                = cellwidth_val,        
  cellheight               = cellheight_val,        
  main                     = "Top 50 differentially expressed genes, L1s, and HERVs",        
  silent                   = TRUE        
)

gt_heatmap <- heatmap_obj$gtable

# title to 7pt plain Arial  
for (i in seq_along(gt_heatmap$grobs)) {  
  if (inherits(gt_heatmap$grobs[[i]], "text")) {  
    gp_current <- gt_heatmap$grobs[[i]]$gp  
    if (!is.null(gp_current$fontsize) && gp_current$fontsize > 8) {  
      gt_heatmap$grobs[[i]]$gp$fontsize   <- 7  
      gt_heatmap$grobs[[i]]$gp$fontface   <- "plain"  
      gt_heatmap$grobs[[i]]$gp$fontfamily <- FONT_FAMILY  
    }  
  }  
}

saveRDS(gt_heatmap, file.path(panel_dir, "panel_F_grob.rds"))        
saveRDS(heatmap_obj, file.path(panel_dir, "panel_F_pheatmap.rds"))

# Save standalone heatmap SVG  
svglite(file.path(panel_dir, "panel_F.svg"),    
        width = FIG_W, height = FIG_H * 0.38)    
grid.draw(gt_heatmap)    
dev.off()

message("Panel F saved.\n")  

message("========================================")        
message("All panels in: ", panel_dir)        
message("========================================")  