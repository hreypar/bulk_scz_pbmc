# =============================================================================      
# SCRIPT 2: ASSEMBLE FIGURE 1 for JCI LETTER      
#   - 17.88 cm × 23.88 cm      
#   - All text 6–8 pt, Arial/Helvetica, no bold      
#   - Primary export: TIFF at 600 ppi      
#   - Panels: A (networks HC+SCZ), C (barplots),    
#             D (boxplot), E (volcano), F (heatmap)    
# =============================================================================

message("====== SCRIPT 2: ASSEMBLE FIGURE 1 ======")

suppressPackageStartupMessages({      
  library(ggplot2)      
  library(patchwork)      
  library(grid)      
  library(grDevices)      
  library(svglite)      
})

# ---- Paths ------------------------------------------------------------------      
panel_dir  <- "results/final_figure/panels"      
output_dir <- "results/final_figure"      
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Constants --------------------------------------------------------------      
FONT_FAMILY <- "Arial"      
FIG_W_CM <- 17.88      
FIG_H_CM <- 23.88      
FIG_W_IN <- FIG_W_CM / 2.54   # 7.039 in      
FIG_H_IN <- FIG_H_CM / 2.54   # 9.402 in

# ---- Load panels ------------------------------------------------------------      
message("Loading panels...")      
pA_hc    <- readRDS(file.path(panel_dir, "panel_A_hc.rds"))      
pA_scz   <- readRDS(file.path(panel_dir, "panel_A_scz.rds"))      
pC_left  <- readRDS(file.path(panel_dir, "panel_C_left.rds"))      
pC_right <- readRDS(file.path(panel_dir, "panel_C_right.rds"))      
pD       <- readRDS(file.path(panel_dir, "panel_D.rds"))      
pE       <- readRDS(file.path(panel_dir, "panel_E.rds"))      
gt_heatmap <- readRDS(file.path(panel_dir, "panel_F_grob.rds"))      
message("Panels loaded.\n")

# ---- Wrap heatmap grob ------------------------------------------------------      
pF <- wrap_elements(full = gt_heatmap)

# ---- Row heights (3 rows) ---------------------------------------------------      
ROW_HEIGHTS <- c(0.40, 0.20, 0.40)

# ---- Remove duplicate legends from volcano ----------------------------------      
pE_noleg <- pE + theme(legend.position = "none")

# ---- Assemble with explicit tags --------------------------------------------    
row1 <- pA_hc + pA_scz +     
  plot_layout(widths = c(1, 1))

row2 <- pC_left + pC_right + pD + pE_noleg +      
  plot_layout(widths = c(0.5, 1.0, 0.85, 1.6))

figure1 <- (      
  row1 /      
    row2 /      
    pF      
) +      
  plot_layout(heights = ROW_HEIGHTS) +      
  plot_annotation(      
    tag_levels = list(c("A", "", "C", "", "D", "E", "F")),    
    theme = theme(      
      plot.tag = element_text(      
        face   = "plain",      
        size   = 8,      
        family = FONT_FAMILY      
      ),      
      plot.margin = margin(4, 3, 3, 3)      
    )      
  )

# ---- Save TIFF ------------------------------      
message("Saving TIFF at 600 ppi...")      
tiff_path <- file.path(output_dir, "Figure1.tiff")      
tiff(      
  filename    = tiff_path,      
  width       = FIG_W_IN,      
  height      = FIG_H_IN,      
  units       = "in",      
  res         = 600,      
  compression = "lzw",      
  type        = "cairo"      
)      
print(figure1)      
dev.off()      
message(sprintf("  TIFF saved: %s (%.1f MB)",       
                tiff_path, file.size(tiff_path) / 1e6))

# ---- Save SVG ----------------------------------------      
message("Saving SVG...")      
svg_path <- file.path(output_dir, "Figure1.svg")      
ggsave(      
  filename = svg_path,      
  plot     = figure1,      
  device   = svglite,      
  width    = FIG_W_IN,      
  height   = FIG_H_IN,      
  units    = "in"      
)      
message(sprintf("  SVG saved: %s", svg_path))

# ---- save PDF ------------------------------------------      
message("Saving PDF...")      
pdf_path <- file.path(output_dir, "Figure1.pdf")      
ggsave(      
  filename = pdf_path,      
  plot     = figure1,      
  device   = cairo_pdf,      
  width    = FIG_W_IN,      
  height   = FIG_H_IN,      
  units    = "in"      
)      
message(sprintf("  PDF saved: %s", pdf_path))

# ---- Save PNG --------------------------------------      
message("Saving PNG...")      
png_path <- file.path(output_dir, "Figure1.png")      
ggsave(      
  filename = png_path,      
  plot     = figure1,      
  device   = "png",      
  width    = FIG_W_IN,      
  height   = FIG_H_IN,      
  units    = "in",      
  dpi      = 600      
)      
message(sprintf("  PNG saved: %s", png_path))

message("========================================")      
message("Figure 1 in: ", output_dir)      
message("=======================================")  