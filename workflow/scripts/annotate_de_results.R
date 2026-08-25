#!/usr/bin/env Rscript

# Annotate differential expression results with GENCODE and TE annotations

library(optparse)  
library(tidyverse)

# Define command line options  
option_list <- list(  
  make_option(c("-i", "--input"), type = "character", default = NULL,  
              help = "Input DE results CSV file", metavar = "FILE"),  
  make_option(c("-o", "--output"), type = "character", default = NULL,  
              help = "Output annotated CSV file", metavar = "FILE"),  
  make_option(c("-g", "--gencode"), type = "character",   
              default = "resources/gencode.v38",  
              help = "Path to GENCODE annotation directory [default: %default]",   
              metavar = "DIR"),  
  make_option(c("-r", "--retro"), type = "character",  
              default = "resources/retro.hg38.v1.rds",  
              help = "Path to retro RDS file [default: %default]",  
              metavar = "FILE"),  
  make_option(c("-t", "--te_meta"), type = "character",  
              default = "resources/TE_annotation.v2.0.tsv",  
              help = "Path to TE metadata TSV file [default: %default]",  
              metavar = "FILE")  
)

opt_parser <- OptionParser(option_list = option_list)  
opt <- parse_args(opt_parser)

# Check required arguments  
if (is.null(opt$input)) {  
  print_help(opt_parser)  
  stop("Input file must be specified with -i/--input", call. = FALSE)  
}

if (is.null(opt$output)) {  
  print_help(opt_parser)  
  stop("Output file must be specified with -o/--output", call. = FALSE)  
}

# Print parameters  
cat("\n=== Annotating DE Results ===\n")  
cat("Input file:", opt$input, "\n")  
cat("Output file:", opt$output, "\n")  
cat("GENCODE directory:", opt$gencode, "\n")  
cat("Retro RDS:", opt$retro, "\n")  
cat("TE metadata:", opt$te_meta, "\n\n")

# Load DE results  
cat("Loading DE results...\n")  
de_results <- read_csv(opt$input, show_col_types = FALSE)  
cat("  Loaded", nrow(de_results), "features\n")

# Load GENCODE v38 metadata  
cat("Loading GENCODE annotations...\n")  
gene_features <- readRDS(file.path(opt$gencode, "metadata.gene_features.rds"))  
gid_tid <- readRDS(file.path(opt$gencode, "metadata.gid_tid.rds"))

# Prepare GENCODE annotations  
cat("Processing GENCODE annotations...\n")

# Fix gene_features: set proper column names from first row  
colnames(gene_features) <- as.character(gene_features[1, ])  
gene_features <- gene_features %>%  
  rownames_to_column("gene_id") %>%  
  filter(gene_id != "gene_id") %>%  
  mutate(across(c(start, end), as.numeric))

# Fix gid_tid: set column name  
colnames(gid_tid) <- "transcript_ids"  
transcripts_per_gene <- gid_tid %>%  
  rownames_to_column("gene_id") %>%  
  filter(gene_id != "tid.txt.gz") %>%  
  mutate(n_transcripts = str_count(transcript_ids, ",") + 1)

# Define primary chromosomes  
primary_chrs <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

gencode_annot <- gene_features %>%  
  select(gene_id, gene_name, chrom, start, end, strand) %>%  
  rename(chr = chrom) %>%  
  left_join(transcripts_per_gene, by = "gene_id") %>%  
  # Prioritize primary assembly over patches/alternative loci  
  mutate(is_primary = chr %in% primary_chrs) %>%  
  group_by(gene_name) %>%  
  arrange(desc(is_primary), gene_id) %>%  
  slice(1) %>%  # Keep first entry per gene_name (primary if exists, else patch)  
  ungroup() %>%  
  select(-is_primary)

cat("  Kept", nrow(gencode_annot), "unique genes from GENCODE\n")

# Load TE annotations  
cat("Loading TE annotations...\n")  
retro <- readRDS(opt$retro)  
meta_retro <- read.table(opt$te_meta, header = TRUE)

# Prepare TE annotations  
cat("Processing TE annotations...\n")  
te_annot <- retro %>%  
  left_join(meta_retro, by = c("locus" = "Locus")) %>%  
  select(  
    locus,  
    chr = chrom,  
    start,  
    end,  
    strand,  
    family,  
    category,  
    te_class = Class,  
    te_coding = TE_CODING,  
    te_type = TE_type,  
    closest_upstream = ClosestUpstream_gn,  
    closest_downstream = ClosestDownstream_gn  
  ) %>%  
  rename(gene_name = locus)

# Annotate DE results  
cat("Merging annotations with DE results...\n")  
de_annotated <- de_results %>%  
  # First try to match with GENCODE (CG features)  
  left_join(gencode_annot, by = "gene_name") %>%  
  # Then fill in TE annotations for LINE/LTR features  
  left_join(  
    te_annot,  
    by = "gene_name",  
    suffix = c("", "_te")  
  ) %>%  
  # Combine the coordinate columns (use GENCODE if available, otherwise TE)  
  mutate(  
    chr = coalesce(chr, chr_te),  
    start = coalesce(start, start_te),  
    end = coalesce(end, end_te),  
    strand = coalesce(strand, strand_te)  
  ) %>%  
  # Clean up duplicate columns  
  select(-ends_with("_te")) %>%  
  # Reorder columns  
  select(  
    gene_id, gene_name, gene_type,   
    chr, start, end, strand,  
    n_transcripts, transcript_ids,  
    # TE-specific annotations (will be NA for CG features)  
    family, category, te_class, te_coding, te_type,  
    closest_upstream, closest_downstream,  
    # DE stats  
    baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, direction  
  )

# Save annotated results  
cat("Writing output...\n")  
write_csv(de_annotated, opt$output)

# Summary statistics  
cat("\n=== Annotation Summary ===\n")  
cat("Total features:", nrow(de_annotated), "\n")  
cat("  Input features:", nrow(de_results), "\n")  
if (nrow(de_annotated) != nrow(de_results)) {  
  cat("  WARNING: Row count mismatch!\n")  
}  
cat("\nBy annotation source:\n")  
cat("  GENCODE genes (with gene_id):", sum(!is.na(de_annotated$gene_id)), "\n")  
cat("  TE features (with family):", sum(!is.na(de_annotated$family)), "\n")  
cat("  Missing coordinates:", sum(is.na(de_annotated$chr)), "\n\n")

cat("By gene_type:\n")  
print(table(de_annotated$gene_type, useNA = "ifany"))

cat("\nDone! Output written to:", opt$output, "\n")  

