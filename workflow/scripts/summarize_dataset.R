#! /usr/bin/env Rscript --vanilla

if (!requireNamespace("scopetools", quietly = TRUE))
    devtools::install_github("nixonlab/scopetools")

suppressPackageStartupMessages({
    library(tidyverse)
    library(scopetools)
})


importCounts <- function(
    files.tele,
    files.star,
    removePAR = TRUE,
    useGeneNames = TRUE,
    gid2gname.rds = NULL
    ) {
    
    # Reorder files to match
    if(!all(names(files.tele) == names(files.star)))
        files.star <- files.star[names(files.tele)]
    
    # Check file order
    sampids <- names(files.tele)
    stopifnot(all(sampids == names(files.tele)))
    stopifnot(all(sampids == names(files.star)))

    # Check that both files exist
    comp <- file.exists(files.star) & file.exists(files.star)
    if(!all(comp)) {
        cat('WARNING: Missing outputs: ', paste(sampids[!comp], sep=','), '\n')
        sampids <- sampids[comp]
        files.star <- files.star[comp]
        files.tele <- files.tele[comp]
    }
    cat(sprintf('Combining %d samples\n', length(sampids)))

    counts.te <- scopetools::load_telescope_reports(
        files.tele, 
        all_locs = scopetools::retro.hg38.v1$locus
    )
    
    feats.te <- scopetools::retro.hg38.v1 %>%
        dplyr::select(gene.name=locus, gene.type=te_class) %>%
        dplyr::mutate_all(as.character)
    
    stopifnot(all(rownames(counts.te) == feats.te$gene.name))
    
    
    counts.cg <- scopetools::load_star_counts(
        files.star,
        stranded="unstranded"
    )

    if(removePAR) {
        is_par <- grepl('_PAR_Y$', rownames(counts.cg))
        if(sum(counts.cg[is_par,])!=0)
            cat('WARNING: PAR counts is not zero\n')
        counts.cg <- counts.cg[!is_par,]
        rm(is_par)
    }
    
    if(useGeneNames) {
        if(!is.null(gid2gname.rds)) {
            if(file.exists(gid2gname.rds)) {
                gid2gname <- tibble::rownames_to_column(readRDS(gid2gname.rds), 'gid')
                names(gid2gname) <- c('gid', 'gname')

                newlabels <- data.frame(gid = rownames(counts.cg)) %>%
                    dplyr::left_join(gid2gname, by = 'gid') %>%
                    dplyr::mutate(unique.gname = make.unique(gname))

                rownames(counts.cg) <- newlabels$unique.gname
                rm(gid2gname, newlabels)
            }
        }
    }
    
    feats.cg <- data.frame(
        gene.name = rownames(counts.cg),
        gene.type = 'CG'
    )
    
    stopifnot(all(colnames(counts.cg) == colnames(counts.te)))
    counts <- as.matrix(rbind(counts.cg, counts.te))
    
    stopifnot(all(colnames(feats.cg) == colnames(feats.te)))
    feats <- rbind(feats.cg, feats.te)
    
    stopifnot(all(rownames(counts) == feats$gene.name))
    
    return(list(counts=counts, feats=feats))
}


################################################################################
#   Parse input
################################################################################
if(exists("snakemake")) {
    addNames <- function(v) {
        names(v) <- stringr::str_split_i(v, '/', 3)
        v
    }
    
    ret <- importCounts(
        addNames(snakemake@input[['telescope_tsv']]),
        addNames(snakemake@input[['star_tsv']]),
        snakemake@params[["removePAR"]],
        snakemake@params[["useGeneNames"]],
        snakemake@params[["gid2gname_rds"]]
    )
    
    saveRDS(ret$counts, snakemake@output[['counts_rds']])
    saveRDS(ret$feats, snakemake@output[['features_rds']])
} else {
    cat("Not using snakemake\n")
    cat("Script is intended to be used within snakemake\n")
    stop()
}
