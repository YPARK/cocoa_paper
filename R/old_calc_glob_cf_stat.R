#!/usr/bin/env Rscript

mu.file <- "result/cocoa/Endo_1.cf_mu.gz"
col.file <- "result/cocoa/Endo_1.mu_cols.gz"
sum.file  <- "result/aggregate/Endo_1.sum.gz"
mean.file  <- "result/aggregate/Endo_1.mean.gz"
sum.col.file  <- "result/aggregate/Endo_1.mu_cols.gz"
row.file <- "data/brain_2018-05-03/features.tsv.gz"
pheno.file <- "result/phenotyped.txt.gz"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 8) {
    q()
}

mu.file      <- argv[1]
col.file     <- argv[2]
sum.file     <- argv[3]
mean.file    <- argv[4]
sum.col.file <- argv[5]
row.file     <- argv[6]
pheno.file   <- argv[7]
out.file     <- argv[8]

################################################################

library(tidyverse)
library(data.table)

################################################################

.fread.melt <- function(x, .name, rows, cols = NULL, col.file = NULL) {

    if(is.null(cols)) {
        cols <- fread(col.file, header = FALSE, col.names = "sample")
        cols[, c("projid", "celltype") := tstrsplit(sample, "_")]
        cols <- cols$projid %>% unlist
    }

    .ret <- fread(x, header = FALSE, col.names = as.character(cols))
    .ret[, gene := rows]
    ret <- melt.data.table(.ret, id.vars = "gene",
                           variable.name = "projid",
                           value.name = .name)
    ret[, projid := as.character(projid)]
    ret[, projid := as.integer(projid)]
    return(ret)
}

run.test <- function(gg, xx, yy, max.mean.na = .3) {

    x.missing <- apply(t(is.na(xx)), 2, mean)
    y.missing <- apply(t(is.na(yy)), 2, mean)

    .valid <- which(x.missing <= max.mean.na & y.missing <= max.mean.na)

    .wilcox.test <- lapply(.valid, function(j) { wilcox.test(yy[j, ], xx[j, ]) })

    pv <- sapply(.wilcox.test, function(x) x$p.value)
    tt <- sapply(.wilcox.test, function(x) x$statistic)

    ret <- gg[.valid, ] %>%
        mutate(wilcox.pv = pv, wilcox.stat = tt) %>% as.data.table

    .t.test <- lapply(.valid, function(j) { t.test(log(1 + yy[j, ]), log(1 + xx[j, ])) })

    pv <- sapply(.t.test, function(x) x$p.value)
    tt <- sapply(.t.test, function(x) x$statistic)
    se <- sapply(.t.test, function(x) x$stderr)

    ret <- ret %>%
        mutate(missing.1 = y.missing[.valid]) %>%
        mutate(missing.0 = x.missing[.valid]) %>%
        mutate(t.pv = pv, t.stat = tt, t.se = se) %>%
        as.data.table
}

.read.mu <- function(...){

    rows <- fread(row.file, header = FALSE, col.names = "gene")
    cols <- fread(col.file, header = FALSE, col.names = "sample")
    cols[, c("projid", "celltype", "pheno.col") := tstrsplit(sample, "_")]
    cols[, projid := as.integer(projid)]

    return(.fread.melt(..., rows = rows$gene, cols = cols$projid))
}

.read.sum <- function(...){

    rows <- fread(row.file, header = FALSE, col.names = "gene")
    cols <- fread(sum.col.file, header = FALSE, col.names = "sample")
    cols[, c("projid", "celltype", "pheno.col") := tstrsplit(sample, "_")]
    cols[, projid := as.integer(projid)]

    return(.fread.melt(..., rows = rows$gene, cols = cols$projid))
}

.pheno.col <-
    basename(mu.file) %>%
    str_split(pattern="[.]", simplify=TRUE) %>%
    (function(x) x[1]) %>%
    str_split(pattern="[_]", simplify=TRUE) %>%
    (function(x) x[2]) %>%
    as.integer

.pheno <- fread(pheno.file, header=TRUE, na.strings="-9") %>%
    select(-TAG) %>%
    distinct %>%
    as.data.frame

.pheno.name <- names(.pheno)[1 + .pheno.col]

.pheno <- .pheno[, c(1, .pheno.col + 1)] %>%
    (function(x) { colnames(x) <- c("projid", "pheno"); x })

stat.dt <-
    .read.mu(mu.file, "cf.mu") %>%
    left_join(.read.sum(sum.file, "tot"), by = c("gene", "projid")) %>%
    left_join(.read.sum(mean.file, "avg"), by = c("gene", "projid")) %>%
    left_join(.pheno) %>%
    as.data.table

test.all.genes <- function(.var, tot.cutoff = 1) {

    .glob.dt <- dcast(stat.dt[tot >= tot.cutoff],
                      gene ~ pheno + projid, value.var = .var)

    gg <- .glob.dt %>% select(gene) %>% as_tibble
    xx <- .glob.dt %>% select(starts_with("0_")) %>% as.matrix
    yy <- .glob.dt %>% select(starts_with("1_")) %>% as.matrix

    run.test(gg, xx, yy)
}

.ret.1 <- test.all.genes("cf.mu") %>% mutate(method = "cf")
.ret.2 <- test.all.genes("avg") %>% mutate(method = "avg")
.ret.3 <- test.all.genes("tot") %>% mutate(method = "tot")

ret <-
    rbind(.ret.1, .ret.2, .ret.3) %>%
    mutate(pheno = .pheno.name) %>%
    as.data.table

fwrite(ret, file = out.file, sep = "\t", row.names = FALSE, col.names = TRUE)
