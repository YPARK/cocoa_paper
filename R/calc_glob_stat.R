#!/usr/bin/env Rscript

ln.mu.file <- "result/cocoa/Microglia_1.ln_resid_mu.gz"
col.file <- "result/cocoa/Microglia_1.mu_cols.gz"
sum.file  <- "result/aggregate/Microglia.sum.gz"
mean.file  <- "result/aggregate/Microglia.mean.gz"
sum.col.file  <- "result/aggregate/Microglia.mu_cols.gz"
row.file <- "data/brain_2018-05-03/features.tsv.gz"
pheno.file <- "result/phenotyped.txt.gz"

argv <- commandArgs(trailingOnly = TRUE)

ln.mu.file   <- argv[1] # "Microglia.ln_resid_mu.gz"
col.file     <- argv[2] # "Microglia.mu_cols.gz"
sum.file     <- argv[3] # "Microglia.sum.gz"
sum.col.file <- argv[4] # "Microglia.mu_cols.gz"
row.file     <- argv[5] # "genes.rows.gz"
pheno.file   <- argv[6] # "clinical.csv"
out.file     <- argv[7]

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

run.t.test <- function(gg, xx, yy) {

    x.obs <- apply(t(!is.na(xx)), 2, mean)
    y.obs <- apply(t(!is.na(xx)), 2, mean)

    .valid <- which(x.obs > .5 & y.obs > .5)

    .t.test <- lapply(.valid, function(j) t.test(yy[j, ], xx[j, ]))
    
    pv <- sapply(.t.test, function(x) x$p.value)
    tt <- sapply(.t.test, function(x) x$statistic)
    se <- sapply(.t.test, function(x) x$stderr)

    gg[.valid, ] %>%
        mutate(pv = pv, t = tt, se = se)

}

.pheno <- fread(pheno.file, header=TRUE) %>%
    select(-TAG) %>%
    distinct

.read.cocoa <- function(...){

    rows <- fread(row.file, header = FALSE, col.names = "gene")
    cols <- fread(col.file, header = FALSE, col.names = "sample")
    cols[, c("projid", "celltype", "pheno.col") := tstrsplit(sample, "_")]
    cols[, projid := as.integer(projid)]

    return(.fread.melt(..., rows = rows$gene, cols = cols$projid))
}


.read.sum <- function(...){

    rows <- fread(row.file, header = FALSE, col.names = "gene")
    cols <- fread(sum.col.file, header = FALSE, col.names = "sample")
    cols[, c("projid", "celltype") := tstrsplit(sample, "_")]
    cols[, projid := as.integer(projid)]

    return(.fread.melt(..., rows = rows$gene, cols = cols$projid))
}



stat.dt <- 
    .read.cocoa(ln.mu.file, "ln.mu") %>%
    left_join(.read.sum(sum.file, "tot"), by = c("gene", "projid")) %>%
    left_join(.read.sum(mean.file, "avg"), by = c("gene", "projid")) %>%
    left_join(.pheno, by = "projid") %>%
    as.data.table

.xx <- stat.dt %>%
    filter(gene == "APOE") %>%
    arrange(ln.mu)

.xx


t.test(.xx[pathoAD==1]$ln.mu, .xx[pathoAD==0]$ln.mu)

ggplot(.xx, aes(x=as.factor(pathoAD), y = ln.mu)) +
    geom_boxplot()

t.test(.xx[pathoAD==1]$avg, .xx[pathoAD==0]$avg)

ggplot(.xx, aes(x=as.factor(pathoAD), y = log(avg))) +
    geom_boxplot()




.glob.dt <- dcast(stat.dt[tot > 0], hgnc + ensg ~ pathoAD + projid, value.var = "ln.mu")

gg <- .glob.dt %>% select(ensg, hgnc) %>% as_tibble
xx <- .glob.dt %>% select(starts_with("0_")) %>% as.matrix
yy <- .glob.dt %>% select(starts_with("1_")) %>% as.matrix

patho.test.dt <- run.t.test(gg, xx, yy)

write_tsv(patho.test.dt, file = out.file)
