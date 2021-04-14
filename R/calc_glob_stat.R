#!/usr/bin/env Rscript

ln.mu.file <- "result/cocoa/1.ln_resid_mu.gz"
ln.mu.sd.file <- "result/cocoa/1.ln_resid_mu_sd.gz"
cf.mu.file <- "result/cocoa/1.cf_mu.gz"
cf.mu.sd.file <- "result/cocoa/1.cf_mu_sd.gz"
col.file <- "result/cocoa/1.mu_cols.gz"
sum.file  <- "result/aggregate/1.sum.gz"
mean.file  <- "result/aggregate/1.mean.gz"
sum.col.file  <- "result/aggregate/1.mu_cols.gz"
row.file <- "data/brain_2018-05-03/features.tsv.gz"
pheno.file <- "result/phenotyped.txt.gz"
pheno.col <- 1

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 12) {
    cat("Missing arguments: ", paste(argv, collapse=", "))
    q()
}

ln.mu.file    <- argv[1] # "Microglia.ln_resid_mu.gz"
ln.mu.sd.file <- argv[2] # "Microglia.ln_resid_mu_sd.gz"
cf.mu.file    <- argv[3] # "Microglia.cf_mu.gz"
cf.mu.sd.file <- argv[4] # "Microglia.cf_mu_sd.gz"
col.file      <- argv[5] # "Microglia.mu_cols.gz"
sum.file      <- argv[6] # "Microglia.sum.gz"
mean.file     <- argv[7] # "Microglia.sum.gz"
sum.col.file  <- argv[8] # "Microglia.mu_cols.gz"
row.file      <- argv[9] # "genes.rows.gz"
pheno.file    <- argv[10] # "clinical.csv"
pheno.col     <- as.integer(argv[11])
out.file      <- argv[12]

tot.cutoff <- 10

################################################################

library(tidyverse)
library(data.table)

`%&%` <- function(a,b) paste0(a,b)

################################################################

.fread.melt <- function(x, .name, rows, cols = NULL, col.file = NULL) {

    if(is.null(cols)) {
        cols <- fread(col.file, header = FALSE, col.names = "sample") %>%
            unlist
    }

    .ret <- fread(x,
                  header = FALSE,
                  col.names = as.character(cols),
                  colClasses = "double")

    .ret[, gene := rows]

    ret <-
        melt.data.table(.ret,
                        id.vars = "gene",
                        variable.name = "sample",
                        value.name = .name)

    ret[, c("projid", "celltype") := tstrsplit(sample, "_")]
    ret[, projid := as.integer(projid)]

    return(ret)
}

run.test <- function(gg, xx, yy) {

    .wilcox.test <- lapply(1:nrow(xx), function(j) { wilcox.test(yy[j, ], xx[j, ]) })

    pv <- sapply(.wilcox.test, function(x) x$p.value)

    ret <- gg %>%
        mutate(wilcox.pv = pv) %>%
        as.data.table
}

test.all.genes <- function(.var, .stat.dt) {

    .glob.dt <- dcast(.stat.dt,
                      gene + celltype ~ pheno + projid, value.var = .var)

    gg <- .glob.dt %>% mutate(gene = gene %&% "_" %&% celltype) %>% select(gene) %>% as_tibble
    xx <- .glob.dt %>% select(starts_with("0_")) %>% as.matrix
    yy <- .glob.dt %>% select(starts_with("1_")) %>% as.matrix

    ret <- run.test(gg, xx, yy)
    ret[, c("gene", "celltype") := tstrsplit(gene, "_")]
    return(ret)
}

.read.cocoa <- function(...){

    rows <- fread(row.file, header = FALSE, col.names = "gene")
    cols <- fread(col.file, header = FALSE, col.names = "sample")
    cols[, c("projid", "celltype") := tstrsplit(sample, "_")]
    cols[, col := projid %&% "_" %&% celltype]

    return(.fread.melt(..., rows = rows$gene, cols = cols$col))
}

.read.sum <- function(...){

    rows <- fread(row.file, header = FALSE, col.names = "gene")
    cols <- fread(sum.col.file, header = FALSE, col.names = "sample")
    cols[, c("projid", "celltype") := tstrsplit(sample, "_")]
    cols[, col := projid %&% "_" %&% celltype]

    return(.fread.melt(..., rows = rows$gene, cols = cols$col))
}

.take.effect <- function(.dt, .cutoff = 1) {

    .dt <- copy(as.data.table(.dt))

    ## adjust potential bias
    .dt[tot >= .cutoff, b := mean(beta), by = .(gene, celltype)]
    .dt[, beta := beta - b]

    ## take into account uncertainty across individuals
    .sd <- .dt[tot >= .cutoff,
               .(min.sd = sd(beta) / sqrt(.N)),
               by = .(gene, celltype)]

    .temp <- copy(.dt) %>%
        filter(tot >= .cutoff) %>%
        left_join(.sd) %>%
        mutate(w = 1/(beta.se^2 + min.sd^2)) %>%
        as.data.table

    .temp1 <- .temp[,
                    .(ADE.beta = sum((2 * pheno - 1) * beta * w) / sum(w),
                      ADE.se = sqrt(1/sum(w))),
                    by = .(gene, celltype)]

    .temp2 <- .temp %>% filter(pheno == 1) %>% as.data.table

    .temp2 <- .temp2[,
                    .(ADD.beta = sum(beta * w) / sum(w),
                      ADD.se = sqrt(1/sum(w))),
                    by = .(gene, celltype)]

    .temp3 <- .temp %>% filter(pheno == 0) %>% as.data.table

    .temp3 <- .temp3[,
                    .(ADC.beta = -sum(beta * w) / sum(w),
                      ADC.se = sqrt(1/sum(w))),
                    by = .(gene, celltype)]

    .ret <- .temp1 %>%
        left_join(.temp2) %>%
        left_join(.temp3) %>%
        filter(!is.na(ADC.beta) | !is.na(ADD.beta)) %>%
        as.data.table

    return(.ret)
}

.pheno <-
    fread(pheno.file, header=TRUE, na.strings="-9") %>%
    select(-TAG) %>%
    distinct %>%
    as.data.frame

.pheno.name <- names(.pheno)[1 + pheno.col]

.pheno <-
    .pheno[, c(1, pheno.col + 1)] %>%
    (function(x) { colnames(x) = c("projid", "pheno"); x }) %>%
    mutate(projid = as.integer(projid))

######################
## read data matrix ##
######################

cf.dt <- 
    .read.cocoa(cf.mu.file, "cf") %>% 
    left_join(.read.cocoa(cf.mu.sd.file, "cf.sd")) %>%
    mutate(ln.cf = log(cf), ln.cf.sd = cf.sd/cf)

stat.dt <-
    .read.cocoa(ln.mu.file, "ln.mu") %>%
    left_join(.read.cocoa(ln.mu.sd.file, "ln.mu.sd")) %>%
    left_join(cf.dt) %>% 
    left_join(.read.sum(sum.file, "tot")) %>%
    left_join(.read.sum(mean.file, "avg")) %>%
    left_join(.pheno) %>%
    filter(!is.na(celltype)) %>% 
    as.data.table

if(length(unique(stat.dt$pheno)) > 2) {
    .max.val <- max(stat.dt$pheno)
    stat.dt <- stat.dt[pheno != 1]
    stat.dt[pheno == .max.val, pheno := 1]
}

tot.1.dt <- stat.dt[pheno == 1,
                  .(tot = sum(tot),
                    sig.cocoa = sd(ln.mu, na.rm=TRUE),
                    sig.tot = sd(log(1 + tot), na.rm=TRUE)),
                  by = .(gene, celltype)]

tot.0.dt <- stat.dt[pheno == 0,
                  .(tot = sum(tot),
                    sig.cocoa = sd(ln.mu, na.rm=TRUE),
                    sig.tot = sd(log(1 + tot), na.rm=TRUE)),
                  by = .(gene, celltype)]

tot.dt <- left_join(tot.1.dt, tot.0.dt,
                    by = c("gene", "celltype"),
                    suffix = c(".1", ".0")) %>%
    as.data.table

.ret.0 <- stat.dt %>%
    mutate(beta = ln.mu, beta.se = ln.mu.sd) %>%
    .take.effect(.cutoff = tot.cutoff)

.ret.null <- stat.dt %>%
    mutate(beta = ln.cf, beta.se = ln.cf.sd) %>%
    .take.effect(.cutoff = tot.cutoff)

.ret.1 <- test.all.genes("ln.mu", stat.dt) %>% rename(pv.cocoa = wilcox.pv)
.ret.2 <- test.all.genes("avg", stat.dt) %>% rename(pv.avg = wilcox.pv)
.ret.3 <- test.all.genes("tot", stat.dt) %>% rename(pv.tot = wilcox.pv)
.ret.4 <- test.all.genes("ln.cf", stat.dt) %>% rename(pv.cf = wilcox.pv)

ret <-
    .ret.0 %>%
    left_join(.ret.null, by=c("gene","celltype"), suffix=c("",".cf")) %>% 
    left_join(.ret.1) %>% 
    left_join(.ret.2) %>% 
    left_join(.ret.3) %>%
    left_join(.ret.4) %>%
    mutate(pheno = .pheno.name) %>%
    left_join(tot.dt) %>%
    as.data.table

fwrite(ret, file = out.file, sep = "\t", row.names = FALSE, col.names = TRUE)

## save top 20 genes per cell type to visualize
.ret <- ret[sign(ADD.beta) == sign(ADC.beta) & pv.cocoa < 1e-2]

top.genes <-
    .ret[order(.ret$pv.cocoa), head(.SD, 20), by = .(celltype)] %>% 
    select(gene) %>% 
    unique %>% 
    left_join(stat.dt, by = "gene") %>%
    rename(label = pheno) %>%
    mutate(pheno = .pheno.name) %>%
    as.data.table

out.data.file <- str_replace(out.file, ".gz", "_data.gz")

fwrite(top.genes, file = out.data.file, sep = "\t", row.names=FALSE, col.names=TRUE)
