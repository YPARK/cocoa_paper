#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 2) {
    q()
}

hdr <- argv[1]      # e.g.,"3_3_0_${rseed}"
out.file <- argv[2] # e.g., "output.txt.gz"

`%&%` <- function(a,b) {
    paste0(a,b)
}

causal.file <- "data/" %&% hdr %&% ".causal.gz"
label.file <- "data/" %&% hdr %&% ".label.gz"

library(data.table)
library(tidyverse)

.fread <- function(...) fread(..., header=FALSE)

.fwrite <- function(...) {
    fwrite(..., row.names=F, sep=" ")
}

.unlist <- function(...) unlist(..., use.names=FALSE)

.read.named <- function(hdr, ext) {
    cols.dt <- .fread("result/" %&% hdr %&% ".mu_cols.gz")
    cols.dt[, c("id", "p1", "p0", "pa", "rseed") := tstrsplit(V1,split="_")]

    .cols <- cols.dt %>%
        select(id) %>%
        unlist %>%
        str_remove("Ind") %>%
        as.integer %>%
        as.character

    data.dt <- .fread("result/" %&% hdr %&% ext,
                      col.names = .cols)
}

lab.dt <- .fread(label.file, col.names="label") %>%
    mutate(id = as.character(1:n()))

calc.stat <- function(.dt, log.trans = FALSE) {

    .lab <-
        tibble(id = colnames(.dt)) %>%
        left_join(lab.dt, by = "id") %>%
        select(label) %>%
        as.matrix

    .mat <- t(.dt) %>% as.matrix

    if(log.trans) {
        .mat <- log(1 + .mat)
    }

    zqtl::calc.qtl.stat(.mat, .lab) %>%
        mutate(q = p.adjust(p.val, "fdr")) %>%
        as.data.table
}

take.power <- function(.stat, causal) {
    nc <- length(causal)
    .out <- .stat[order(.stat$p.val)]
    .out[, lab := as.integer(x.col %in% causal)]
    .out[, d := 1:nrow(.out)]
    .out <- .out[, .(lab, d)]
    .out[, n1 := cumsum(lab)]
    .out[, n0 := cumsum(1 - lab)]
    .out[, fdr := n0 / d]
    .out[, tdr := n1 / d]
    .out[, power := n1 / nc]
    return(.out)
}

take.summary <- function(.stat) {
    .ret <- .stat[order(.stat$power, decreasing=TRUE),
                  head(.SD, 1),
                  by = .(fdr)]

    .ret %>%
        arrange(desc(fdr), power) %>%
        select(fdr, tdr, power) %>%
        mutate(f1 = 2 * (1-fdr) * power / (1-fdr + power))
}

power.summary <- function(.stat, causal) {
    take.power(.stat, causal) %>%
        take.summary
}

sum.dt <- .read.named(hdr, ".sum.gz")
mean.dt <- .read.named(hdr, ".mean.gz")

mu.dt <- .read.named(hdr, ".ln_obs_mu.gz")
cocoa.dt <- .read.named(hdr, ".ln_resid_mu.gz")

causal <- .fread(causal.file) %>% .unlist

cocoa.stat <- calc.stat(cocoa.dt,FALSE) %>% power.summary(causal) %>% mutate(method="cocoa")

mu.stat <- calc.stat(mu.dt,FALSE) %>% power.summary(causal) %>% mutate(method="mu")

sum.stat <- calc.stat(sum.dt,TRUE) %>% power.summary(causal) %>% mutate(method="tot")

mean.stat <- calc.stat(mean.dt,TRUE) %>% power.summary(causal) %>% mutate(method="avg")

out.stat <- rbind(cocoa.stat,
                  mu.stat,
                  mean.stat,
                  sum.stat) %>%
    mutate(sim = hdr) %>%
    separate("sim", c("p1","p0","pa","rseed"), sep="[_]")

.fwrite(out.stat, out.file)
