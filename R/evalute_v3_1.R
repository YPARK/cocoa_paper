#!/usr/bin/env Rscript

hdr <- "3_3_8_1"
out.file <- "temp.txt.gz"

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

    data.dt <- .fread("result/" %&% hdr %&% ext, col.names = .cols)
}

lab.dt <- .fread(label.file, col.names="label") %>%
    mutate(id = as.character(1:n()))

calc.stat <- function(.dt) {

    .lab <-
        tibble(id = colnames(.dt)) %>%
        left_join(lab.dt, by = "id") %>%
        select(label) %>%
        unlist

    .mat <- as.matrix(.dt)

    .stat <- apply(t(.mat), 2, function(x) {
        .w <- wilcox.test(x[.lab==1], x[.lab==0])
        data.table(pv = .w$p.value)
    }) %>%
        do.call(what=rbind) %>%
        cbind(x.col = 1:nrow(.mat)) %>%
        na.omit
}

take.effect <- function(.mu.in, .sd.in) {

    .lab <-
        tibble(id = colnames(.mu.in)) %>%
        left_join(lab.dt, by = "id")

    .add.gene <- function(.dt, ...) {
        .dt %>% mutate(x.col = 1:n()) %>%
            melt(id.vars="x.col", variable.name="id", ...)
    }

    .dt <-
        .add.gene(.mu.in, value.name = "mu") %>%
        left_join(.add.gene(.sd.in, value.name = "se")) %>% 
        left_join(.lab)

    .sd <- function(...) sd(..., na.rm=TRUE)

    ## .sd.dt <- .dt[, .(min.sd = sd(mu)/sqrt(.N)), by = .(x.col)]

    min.sd <- .sd(unlist(.mu.in))

    .dt <- .dt %>%
        mutate(w = 1/(se^2 + min.sd^2))

    ret <-
        .dt[,
            .(beta = sum((2 * label - 1) * mu * w) / sum(w),
              se = sqrt(1/sum(w))),
            by = .(x.col)] %>% 
        mutate(z = beta / se) %>% 
        mutate(pv = 2 * pnorm(abs(z), lower.tail = FALSE))

    return(ret)
}

take.power <- function(.stat, causal) {

    nc <- length(causal)
    .out <- .stat[order(pv)]
    .out[, lab := as.integer(x.col %in% causal)]
    .out[, d := 1:nrow(.out)]
    .out <- .out[, .(lab, d)]
    .out[, n1 := cumsum(lab)]
    .out[, n0 := cumsum(1 - lab)]
    .out[, fdr := n0 / d]
    .out[, tdr := n1 / d]
    .out[, prec := n1 / d]
    .out[, power := n1 / nc]
    .out[, rec := n1 / nc]

    .out <- .out[order(prec, decreasing=TRUE), head(.SD, 1), by = .(rec)]

    return(.out)
}

take.summary <- function(.stat) {

    .ret <-
        .stat[order(.stat$rec),
              head(.SD, 1),
              by = .(fdr)]
    
    .auprc <- DescTools::AUC(x = .ret$rec,
                             y = .ret$prec,
                             from = 0, to = 1)

    .ret <-
        .ret %>% arrange(fdr)

    .auroc <- DescTools::AUC(x = .ret$fdr,
                             y = .ret$power,
                             from = 0, to = 1)

    .f1 <- .ret %>%
        mutate(f1 = 2 * (1-fdr) * power / pmax(1-fdr + power, 1e-10)) %>%
        select(f1) %>% 
        slice(which.max(f1)) %>%
        unlist

    .power10 <- 0

    if(nrow(.ret[fdr < .1]) > 1)
        .power10 <- .ret[fdr < .1, .(power)] %>% unlist %>% max

    data.table(f1 = .f1, auprc = .auprc, auroc = .auroc, power10 = .power10)
}

power.summary <- function(.stat, causal) {
    take.power(.stat, causal) %>%
        take.summary
}

sum.dt <- .read.named(hdr, ".sum.gz")
mean.dt <- .read.named(hdr, ".mean.gz")
mu.dt <- .read.named(hdr, ".ln_obs_mu.gz")

cocoa.dt <- .read.named(hdr, ".ln_resid_mu.gz")
cocoa.sd.dt <- .read.named(hdr, ".ln_resid_mu_sd.gz")
cf.dt <- .read.named(hdr, ".cf_mu.gz")

causal <- .fread(causal.file) %>% .unlist

cocoa.stat <- calc.stat(cocoa.dt) %>% power.summary(causal) %>% mutate(method="cocoa")

cocoa.ade.stat <-
    take.effect(.mu.in = cocoa.dt, .sd.in = cocoa.sd.dt) %>%
    power.summary(causal) %>% 
    mutate(method="cocoa.ade")

cf.stat <- calc.stat(cf.dt) %>% power.summary(causal) %>% mutate(method="cf")
mu.stat <- calc.stat(mu.dt) %>% power.summary(causal) %>% mutate(method="mu")
mean.stat <- calc.stat(mean.dt) %>% power.summary(causal) %>% mutate(method="avg")
sum.stat <- calc.stat(sum.dt) %>% power.summary(causal) %>% mutate(method="tot")

################################################################

out.stat <- rbind(cocoa.stat,
                  cocoa.ade.stat,
                  cf.stat,
                  mu.stat,
                  sum.stat,
                  mean.stat) %>%
    mutate(sim = hdr) %>%
    select(method, sim, f1, power10, starts_with("au")) %>% 
    separate("sim", c("p1","p0","pa","rseed"), sep="[_]")

.fwrite(out.stat, out.file)
