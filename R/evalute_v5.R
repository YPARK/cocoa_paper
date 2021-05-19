#!/usr/bin/env Rscript

hdr <- "3_3_5_0_2"
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
mtx.file <- "data/" %&% hdr %&% ".mtx.gz"
label.file <- "data/" %&% hdr %&% ".label.gz"
ind.file <- "data/" %&% hdr %&% ".ind.gz"

library(data.table)
library(tidyverse)

.fread <- function(...) fread(..., header=FALSE)

.fwrite <- function(...) {
    fwrite(..., row.names=F, sep=" ", na="NA", quote=FALSE)
}

.unlist <- function(...) unlist(..., use.names=FALSE)

.read.named <- function(hdr, ext) {
    cols.dt <- .fread(hdr %&% ".mu_cols.gz")
    cols.dt[, c("id", "p1", "p0", "pa", "pf", "rseed") := tstrsplit(V1,split="_")]

    .cols <- cols.dt %>%
        select(id) %>%
        unlist %>%
        str_remove("Ind") %>%
        as.integer %>%
        as.character

    data.dt <- .fread(hdr %&% ext, col.names = .cols)
}

lab.dt <- .fread(label.file, col.names="label") %>%
    mutate(id = as.character(1:n()))

causal <- .fread(causal.file) %>% .unlist

calc.stat <- function(.dt, .cutoff = -4, .perm = FALSE) {

    .lab <-
        tibble(id = colnames(.dt)) %>%
        left_join(lab.dt, by = "id") %>%
        select(label) %>%
        unlist

    .mat <- as.matrix(.dt)

    .mat[.mat < .cutoff] <- NA

    .stat <- apply(t(.mat), 2, function(x) {
        if(.perm) { .lab <- sample(.lab) }
        .w <- wilcox.test(x[.lab==1], x[.lab==0])
        data.table(pv = .w$p.value)
    }) %>%
        do.call(what=rbind) %>%
        cbind(x.col = 1:nrow(.mat)) %>%
        na.omit
}

take.efdr <- function(.dt, causal) {

    .stat <- calc.stat(.dt)
    .stat.0 <- calc.stat(.dt, .perm = TRUE)
    .pv <- qvalue::empPvals(-log10(.stat$pv), -log10(.stat.0$pv))
    .xx <- .stat$x.col
    .efdr <- function(.cutoff) {
        .qv <- qvalue::qvalue(.pv, fdr.level = .cutoff)
        .num <- sum(.qv$qvalues < .cutoff & !(.xx %in% causal))
        .denom <- pmax(sum(.qv$qvalues < .cutoff), 1)
        return(.num/.denom)
    }

    .edr <- function(.cutoff) {
        .qv <- qvalue::qvalue(.pv, fdr.level = .cutoff)
        mean(.qv$qvalues < .cutoff)
    }

    data.table(efdr01 = .efdr(.01),
               efdr05 = .efdr(.05),
               efdr10 =  .efdr(.10),
               edr01 = .edr(.01),
               edr05 = .edr(.05),
               edr10 =  .edr(.10))
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
    
    .auprc <- DescTools::AUC(x = c(0, .ret$rec, 1),
                             y = c(0, .ret$prec, 0),
                             from = 0, to = 1)

    .ret <-
        .ret %>% arrange(fdr)

    .auroc <- DescTools::AUC(x = c(0, .ret$fdr, 1),
                             y = c(0, .ret$power, 1),
                             from = 0, to = 1)

    .power10 <- 0

    if(nrow(.ret[fdr < .1]) > 1)
        .power10 <- .ret[fdr < .1, .(power)] %>% unlist %>% max

    .power5 <- 0

    if(nrow(.ret[fdr < .05]) > 1)
        .power5 <- .ret[fdr < .05, .(power)] %>% unlist %>% max

    .power1 <- 0

    if(nrow(.ret[fdr < .01]) > 1)
        .power1 <- .ret[fdr < .01, .(power)] %>% unlist %>% max

    data.table(auprc = .auprc, auroc = .auroc, power01 = .power1, power05 = .power5, power10 = .power10)
}

power.summary <- function(.stat, causal) {
    take.power(.stat, causal) %>%
        take.summary
}

################################################################

dat1 <- .read.named("1/" %&% hdr, ".ln_resid_mu.gz")
dat10 <- .read.named("10/" %&% hdr, ".ln_resid_mu.gz")
dat50 <- .read.named("50/" %&% hdr, ".ln_resid_mu.gz")
dat100 <- .read.named("100/" %&% hdr, ".ln_resid_mu.gz")
dat200 <- .read.named("200/" %&% hdr, ".ln_resid_mu.gz")

################################################################
stat.1 <- calc.stat(dat1) %>% power.summary(causal) %>% mutate(knn=1)
stat.10 <- calc.stat(dat10) %>% power.summary(causal) %>% mutate(knn=10)
stat.50 <- calc.stat(dat50) %>% power.summary(causal) %>% mutate(knn=50)
stat.100 <- calc.stat(dat100) %>% power.summary(causal) %>% mutate(knn=100)
stat.200 <- calc.stat(dat200) %>% power.summary(causal) %>% mutate(knn=200)

out.1 <- rbind(stat.1,
               stat.10,
               stat.50,
               stat.100,
               stat.200)

efdr.1 <- take.efdr(dat1, causal) %>% mutate(knn = 1)
efdr.10 <- take.efdr(dat10, causal) %>% mutate(knn = 10)
efdr.50 <- take.efdr(dat50, causal) %>% mutate(knn = 50)
efdr.100 <- take.efdr(dat100, causal) %>% mutate(knn = 100)
efdr.200 <- take.efdr(dat200, causal) %>% mutate(knn = 200)

out.2 <- rbind(efdr.1,
               efdr.10,
               efdr.50,
               efdr.100,
               efdr.200)

################################################################

out.stat <-
    out.1 %>%
    left_join(out.2) %>% 
    mutate(sim = hdr) %>%
    select(knn, sim, starts_with("power"), starts_with("au"),
           starts_with("efdr"), starts_with("edr")) %>% 
    separate("sim", c("p1","p0","pa","pf","rseed"), sep="[_]")

.fwrite(out.stat, out.file)
