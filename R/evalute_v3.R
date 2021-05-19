#!/usr/bin/env Rscript

hdr <- "3_3_3_1"
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

    data.table(auprc = .auprc, auroc = .auroc, power10 = .power10)
}

power.summary <- function(.stat, causal) {
    take.power(.stat, causal) %>%
        take.summary
}

## adjust PCs
.adj.data <- function(.dt, .svd, kk, take.log=TRUE) {

    if(take.log) {
        .mat <- log(1 + t(as.matrix(.dt)))
    } else {
        .mat <- t(as.matrix(.dt))
    }

    .fun <- function(x) {
        .lm <- lm(x ~ .svd$v[, kk, drop=FALSE])
        residuals(.lm)
    }

    .temp <- apply(.mat, 2, .fun)
    
    ret <- as.data.table(t(.temp))
    colnames(ret) <- colnames(sum.dt)
    return(ret)
}

sum.dt <- .read.named(hdr, ".sum.gz")
causal <- .fread(causal.file) %>% .unlist

mean.dt <- .read.named(hdr, ".mean.gz")
mu.dt <- .read.named(hdr, ".ln_obs_mu.gz")
cocoa.dt <- .read.named(hdr, ".ln_resid_mu.gz")
cf.dt <- .read.named(hdr, ".cf_mu.gz")

################################################################
## RUV to find factors that are less to be causal
.svd <- rsvd::rsvd(log(1 + as.matrix(sum.dt)), k=10)
.cor.pc <- lapply(1:10, function(k) cor.test(.svd$v[,k], lab.dt$label))
pv.cutoff <- .05
.nc.pc <- which(sapply(.cor.pc, function(x) x$p.value > pv.cutoff))

sum.pc.safe.dt <- .adj.data(sum.dt, .svd, .nc.pc)
sum.pc.full.dt <- .adj.data(sum.dt, .svd, 1:10)

cocoa.pc.safe.dt <- .adj.data(cocoa.dt, .svd, .nc.pc, take.log=FALSE)
cocoa.pc.full.dt <- .adj.data(cocoa.dt, .svd, 1:10, take.log=FALSE)

################################################################
cocoa.stat <- calc.stat(cocoa.dt) %>% power.summary(causal) %>% mutate(method="cocoa")
cf.stat <- calc.stat(cf.dt) %>% power.summary(causal) %>% mutate(method="cf")
mu.stat <- calc.stat(mu.dt) %>% power.summary(causal) %>% mutate(method="mu")
mean.stat <- calc.stat(mean.dt) %>% power.summary(causal) %>% mutate(method="avg")
sum.stat <- calc.stat(sum.dt) %>% power.summary(causal) %>% mutate(method="tot")

sum.pc.safe.stat <- calc.stat(sum.pc.safe.dt) %>%
    power.summary(causal) %>% mutate(method="tot.pc.safe")

sum.pc.full.stat <- calc.stat(sum.pc.full.dt) %>%
    power.summary(causal) %>% mutate(method="tot.pc.full")

cocoa.pc.safe.stat <- calc.stat(cocoa.pc.safe.dt) %>%
    power.summary(causal) %>% mutate(method="cocoa.pc.safe")

cocoa.pc.full.stat <- calc.stat(cocoa.pc.full.dt) %>%
    power.summary(causal) %>% mutate(method="cocoa.pc.full")

################################################################

clean.mu.dt <- .fread("data/" %&% hdr %&% ".mu-clean.gz")
covar.dt <- .fread("data/" %&% hdr %&% ".covar.gz")

mm.true <- clean.mu.dt %>% as.matrix

.y <- t(log(1 + as.matrix(sum.dt)))
.x <- t(as.matrix(covar.dt))
.lm <- lm(.y ~ .x)
cf.true <- t(predict(.lm)) %>% as.matrix

mm.cocoa <- cocoa.dt %>% as.matrix
mm.cf <- cf.dt %>% as.matrix
mm.mu <- mu.dt %>% as.matrix
mm.mean <- mean.dt %>% as.matrix
mm.sum <- sum.dt %>% as.matrix

mm.sum.pc.safe <- sum.pc.safe.dt %>% as.matrix
mm.sum.pc.full <- sum.pc.full.dt %>% as.matrix

mm.cocoa.pc.safe <- cocoa.pc.safe.dt %>% as.matrix
mm.cocoa.pc.full <- cocoa.pc.full.dt %>% as.matrix

causal.cor <- function(.hat, .true) {
    sapply(causal,
           function(r) {
               cor(.hat[r, ], .true[r, ], method="spearman")
           }) %>%
        mean
}           

noncausal.cor <- function(.hat, .true) {
    sapply(setdiff(1:nrow(.hat), causal),
           function(r) {
               cor(.hat[r, ], .true[r, ], method="spearman")
           }) %>%
        mean
}           

.methods <- c("cocoa", "cf", "mu", "tot", "avg",
              "tot.pc.safe", "tot.pc.full",
              "cocoa.pc.safe", "cocoa.pc.full")

.mat.list <- list(mm.cocoa, mm.cf, mm.mu, mm.sum, mm.mean,
                  mm.sum.pc.safe, mm.sum.pc.full,
                  mm.cocoa.pc.safe, mm.cocoa.pc.full)

rr <- sapply(.mat.list, causal.cor, .true = mm.true)
r0 <- sapply(.mat.list, noncausal.cor, .true = mm.true)
s0 <- sapply(.mat.list, noncausal.cor, .true = cf.true)

.dt <- data.table(method = .methods, r = rr, r0 = r0, s0 = s0)

################################################################

out.stat <- rbind(cocoa.stat,
                  cf.stat,
                  mu.stat,
                  sum.stat,
                  mean.stat,
                  sum.pc.safe.stat,
                  sum.pc.full.stat,
                  cocoa.pc.safe.stat,
                  cocoa.pc.full.stat) %>%
    left_join(.dt) %>% 
    mutate(sim = hdr) %>%
    select(method, sim, power10, starts_with("au"), r, r0, s0) %>% 
    separate("sim", c("p1","p0","pa","rseed"), sep="[_]")

.fwrite(out.stat, out.file)
