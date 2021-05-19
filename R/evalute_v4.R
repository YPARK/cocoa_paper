#!/usr/bin/env Rscript

hdr <- "3_5_5_0_1"
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
    cols.dt <- .fread("result/" %&% hdr %&% ".mu_cols.gz")
    cols.dt[, c("id", "p1", "p0", "pa", "pf", "rseed") := tstrsplit(V1,split="_")]

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

causal <- .fread(causal.file) %>% .unlist

#' @param mtx.file e.g., "data/3_3_5_3_1.mtx.gz"
#' @param label.file e.g., "data/3_3_5_3_1.label.gz"
#' @param ind.file e.g., "data/3_3_5_3_1.ind.gz"
run.MAST <- function(mtx.file, 
                     label.file,
                     ind.file) {

    library(MAST)

    Y <- Matrix::readMM(mtx.file) %>% as.matrix
    ## Y.int <- apply(Y, 2, function(x) { storage.mode(x) <- 'integer'; x})

    colnames(Y) <- 1:ncol(Y)

    lab <- fread(label.file, header=FALSE) %>% unlist
    ind <- fread(ind.file, header=FALSE) %>% unlist %>%
        str_remove_all("Ind") %>% as.integer

    grp <- lab[ind]

    mm <- nrow(Y)
    nn <- ncol(Y)

    cdat <- data.frame(condition = grp)
    rdat <- data.frame(primerid = 1:mm)

    .dat <- FromMatrix(Y, cdat, rdat, check_sanity=FALSE)
    colData(.dat)$condition <- factor(colData(.dat)$condition)

    .take.dt <- function(.zlm.out) {

        out <- summary(.zlm.out, doLRT='condition1')
        .dt <- out$datatable
        out.dt <-
            .dt[contrast=='condition1' & component=='H',.(primerid, `Pr(>Chisq)`)] %>% 
            dplyr::rename(pv = `Pr(>Chisq)`) %>%
            dplyr::rename(x.col = primerid) %>%
            as.data.table
    }

    zlm(~condition, .dat) %>%
        .take.dt %>% 
        mutate(method = "MAST")
}

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
    .pv <- .stat$pv
    .xx <- .stat$x.col
    .efdr <- function(.cutoff) {
        .qv <- p.adjust(.pv, "fdr")
        .num <- sum(.qv < .cutoff & !(.xx %in% causal))
        .denom <- pmax(sum(.qv < .cutoff), 1)
        return(.num/.denom)
    }

    .edr <- function(.cutoff) {
        .qv <- p.adjust(.pv, "fdr")
        mean(.qv < .cutoff)
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
## run MAST to compare

.mast <- run.MAST(mtx.file,
                  label.file,
                  ind.file)

mast.efdr <- function(.stat, causal) {
    
    .pv <- .stat$pv
    .xx <- .stat$x.col
    .efdr <- function(.cutoff) {
        .qv <- p.adjust(.pv, "fdr")
        .num <- sum(.qv < .cutoff & !(.xx %in% causal))
        .denom <- pmax(sum(.qv < .cutoff), 1)
        return(.num/.denom)
    }
    .edr <- function(.cutoff) {
        .qv <- p.adjust(.pv, "fdr")
        mean(.qv < .cutoff)
    }
    
    data.table(efdr01 = .efdr(.01),
               efdr05 = .efdr(.05),
               efdr10 =  .efdr(.10),
               edr01 = .edr(.01),
               edr05 = .edr(.05),
               edr10 =  .edr(.10))
}

mast.efdr <- .mast[, as.list(mast.efdr(.SD, causal)),
                   by = .(method)]

mast.stat <- .mast[, as.list(power.summary(.SD, causal)),
                   by = .(method)]

################################################################

sum.dt <- .read.named(hdr, ".sum.gz")
mean.dt <- .read.named(hdr, ".mean.gz")
mu.dt <- .read.named(hdr, ".ln_obs_mu.gz")
cocoa.dt <- .read.named(hdr, ".ln_resid_mu.gz")
cf.dt <- .read.named(hdr, ".cf_mu.gz")

################################################################
cocoa.stat <- calc.stat(cocoa.dt) %>% power.summary(causal) %>% mutate(method="cocoa")
cf.stat <- calc.stat(cf.dt) %>% power.summary(causal) %>% mutate(method="cf")
mu.stat <- calc.stat(mu.dt) %>% power.summary(causal) %>% mutate(method="mu")
mean.stat <- calc.stat(mean.dt) %>% power.summary(causal) %>% mutate(method="avg")
sum.stat <- calc.stat(sum.dt) %>% power.summary(causal) %>% mutate(method="tot")

out.1 <- rbind(cocoa.stat,
               cf.stat,
               mu.stat,
               sum.stat,
               mean.stat,
               mast.stat)

cocoa.efdr <- take.efdr(cocoa.dt, causal) %>% mutate(method = "cocoa")
cf.efdr <- take.efdr(cf.dt, causal) %>% mutate(method = "cf")
mu.efdr <- take.efdr(mu.dt, causal) %>% mutate(method = "mu")
mean.efdr <- take.efdr(mean.dt, causal) %>% mutate(method = "avg")
sum.efdr <- take.efdr(sum.dt, causal) %>% mutate(method = "tot")

out.2 <- rbind(cocoa.efdr,
               cf.efdr,
               mu.efdr,
               sum.efdr,
               mean.efdr,
               mast.efdr)

################################################################

out.stat <- out.1 %>%
    left_join(out.2) %>% 
    mutate(sim = hdr) %>%
    select(method, sim, starts_with("power"), starts_with("au"),
           starts_with("efdr"), starts_with("edr")) %>% 
    separate("sim", c("p1","p0","pa","pf","rseed"), sep="[_]")

.fwrite(out.stat, out.file)
