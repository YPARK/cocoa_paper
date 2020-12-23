argv <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(data.table)
library(patchwork)
source("Util.R")

ln.mu.file <- "result/step4/cfa/Microglia.ln_resid_mu.gz"
ln.sd.file <- "result/step4/cfa/Microglia.ln_resid_mu_sd.gz"
col.file <- "result/step4/cfa/Microglia.mu_cols.gz"

sum.file <- "result/step4/aggregate/Microglia.sum.gz"
sum.col.file <- "result/step4/aggregate/Microglia.mu_cols.gz"

row.file <- "result/step1/merged.rows.gz"
pheno.file <- "data/ROSMAP_clinical.csv"




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
}

rows <- fread(row.file, header = FALSE, col.names = "gene")
rows[, c("ensg", "hgnc") := tstrsplit(gene, "_")]

## Reading the aggregate results
tot.dt <- .fread.melt(sum.file, "tot", rows$gene, col.file = sum.col.file)

## Reading the counterfactually-adjusted results
cols <- fread(col.file, header = FALSE, col.names = "sample")
cols[, c("projid", "celltype") := tstrsplit(sample, "_")]
cols[, projid := as.integer(projid)]

.read.cfa <- function(...)
    .fread.melt(..., rows = rows$gene, cols = cols$projid)

## Join with phenotype information
.age.code <- function(x) {
    if(nchar(x) == 0) return(NA)
    if(x == "90+") return(2);
    if(as.numeric(x) > 80) return(1);
    return(0)
}

.apoe.code <- function(x) {
    if(is.na(x)) return(NA) 
    if(x == 44) return(2)
    if(x %in% c(24, 34)) return(1)
    return(0) 
}

.ad.dx <- function(b) {
    ret <- rep(NA, length(b))
    ret[!is.na(b)] <- 0
    ret[b >= 3] <- 1
    return(ret)
}

.pheno <- fread(pheno.file, header=TRUE)
.pheno[, age.death := sapply(age_death, .age.code)]
.pheno[, age.dx := sapply(age_first_ad_dx, .age.code)]
.pheno[, apoe.e4 := sapply(apoe_genotype, .apoe.code)]
.pheno[, pathoAD := sapply(braaksc, .ad.dx)]

stat.dt <- tot.dt %>% 
    merge(.read.cfa(ln.mu.file, "ln.mu")) %>% 
    merge(.read.cfa(ln.sd.file, "ln.sd")) %>% 
    merge(rows, by = "gene") %>%
    left_join(.pheno, by = "projid")

################################################################

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

.plot.gene.view <- function(g, tot.stat.dt) {

    .plot <- function(.dt, .dt.n, .x.lab, .title,
                      .aes.1, .aes.2, .aes.3,
                      is.log.y = TRUE,
                      fill.colours = c("gray80", "#2B8CBE")) {

        p1 <-
            ggplot(.dt.n, .aes.1) +
            ggtitle(.title) +
            theme_classic() +
            theme(axis.title = element_text(size=6)) +
            theme(axis.text = element_text(size=6)) +
            theme(axis.title.x = element_blank()) +
            theme(axis.ticks = element_blank()) +
            theme(axis.text = element_blank()) +
            geom_bar(stat = "identity", fill = "gray60", width = .1) +
            ggrepel::geom_text_repel(aes(label = n), size = 2) +
            ylab("#ind")

        p2 <-
            ggplot(.dt, .aes.2) +
            theme_classic() +
            theme(axis.title = element_text(size=6)) +
            theme(axis.text = element_text(size=6)) +
            theme(axis.title.x = element_blank()) +
            geom_boxplot(width = .2, outlier.size=0, outlier.stroke=0) + 
            scale_fill_manual(values = fill.colours, guide = FALSE) +
            ylab("adjusted")
        
        .y.2 <- .dt %>% select(as_label(.aes.2$y)) %>% unlist
        .y.2.lim <- boxplot.stats(.y.2)$stats[c(1,5)]
        .y.scale <- scale_y_continuous(limits = .y.2.lim*1.01)

        .lab <- function(x) format(exp(x), digits=2, scientific=TRUE)

        if(is.log.y)
            .y.scale <- scale_y_continuous(labels = .lab, limits = .y.2.lim*1.01)

        p2 <- p2 + .y.scale

        .tot.med <- .dt %>% select(as_label(.aes.3$y)) %>% unlist %>% median

        p3 <-
            ggplot(.dt, .aes.3) +
            theme_classic() +
            theme(axis.title = element_text(size=6)) +
            theme(axis.text = element_text(size=6)) +
            geom_boxplot(width = .2, outlier.size=0, outlier.stroke=0) + 
            scale_fill_manual(values = fill.colours, guide = FALSE) +
            ylab("unadjusted") +
            xlab(.x.lab)
        
        .y.3 <- .dt %>% select(as_label(.aes.3$y)) %>% unlist
        .y.3.lim <- boxplot.stats(.y.3)$stats[c(1,5)]
        .y.scale <- scale_y_continuous(limits = .y.3.lim*1.01)

        if(is.log.y)
            .y.scale <- scale_y_continuous(labels = .lab, limits = .y.3.lim*1.01)

        p3 <- p3 + .y.scale

        plt <-
            (p1/p2/p3) + plot_layout(heights = c(1, 6, 3))

        return(plt)
    }

    .dt <- tot.stat.dt[hgnc == g & tot > 0]

    .dt[, ln.tot := log(1 + tot)]
    .dt[, apoe.e4 := factor(apoe.e4, c(0, 1, 2), c("none", "E4", "E4E4"))]
    .dt[, msex := factor(msex, c(0, 1), c("Female", "Male"))]
    .dt[, age.death := factor(.dt$age.death, c(0, 1, 2), c("[, 80)","[80,90]","[90, )"))]
    .dt[, educ := factor(round(as.integer(.dt$educ)/10),
                         c(0, 1, 2, 3),
                         c("[, 5)", "[5, 15)", "[15, 25)", "[25, )"))]

    p1 <- .plot(.dt, .dt[, .(n = .N), by = .(apoe.e4)], "APOE risk genotype", g, 
                aes(x = as.factor(apoe.e4), y = n),
                aes(x = as.factor(apoe.e4), y = ln.mu, fill = as.factor(pathoAD)),
                aes(x = as.factor(apoe.e4), y = ln.tot, fill = as.factor(pathoAD)))

    p2 <- .plot(.dt, .dt[, .(n = .N), by = .(msex)], "Sex", g, 
                aes(x = as.factor(msex), y = n),
                aes(x = as.factor(msex), y = ln.mu, fill = as.factor(pathoAD)),
                aes(x = as.factor(msex), y = ln.tot, fill = as.factor(pathoAD)))
    
    p3 <- .plot(.dt, .dt[, .(n = .N), by = .(braaksc)], "Braak Stage", g,
                aes(x = as.factor(braaksc), y = n),
                aes(x = as.factor(braaksc), y = ln.mu, fill = as.factor(pathoAD)),
                aes(x = as.factor(braaksc), y = ln.tot, fill = as.factor(pathoAD)))

    p4 <- .plot(.dt, .dt[, .(n = .N), by = .(age.death)], "Age at death", g, 
                aes(x = as.factor(age.death), y = n),
                aes(x = as.factor(age.death), y = ln.mu, fill = as.factor(pathoAD)),
                aes(x = as.factor(age.death), y = ln.tot, fill = as.factor(pathoAD)))
    
    p5 <- .plot(.dt, .dt[, .(n = .N), by = .(educ)], "Education (yrs)", g, 
                aes(x = as.factor(educ), y = n),
                aes(x = as.factor(educ), y = ln.mu, fill = as.factor(pathoAD)),
                aes(x = as.factor(educ), y = ln.tot, fill = as.factor(pathoAD)))

    return(p1 | p2 | p3 | p4 | p5)
}

.glob.dt <- dcast(stat.dt[tot > 0], hgnc + ensg ~ pathoAD + projid, value.var = "ln.mu")

gg <- .glob.dt %>% select(ensg, hgnc) %>% as_tibble
xx <- .glob.dt %>% select(starts_with("0_")) %>% as.matrix
yy <- .glob.dt %>% select(starts_with("1_")) %>% as.matrix

patho.test.dt <- run.t.test(gg, xx, yy)



## ggplot(patho.test.dt


g <- "PTPRG"
