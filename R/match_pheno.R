#!/usr/bin/env Rscript

col.file <- "result/temp/Per_1.cols.gz"
pheno.file <- "result/phenotyped.txt.gz"
pheno.col <- 1

argv <- commandArgs(trailingOnly = TRUE)

col.file <- argv[1]
pheno.file <- argv[2]
pheno.col <- as.integer(argv[3])
out.file <- argv[4]

library(tidyverse)
library(data.table)

col.dt <- fread(col.file, header=FALSE, col.names = "TAG")
pheno.dt <- fread(pheno.file, header=TRUE) %>% distinct
out.dt <- left_join(col.dt, pheno.dt, by = "TAG")
.pheno <- colnames(pheno.dt)[2 + pheno.col]

out.dt <- out.dt %>%
    select(all_of(.pheno))

fwrite(out.dt, file = out.file, col.names = FALSE, row.names=FALSE)
