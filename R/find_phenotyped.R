#!/usr/bin/env Rscript

## pheno.file <- "data/brain_2018-05-03/phenotypes.csv"
## meta.file <- "data/brain_2018-05-03/filtered_column_metadata.txt.gz"
argv <- commandArgs(trailingOnly = TRUE)

pheno.file <- argv[1]
meta.file <- argv[2]
out.file <- argv[3]

library(tidyverse)
library(data.table)

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

.pheno <- fread(pheno.file, header=TRUE)

.pheno[, age.death := sapply(age_death, .age.code)]
.pheno[, apoe.e4 := sapply(apoe_genotype, .apoe.code)]

.cols <- fread(meta.file, header=TRUE)
.cols <- left_join(.cols, .pheno, by = "projid")

.out <- .cols[, .(TAG, projid, pathoAD, msex, age.death, apoe.e4)]

fwrite(.out, out.file, col.names=TRUE, row.names=FALSE, sep="\t")
