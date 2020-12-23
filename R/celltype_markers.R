argv <- commandArgs(trailingOnly = TRUE)

marker.file <- argv[1] # e.g., "data/PsychENCODE.marker"
row.file <- argv[2] # e.g., "result/merged.rows.gz"
out.file <- argv[3] # e.g., "output.txt.gz"

`%&%` <- function(a,b) paste0(a,b)

library(tidyverse)
library(data.table)

genes.dt <- fread(row.file, header=FALSE, col.names = "gene")
genes.dt[, c("ensg","hgnc") := tstrsplit(gene, "_", fixed=TRUE)]

ct.dt <- fread(marker.file, col.names =  c("hgnc", "celltype"))

annot.dt <- left_join(ct.dt, genes.dt, by="hgnc") %>%
    na.omit

.dt <- annot.dt[, .(gene, celltype)]
fwrite(.dt, out.file, sep = " ", col.names = FALSE)

