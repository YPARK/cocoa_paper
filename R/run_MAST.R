## mtx.file <- "data/3_3_5_3_1.mtx.gz"
## label.file <- "data/3_3_5_3_1.label.gz"
## ind.file <- "data/3_3_5_3_1.ind.gz"

run.MAST <- function(mtx.file, 
                     label.file,
                     ind.file) {

    library(data.table)
    library(tidyverse)
    library(MAST)

    Y <- Matrix::readMM(mtx.file) %>% as.matrix
    ## Y.int <- apply(Y, 2, function(x) { storage.mode(x) <- 'integer'; x})

    colnames(Y.int) <- 1:ncol(Y.int)

    lab <- fread(label.file, header=FALSE) %>% unlist
    ind <- fread(ind.file, header=FALSE) %>% unlist %>%
        str_remove_all("Ind") %>% as.integer

    grp <- lab[ind]

    mm <- nrow(Y.int)
    nn <- ncol(Y.int)

    cdat <- data.frame(condition = grp)
    rdat <- data.frame(primerid = 1:mm)

    .dat <- FromMatrix(Y.int, cdat, rdat, check_sanity=FALSE)

    pc <- (rsvd::rpca(t(assay(.dat)), retx=TRUE, k=3)$x) %>%
        as.data.table

    colnames(pc) <- paste0("pc", 1:ncol(pc))
    colData(.dat)$condition <- factor(colData(.dat)$condition)
    colData(.dat) <- cbind(colData(.dat), pc)

    .take.dt <- function(.zlm.out) {

        out <- summary(.zlm.out, doLRT='condition1')
        .dt <- out$datatable
        out.dt <-
            .dt[contrast=='condition1' & component=='H',.(primerid, `Pr(>Chisq)`)]
        dplyr::rename(pv = `Pr(>Chisq)`) %>%
            dplyr::rename(x.col = primerid) %>%
            as.data.table
    }

    .dt3 <- zlm(~condition + pc1 + pc2 + pc3, .dat) %>%
        .take.dt %>% 
        mutate(method = "MAST-3")

    .dt2 <- zlm(~condition + pc1 + pc2, .dat) %>%
        .take.dt %>% 
        mutate(method = "MAST-2")

    .dt1 <- zlm(~condition + pc1, .dat) %>%
        .take.dt %>% 
        mutate(method = "MAST-1")

    .dt0 <- zlm(~condition, .dat) %>%
        .take.dt %>% 
        mutate(method = "MAST-0")

    rbind(.dt0, .dt1, .dt2, .dt3)
}
