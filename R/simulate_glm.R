#!/usr/bin/env Rscript

##############
## examples ##
##############

nind <- 40
ngene <- 1000
ncausal <- 5
ncovar <- 3
ncell.ind <- 10
pve.1 <- 0.3
pve.c <- 0.3
pve.a <- 0.0
rseed <- 13
rho.a <- 2
rho.b <- 2

#######################
## command arguments ##
#######################

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 12)
{
    q()
}

nind <- as.integer(argv[1])
ngene <-  as.integer(argv[2])
ncausal <-  as.integer(argv[3])
ncovar <- as.integer(argv[4])
ncell.ind <- as.integer(argv[5])

pve.1 <- as.numeric(argv[6])
pve.c <- as.numeric(argv[7])
pve.a <- as.numeric(argv[8])

rseed <- as.numeric(argv[9])

rho.a <- as.numeric(argv[10])
rho.b <- as.numeric(argv[11])

out <- argv[12]

set.seed(rseed)

#' write data file
.write <- function(...) {
    write.table(..., row.names=F, col.names=F, sep=" ")
}

#' simple concat
`%&%` <- function(a,b) {
    paste0(a,b)
}

#' random standard Normal
#' @param n1
#' @param n2
.rnorm <- function(n1, n2) {
    matrix(rnorm(n1 * n2), nrow = n1, ncol = n2)
}

#' just a zero matrix
#' @param n1
#' @param n2
.zero <- function(n1, n2) {
    matrix(0, nrow = n1, ncol = n2)
}

.sigmoid <- function(x) {
    1/(1 + exp(-x))
}

.scale <- function(x) {
    .sd <- pmax(apply(x, 2, sd), 1e-8)
    sweep(x, 2, .sd, `/`)
}

#' sample model parameters
#' @param nind number of individuals/samples
#' @param ncovar number of other covariates
#' @param ngenes number of genes/features
#' @param ncausal number of causal genes
sample.seed.data <- function(nind, ncovar, pve) {

    if(ncovar > 0) {
        xx <- .rnorm(nind, ncovar) # covariates
    } else {
        xx <- .zero(nind, 1) # empty
    }

    ## assignment mechanism
    .delta <- .rnorm(ncol(xx), 1) / sqrt(ncol(xx))

    true.logit <- .rnorm(nind, 1)

    logit <- .scale(xx %*% .delta) * sqrt(pve)
    logit <- logit + true.logit * sqrt(1 - pve)

    ww <- rbinom(prob=.sigmoid(logit), n=nind, size=1)

    list(w = ww, lib = true.logit, x = xx)
}

.param <- sample.seed.data(nind, ncovar, pve.a)


nn <- length(.param$w)

xx <- .param$x   # covariates
ww <- .param$w   # stochastic assignment

causal <- sample(ngene, ncausal) # causal genes

ln.mu.w <- sapply(1:ngene, function(j) sample(ww) * rnorm(1))
ln.mu.w[, causal] <- sapply(1:ncausal, function(j) ww * rnorm(1))
ln.mu.x <- xx %*% .rnorm(ncol(xx), ngene)
ln.mu.eps <- .rnorm(nn, ngene)

ln.mu <- .scale(ln.mu.w) * sqrt(pve.1) +
    .scale(ln.mu.x) * sqrt(pve.c) +
    .scale(ln.mu.eps) * sqrt(1 - pve.1 - pve.c)

ln.mu <- .scale(ln.mu)

resid.ln.mu <- residuals(lm(ln.mu ~ ln.mu.x))

mu <- exp(ln.mu)
resid.mu <- exp(resid.ln.mu)

rr <- rgamma(ncell.ind * nn, shape=rho.a, scale=1/rho.b)

cells <- 1:length(rr)

.write(t(mu), file = gzfile(out %&% ".mu.gz"))
.write(t(resid.mu), file = gzfile(out %&% ".mu-clean.gz"))
.write(t(xx), file = gzfile(out %&% ".covar.gz"))
.write(data.frame(ww), file = gzfile(out %&% ".label.gz"))
.write(data.frame(cells), file = gzfile(out %&% ".cols.gz"))
.write(data.frame(rr), file = gzfile(out %&% ".rho.gz"))
.write(data.frame(sort(causal)), file = gzfile(out %&% ".causal.gz"))
