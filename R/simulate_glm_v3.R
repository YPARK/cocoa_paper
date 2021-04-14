#!/usr/bin/env Rscript

##############
## examples ##
##############

nind <- 40
ngene <- 1000
ncausal <- 5
ncovar.conf <- 1
ncovar.batch <- 3
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

if(length(argv) != 13)
{
    q()
}

nind <- as.integer(argv[1])
ngene <-  as.integer(argv[2])
ncausal <-  as.integer(argv[3])
ncovar.conf <- as.integer(argv[4])
ncovar.batch <- as.integer(argv[5])
ncell.ind <- as.integer(argv[6])

pve.1 <- as.numeric(argv[7])
pve.c <- as.numeric(argv[8])
pve.a <- as.numeric(argv[9])

rseed <- as.numeric(argv[10])

rho.a <- as.numeric(argv[11])
rho.b <- as.numeric(argv[12])

out <- argv[13]

dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)

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


#' random from collection
#' @param n1
#' @param n2
.rand <- function(n1, n2, .sample = c(-1, 1)) {
    matrix(sample(.sample, n1 * n2, TRUE), nrow = n1, ncol = n2)
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
#' @param ncovar.conf number of covar shared 
#' @param ncovar.batch number of covar on mu, batch effect
#' @param ngenes number of genes/features
#' @param ncausal number of causal genes
sample.seed.data <- function(nind, ncovar.conf, ncovar.batch, pve) {

    if(ncovar.conf > 0) {
        xx <- .rnorm(nind, ncovar.conf) # covariates
    } else {
        xx <- .zero(nind, 1) # empty
    }

    if(ncovar.batch > 0) {
        xx.mu <- .rnorm(nind, ncovar.batch) # covariates
    } else {
        xx.mu <- .zero(nind, 1) # empty
    }

    ## Biased assignment mechanism
    .delta <- .rnorm(ncol(xx), 1) / sqrt(ncol(xx))

    true.logit <- .rnorm(nind, 1)

    logit <- .scale(xx %*% .delta) * sqrt(pve)
    logit <- logit + true.logit * sqrt(1 - pve)

    ww <- rbinom(prob=.sigmoid(logit), n=nind, size=1)

    list(w = ww, lib = true.logit, x = xx, x.mu = xx.mu)
}

.param <- sample.seed.data(nind, ncovar.conf, ncovar.batch, pve.a)

nn <- length(.param$w)

xx <- .param$x       # covariates shared
xx.mu <- .param$x.mu # covariates on mu only 
ww <- .param$w       # stochastic assignment

causal <- sample(ngene, ncausal) # causal genes

#######################
## Treatment effects ##
#######################

ln.mu.w <- sapply(1:ngene, function(j) sample(ww) * rnorm(1))
ln.mu.w[, causal] <- sapply(1:ncausal, function(j) ww * rnorm(1))

#########################
## confounding effects ##
#########################

ln.mu.x <- xx %*% .rand(ncol(xx), ngene)

#########################
## other batch effects ##
#########################

.batch <- xx.mu %*% .rand(ncol(xx.mu), ngene)

ln.mu.x <- ln.mu.x + .batch

########################
## unstructured noise ##
########################

ln.mu.eps <- .rnorm(nn, ngene)

#############################
## combine all the effects ##
#############################

ln.mu <- .scale(ln.mu.w) * sqrt(pve.1) +
    .scale(ln.mu.x) * sqrt(pve.c) +
    .scale(ln.mu.eps) * sqrt(1 - pve.1 - pve.c)

ln.mu <- .scale(ln.mu)

mu <- exp(ln.mu)

##########################
## unconfounded signals ##
##########################

clean.ln.mu <- .scale(ln.mu.w) * sqrt(pve.1) +
    .scale(ln.mu.eps) * sqrt(1 - pve.1)

clean.mu <- exp(clean.ln.mu)

######################
## sequencing depth ##
######################

rr <- rgamma(ncell.ind * nn, shape=rho.a, scale=1/rho.b)

cells <- 1:length(rr)

.write(t(mu), file = gzfile(out %&% ".mu.gz"))
.write(t(clean.mu), file = gzfile(out %&% ".mu-clean.gz"))
.write(t(xx), file = gzfile(out %&% ".covar.gz"))
.write(data.frame(ww), file = gzfile(out %&% ".label.gz"))
.write(data.frame(cells), file = gzfile(out %&% ".cols.gz"))
.write(data.frame(rr), file = gzfile(out %&% ".rho.gz"))
.write(data.frame(sort(causal)), file = gzfile(out %&% ".causal.gz"))
