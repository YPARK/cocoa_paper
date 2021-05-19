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
pve.f <- 0.0

rseed <- 13
rho.a <- 2
rho.b <- 2

#######################
## command arguments ##
#######################

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 14)
{
    q()
}

nind <- as.integer(argv[1])
ngene <-  as.integer(argv[2])
ncausal <-  as.integer(argv[3])
ncovar.conf <- as.integer(argv[4])
ncovar.batch <- as.integer(argv[5])
ncell.ind <- as.integer(argv[6])

pve.1 <- as.numeric(argv[7]) # disease variability
pve.c <- as.numeric(argv[8]) # confounder variability on expression
pve.a <- as.numeric(argv[9]) # confounder variability on label
pve.f <- as.numeric(argv[10]) # feedback variability

rseed <- as.numeric(argv[11])

rho.a <- as.numeric(argv[12])
rho.b <- as.numeric(argv[13])

out <- argv[14]

dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)

set.seed(rseed)

stopifnot((pve.1 + pve.c + pve.f) < 1)

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
#' @param ncovar.batch number of covar on lambda, batch effect
#' @param ngenes number of genes/features
#' @param ncausal number of causal genes
sample.seed.data <- function(nind, ncovar.conf, ncovar.batch, pve) {

    if(ncovar.conf > 0) {
        xx <- .rnorm(nind, ncovar.conf) # covariates
    } else {
        xx <- .zero(nind, 1) # empty
    }

    if(ncovar.batch > 0) {
        xx.lambda <- .rnorm(nind, ncovar.batch) # covariates
    } else {
        xx.lambda <- .zero(nind, 1) # empty
    }

    ## Biased assignment mechanism
    .delta <- .rnorm(ncol(xx), 1) / sqrt(ncol(xx))

    true.logit <- .rnorm(nind, 1)

    logit <- .scale(xx %*% .delta) * sqrt(pve)
    logit <- logit + true.logit * sqrt(1 - pve)

    ww <- rbinom(prob=.sigmoid(logit), n=nind, size=1)

    list(w = ww, lib = true.logit, x = xx, x.lambda = xx.lambda)
}

.param <- sample.seed.data(nind, ncovar.conf, ncovar.batch, pve.a)

nn <- length(.param$w)

xx <- .param$x               # covariates shared --> confounding
xx.lambda <- .param$x.lambda # covariates on lambda only --> non-confounding
ww <- .param$w               # stochastic assignment

causal <- sample(ngene, ncausal) # causal genes

#######################
## Treatment effects ##
#######################

sample.w.rand <- function(j) {
    ## make sure that we don't create unlabeled positive genes
    r <- rnorm(length(ww))
    r <- .scale(matrix(residuals(lm(r ~ ww)), ncol=1))
    return(matrix(r * rnorm(1), ncol=1))
}

ln.lambda.w <- sapply(1:ngene, sample.w.rand)
ln.lambda.w[, causal] <- sapply(1:ncausal, function(j) ww * rnorm(1))

#########################
## confounding effects ##
#########################

ln.lambda.x <- xx %*% .rand(ncol(xx), ngene)

#########################
## other batch effects ##
#########################

.batch <- xx.lambda %*% .rand(ncol(xx.lambda), ngene)

ln.lambda.x <- ln.lambda.x + .batch

#########################
## Add feedback effect ##
#########################

.u <- rsvd::rsvd(ln.lambda.w[, causal, drop = FALSE], k=1)$u
ln.feedback <- .u %*% .rnorm(1, ngene)

########################
## unstructured noise ##
########################

ln.lambda.eps <- .rnorm(nn, ngene)

#############################
## combine all the effects ##
#############################

ln.lambda <-
    .scale(ln.lambda.w) * sqrt(pve.1) +
    .scale(ln.lambda.x) * sqrt(pve.c) +
    .scale(ln.feedback) * sqrt(pve.f) +
    .scale(ln.lambda.eps) * sqrt(1 - pve.1 - pve.c - pve.f)

lambda <- exp(ln.lambda)

##########################
## unconfounded signals ##
##########################

clean.ln.lambda <- .scale(ln.lambda.w) * sqrt(pve.1) +
    .scale(ln.lambda.eps) * sqrt(1 - pve.1)

clean.lambda <- exp(clean.ln.lambda)

######################
## sequencing depth ##
######################

rr <- rgamma(ncell.ind * nn, shape=rho.a, scale=1/rho.b)

cells <- 1:length(rr)

.write(t(lambda), file = gzfile(out %&% ".lambda.gz"))
.write(t(clean.lambda), file = gzfile(out %&% ".lambda-clean.gz"))
.write(t(xx), file = gzfile(out %&% ".covar.gz"))
.write(data.frame(ww), file = gzfile(out %&% ".label.gz"))
.write(data.frame(cells), file = gzfile(out %&% ".cols.gz"))
.write(data.frame(rr), file = gzfile(out %&% ".rho.gz"))
.write(data.frame(sort(causal)), file = gzfile(out %&% ".causal.gz"))
