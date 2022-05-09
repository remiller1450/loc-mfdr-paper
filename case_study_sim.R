####################################################################################################
#
#                                                                                     
#   Filename    :	case_study_sim.R										  
#   Input data files  :    ---                                                        
#   Output data files :    table5.pdf
#
#   Required R packages :  ggplot2, ncvreg, Matrix, survival, covTest, selectiveInference
#                          Rccp, gridExtra, locfdr, grid, reshape2, hdi, knockoff
#
#
####################################################################################################


library(ggplot2)
library(ncvreg)
library(Matrix)
library(survival)
library(covTest)
library(selectiveInference)
library(gridExtra)
library(locfdr)
library(grid)
library(reshape2)
library(hdi)
library(knockoff)

## Sets the working directory to current location of this folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Sources functions used in the simulations and the construction of figures
source("functions.R")

## Load TCGA BRCA1 dataset
TCGA <- readRDS("case_study_data//bcTCGA.rds")

## Standardize X matrix
XX <- std(TCGA$X)

## Find 5 max beta's (used to make simulated beta realistic)
init.fit = cv.ncvreg(X = XX, y = TCGA$y, penalty = "lasso")
sort(abs(init.fit$fit$beta[-1,which(init.fit$lambda == init.fit$lambda.min)]), decreasing = TRUE)[1:5]

## Seed
set.seed(133)

nrep = 100
tp = tot = matrix(NA, nrow = nrep, ncol = 4)

for(i in 1:nrep){

## 20 causal features
cause_id = sample(1:ncol(XX), size = 20, replace = FALSE)
beta = rep(0, ncol(XX))

## chosen to reflect slightly larger effects than present in real data
beta[cause_id] = c(rep(0.4, 10), rep(0.2, 10))

## Outcome
y = XX %*% beta + rnorm(nrow(XX), 0, sd = 1)

## Model fit
cv.fit <- cv.ncvreg(XX, y, penalty = "lasso")

## mfdr selections
cv.mfdr.res = local_mfdr(cv.fit$fit, cv.fit$lambda.min)
sel_id = which(cv.mfdr.res$mfdr < 0.1)

### CV selections
cv_sel_id = which(cv.fit$fit$beta[-1,which(cv.fit$lambda == cv.fit$lambda.min)] != 0)


### Univariate testing
tstat <- numeric(ncol(XX))
for (j in 1:ncol(XX)){
  fit.lm <- lm(y ~ XX[,j])
  tstat[j] <- summary(fit.lm)$coefficients[2,3]
}
zstat <- qnorm(pt(tstat,nrow(XX) - 2))
zstat[is.infinite(zstat)] <- tstat[is.infinite(zstat)] ## Make infinite z-stats into their pre-transformed t-stat

univariate.res <- locfdr(zstat, nulltype = 0, plot = 0, df = 10)
uni_sel_id = which(univariate.res$fdr < 0.1)


### Knock-off Filter
#suppressWarnings(knres <- knockoff.filter(X = XX, y = y, fdr = .1))

### selective inf
fit.lar <- tryCatch(lar(XX,y, maxsteps = 26), error=function(e) NULL)
sig.est <- tryCatch(estimateSigma(XX, y), error=function(e) NULL)
si.res <- tryCatch(larInf(fit.lar, sigma = sig.est$sigmahat, k = 25), error=function(e) NULL)

if (!is.null(si.res)){
  step.pv <- forwardStop(si.res$pv)
  
  if (step.pv == 0){
    si.tp = 0
  } else if (step.pv >= 1){
    si.tp = sum(si.res$vars[1:step.pv] %in% cause_id)
  }
} else {
    si.tp = 0
    step.pv = 0
}
  
  
tp[i,] = c(sum(cause_id %in% sel_id), sum(cause_id %in% cv_sel_id), sum(cause_id %in% uni_sel_id), si.tp)
tot[i,] = c(length(sel_id), length(cv_sel_id), sum(univariate.res$fdr < 0.1), step.pv)
}

tp[is.na(tp)] = 0
tot[is.na(tot)] = 0

res = list(tp, tot)
save(res, file = "case_study_sim.RData")

apply(tp, 2, mean)
apply(tot, 2, mean)



### Timings
cv_start = Sys.time()
cv.fit <- cv.ncvreg(XX, y, penalty = "lasso")
cv_stop = Sys.time()
cv_start - cv_stop

mfdr_start = Sys.time()
loc.mfdr <- local_mfdr(cv.fit$fit, cv.fit$lambda.min)
mfdr_stop = Sys.time()
mfdr_start - mfdr_stop

lar_start = Sys.time()
fit.lar <- lar(XX,y, maxsteps = 30)
lar_stop = Sys.time()
lar_start - lar_stop

si_start = Sys.time()
sig.est <- estimateSigma(XX, y)
si.res <- larInf(fit.lar, sigma = sig.est$sigmahat)
si_stop = Sys.time()
si_start - si_stop
