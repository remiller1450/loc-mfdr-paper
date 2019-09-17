####################################################################################################
#
#   Filename    :	loc_mfdr_reproduce.R				
#
#   Output data files :    Fig1.png, Fig2.png, Fig3.png, Fig4.png, Fig5.png
#			                     table1.pdf, table2.pdf, table3.pdf, table4.pdf 
#
#   Required packages (CRAN) :  ggplot2, Matrix, survival, covTest, selectiveInference, devtools,
#                               Rccp, gridExtra, locfdr, grid, reshape2, hdi, knockoff, ashr,
#                               dplyr, pROC
#
#   Required packages (github) : ncvreg
#
####################################################################################################

library(ggplot2)
library(Matrix)
library(survival)
library(covTest)
library(selectiveInference)
library(Rccp)
library(gridExtra)
library(locfdr)
library(grid)
library(reshape2)
library(hdi)
library(knockoff)
library(ashr)
library(devtools)
library(dplyr)
library(pROC)

### Install latest version of ncvreg
#devtools::install_github("pbreheny/ncvreg")
library(ncvreg)

## Sets the working directory to current location of this folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Sources functions used in the simulations and the construction of figures
source("functions.R")

########################################################################################
#
#                           Code to produce Figure 1
#
########################################################################################

### Set Seed
set.seed(13351) 

### Setup parameters
n <- 200
p <- 100
bb <- c(3.5,0,-3.5,0)

### Generate data using funs sourced
D1 <- genData(n, J=2, J0=2, K=2, K0=1, rho=0, rho.g=0.83, beta=bb)  #### 9 correlated vars for each true (type B)
D2 <- genData(n, p - 4, rho=0, beta=0, corr="auto")
X <- cbind(D1$X, .6*D2$X)
y <- D1$y

## Fit cross-validated model
fit <- cv.ncvreg(X,y, penalty = "lasso", returnX = TRUE)
lseq <- fit$lambda

## Find loc mfdrs for all lambda values
locfdr.res <- matrix(NA, nrow = ncol(X), ncol = length(lseq))
sel.mat <- matrix(NA, nrow = ncol(X), ncol = length(lseq))
for(i in 1:length(lseq)){
  temp <- locmfdr.plot1(fit$fit, lambda = lseq[i], number = ncol(X))
  locfdr.res[,i] <- temp$locfdr
  sel.mat[,i] <- (temp$selected == '*')
}

## Setup color pal for plotting
cols <- numeric(p)
cids <- pal(3)
cols[c(1,3)] <- cids[1]
cols[c(2,4)] <- cids[2]
cols[5:p] <- cids[3]

### Format results for plotting
colnames(locfdr.res) <- lseq
wide <- data.frame(var = rownames(temp), col = cols, locfdr.res)
long <- melt(wide, id.vars = c("var","col"))
long$lambda<- -as.numeric(substr(long$variable,2,10))
wide.sel <- data.frame(var = rownames(temp), sel.mat)
long.sel <- melt(wide.sel, id.vars = c('var'))
long <- data.frame(long, sel = long.sel[,3])
nexts <- which(diff(long[order(long$var),]$sel) == 1) + 1
nexts.pos <- numeric(nrow(long))
nexts.pos[nexts] <- 1
long <- long[order(long$var),]

#### Get Univariate Testing results
z <- numeric(p)
for(i in 1:p){
  t <- summary(lm(y ~ X[,i]))$coefficients[2,3]
  z[i] <- qnorm(pt(t,p-2))
}
d <- density(z)
dd <- approxfun(d$x, d$y)
fdr.u <- dnorm(z)/dd(z)
fdr.u[fdr.u > 1] <- 1
fdr.points <- data.frame(value = fdr.u, lambda = -(max(lseq)+.25), col = cols)

## Upper right panel plot
p1 <- ggplot() +
  geom_line(data = rbind(long[long$sel,],long[nexts,]),aes(x=-lambda, y=value, group=var, col = factor(col)), size=1, linetype = "solid") + 
  geom_line(data = rbind(long[!long$sel,], long[nexts,]), aes(x=-lambda, y=value, group=var, col = factor(col)), size=1, linetype = "dashed") +
  geom_vline(xintercept=fit$lambda.min, linetype = "dotted") + 
  scale_x_reverse(limits=c((max(lseq)+.25),0)) +
  geom_point(data = fdr.points, mapping = aes(x = -lambda, y = value, col = factor(col)), shape = 2, size = 1.5) +
  scale_color_manual(name = "Feature Type", values = cids[c(1,3,2)], labels = c("Causal", "Correlated", "Noise")) +
  labs(x = expression(lambda), y = expression(widehat(mfdr))) + theme(legend.position="bottom") + ggtitle("")

## Lower left panel plot 
long.coef <- melt(fit$fit$beta[-1,])
coef.path <- data.frame(long.coef, col = cols)
p2 <- ggplot() + geom_line(data = coef.path, aes(x = Var2, y = value, group = Var1, col = factor(col)), size = 1) + scale_x_reverse() +
  scale_color_manual(name = "Feature Type", values = cids[c(1,3,2)], labels = c("Causal", "Correlated", "Noise")) +
  geom_vline(xintercept=fit$lambda.min, linetype = "dotted") + 
  labs(x = expression(lambda), y = "Coefficient Estimate") + ggtitle("")


### Reset for upper panel plots using larger signal
set.seed(31236) ## Get the same X matrix using genData, but use a different beta vector
bb <- c(4,0,-4,0)

## Generate data
D1 <- genData(n, J=2, J0=2, K=2, K0=1, rho=0, rho.g=0.8, beta=bb)  #### 9 correlated vars for each true (type B)
D2 <- genData(n, p - 4, rho=0, beta=0, corr="auto")
X <- cbind(D1$X, .585*D2$X)
y <- D1$y
fit <- cv.ncvreg(X,y, penalty = "lasso", returnX = TRUE)
lseq <- fit$lambda

## Local mfdr for each lambda
locfdr.res <- matrix(NA, nrow = ncol(X), ncol = length(lseq))
sel.mat <- matrix(NA, nrow = ncol(X), ncol = length(lseq))
for(i in 1:length(lseq)){
  temp <- locmfdr.plot1(fit$fit, lambda = lseq[i], number = ncol(X))
  locfdr.res[,i] <- temp$locfdr
  sel.mat[,i] <- (temp$selected == '*')
}


### Format results for plotting
colnames(locfdr.res) <- lseq
wide <- data.frame(var = rownames(temp), col = cols, locfdr.res)
long <- melt(wide, id.vars = c("var","col"))
long$lambda<- -as.numeric(substr(long$variable,2,10))

wide.sel <- data.frame(var = rownames(temp), sel.mat)
long.sel <- melt(wide.sel, id.vars = c('var'))
long <- data.frame(long, sel = long.sel[,3])

nexts <- which(diff(long[order(long$var),]$sel) == 1) + 1
nexts.pos <- numeric(nrow(long))
nexts.pos[nexts] <- 1
long <- long[order(long$var),]

#### Univariate results
z <- numeric(p)
for(i in 1:p){
  t <- summary(lm(y ~ X[,i]))$coefficients[2,3]
  z[i] <- qnorm(pt(t,p-2))
}
d <- density(z)
dd <- approxfun(d$x, d$y)
fdr.u <- dnorm(z)/dd(z)
fdr.u[fdr.u > 1] <- 1
fdr.points <- data.frame(value = fdr.u, lambda = -(max(lseq)+.25), col = cols)

## Upper Right plot
p3 <- ggplot() +
  geom_line(data = rbind(long[long$sel,],long[nexts,]),aes(x=-lambda, y=value, group=var, col = factor(col)), size=1, linetype = "solid") + 
  geom_line(data = rbind(long[!long$sel,], long[nexts,]), aes(x=-lambda, y=value, group=var, col = factor(col)), size=1, linetype = "dashed") +
  geom_vline(xintercept=fit$lambda.min, linetype = "dotted") + 
  scale_x_reverse(limits=c((max(lseq)+.25),0)) +
  geom_point(data = fdr.points, mapping = aes(x = -lambda, y = value, col = factor(col)), shape = 2, size = 1.5) +
  scale_color_manual(name = "Feature Type", values = cids[c(1,3,2)], labels = c("Causal", "Correlated", "Noise")) +
  labs(x = expression(lambda), y = expression(widehat(mfdr))) + theme(legend.position="bottom") + ggtitle("mfdr Path")

## Upper Left plot
long.coef <- melt(fit$fit$beta[-1,])
coef.path <- data.frame(long.coef, col = cols)
p4 <- ggplot() + geom_line(data = coef.path, aes(x = Var2, y = value, group = Var1, col = factor(col)), size = 1) + scale_x_reverse() +
  scale_color_manual(name = "Feature Type", values = cids[c(1,3,2)], labels = c("Causal", "Correlated", "Noise")) +
  geom_vline(xintercept=fit$lambda.min, linetype = "dotted") + 
  labs(x = expression(lambda), y = "Coefficient Estimate") + ggtitle("Coefficient Path")

## Setup common legend for all plots
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(p1)


## Assemble the 4 plots into a single figure and store
png("figures_tables/Fig1.png", h=7, w=8.5, units = 'in', res = 300)
grid.arrange(arrangeGrob(p4 + theme(legend.position="none"),
                         p3 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"),
                         p1 + theme(legend.position="none"),
                         nrow=2),
             mylegend, nrow=2,heights=c(10, 1))
dev.off()

########################################################################################
#
#           Simulations needed for Figure 1, Table 1, Figure 2
#
########################################################################################

## There are 6 different simulations that need to run to produce Figure 2 and Table 1
## The code below will run those simulations, saving the result after each simulation


############################################
##### Cox reg "assumptions met" scenario
############################################ 

### Sim Seed
set.seed(12345)

### Setup sim parameters
n <- 1000
p <- 600
t <- 60
beta <- c(rep(.15,t/2), rep(-.15,t/2), rep(0,p-t))

### Variable IDs for sim type
id.A <- which(beta != 0)
id.B <- NA   ## No correlated vars in this scenario
id.C <- which(beta == 0)

## Number of simulation replications
nreps <- 200

## Setup objects to store results
res.lassomfdr <- res.lassocv <- res.lassocv1se <- res.uni <- matrix(NA, nrow = p, ncol = nreps)
res.lassocv.kernel <- res.lassocv1se.kernel <- res.uniash <- matrix(NA, nrow = p, ncol = nreps)

for (i in 1:nreps){
  
  ## Generate Survival Data
  D1 <- genData(n, p, beta=beta, family = "survival")
  X <- D1$X
  
  ## 10% random censoring
  c <- rbinom(n, 1, p = .9)
  ct <- cbind(runif(sum(c == 0), min = 0, max = D1$y[which(c == 0)]), D1$y[which(c == 0)])
  D1$y[which(c == 0)] = ct[,1]
  y <- cbind(D1$y, c)
  
  ## Univariate locfdr test stats
  zstat <- tstat <- rep(NA, ncol(X))
  for (j in 1:ncol(X)){
    fit.lm <- coxph(Surv(y) ~ X[,j])
    zstat[j] <- summary(fit.lm)$coefficients[1,4]
  }
  
  ## Univariate using locfdr
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)$fdr
  
  ## Univariate using ashr
  univariate.res.ash <- get_lfdr(ash(beta = zstat, se = rep(1, length(zstat))))
  
  ### Fit lasso model
  cv.fit <- cv.ncvsurv(X,y,penalty = "lasso", returnX = TRUE)
  fit <- cv.fit$fit
  
  ## positions of lambda values of interest
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv1se.lam <- min(which(cv.fit$cve - cv.fit$cvse <  min(cv.fit$cve)))
  cv.lam <- which(cv.fit$fit$lambda == cv.fit$lambda.min)
  
  ## local mfdr using model/lambdas from above
  locfdr.cv  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "ashr")$mfdr
  locfdr.lam <- ncvreg::local_mfdr(fit, fit$lambda[mfdr.lam], method = "ashr")$mfdr
  locfdr.cv1se <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "ashr")$mfdr
  
  locfdr.cv.k  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "kernel")$mfdr
  locfdr.cv1se.k <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "kernel")$mfdr
  
  
  ## Store iteration results (local fdrs)
  res.lassomfdr[,i] <- locfdr.lam
  res.lassocv[,i] <- locfdr.cv
  res.lassocv1se[,i] <- locfdr.cv1se
  res.uni[,i] <- univariate.res
  res.uniash[,i] <- univariate.res.ash
  res.lassocv.kernel[,i] <- locfdr.cv.k
  res.lassocv1se.kernel[,i] <- locfdr.cv1se.k
}

## Save results
save(res.lassomfdr,res.lassocv, res.lassocv1se, res.uni, res.uniash, res.lassocv.kernel, res.lassocv1se.kernel,
     id.A, id.B, id.C, file = "simres/cox_easy.RData")


############################################
##### Cox reg "assumptions violated" scenario
############################################ 

### Sim Seed
set.seed(12345)

### Setup
n <- 200
p <- 600
t <- 60
bb <- numeric(60)
bb[(0:5)*10+1] <- c(.6, -.6, .5, -.5, .4, -.4)  #### 6 true variables (type A)

## Ids of var types
id.A <- which(bb != 0)
id.B <- 1:t
id.B <- id.B[-id.A]
id.C <- (t+1):p


## Number of simulation replications
nreps <- 200

## Setup objects to store results
res.lassomfdr <- res.lassocv <- res.lassocv1se <- res.uni <- matrix(NA, nrow = p, ncol = nreps)
res.lassocv.kernel <- res.lassocv1se.kernel <- res.uniash <- matrix(NA, nrow = p, ncol = nreps)

for(i in 1:nreps){
  
  ## Generate correlated survival data
  D1 <- genData(n, J=6, J0=6, K=10, K0=1, rho=0, rho.g=0.5, beta=bb, family = "survival")  ## A and B variables
  D2 <- genData(n, p - t, rho=0.8, beta=0, corr="auto") ## Noise "C" variables 
  X <- cbind(D1$X, D2$X)
  
  ## Randomly censor 10% of outcomes
  c <- rbinom(n, 1, p = .9)
  ct <- cbind(runif(sum(c == 0), min = 0, max = D1$y[which(c == 0)]), D1$y[which(c == 0)])
  D1$y[which(c == 0)] = ct[,1]
  y <- cbind(D1$y, c)
  
  #### Univariate locfdr
  zstat <- tstat <- rep(NA, ncol(X))
  for (j in 1:ncol(X)){
    fit.lm <- coxph(Surv(y) ~ X[,j])
    zstat[j] <- summary(fit.lm)$coefficients[1,4]
  }

  ## Univariate using locfdr
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)$fdr
  
  ## Univariate using ashr
  univariate.res.ash <- get_lfdr(ash(beta = zstat, se = rep(1, length(zstat))))
  
  ### Fit lasso model
  cv.fit <- cv.ncvsurv(X,y,penalty = "lasso", returnX = TRUE)
  fit <- cv.fit$fit
  
  ## positions of lambda values of interest
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv1se.lam <- min(which(cv.fit$cve - cv.fit$cvse <  min(cv.fit$cve)))
  cv.lam <- which(cv.fit$fit$lambda == cv.fit$lambda.min)
  
  ## local mfdr using model/lambdas from above
  locfdr.cv  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam],  method = "ashr")$mfdr
  locfdr.lam <- ncvreg::local_mfdr(fit, fit$lambda[mfdr.lam],  method = "ashr")$mfdr
  locfdr.cv1se <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam],  method = "ashr")$mfdr
  
  locfdr.cv.k  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "kernel")$mfdr
  locfdr.cv1se.k <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "kernel")$mfdr
  
  ## Store iteration results (local fdrs)
  res.lassomfdr[,i] <- locfdr.lam
  res.lassocv[,i] <- locfdr.cv
  res.lassocv1se[,i] <- locfdr.cv1se
  res.uni[,i] <- univariate.res
  res.uniash[,i] <- univariate.res.ash
  res.lassocv.kernel[,i] <- locfdr.cv.k
  res.lassocv1se.kernel[,i] <- locfdr.cv1se.k
}

## Save simulation results
save(res.lassomfdr,res.lassocv, res.lassocv1se, res.uni, res.uniash, res.lassocv.kernel, res.lassocv1se.kernel,
     id.A, id.B, id.C, file = "simres/cox_hard.RData")



####################################################
##### Logistic reg "assumptions met" scenario
####################################################

### Sim Seed
set.seed(12345)

### Setup sim parameters
n <- 1000
p <- 600
t <- 60
beta <- c(rep(.25,t/2), rep(-.25,t/2), rep(0,p-t))

### Variable IDs for sim type
id.A <- which(beta != 0)
id.B <- NA   ## No correlated vars in this scenario
id.C <- which(beta == 0)

## Number of simulation replications
nreps <- 200

## Setup objects to store results
res.lassomfdr <- res.lassocv <- res.lassocv1se <- res.uni <- matrix(NA, nrow = p, ncol = nreps)
res.lassocv.kernel <- res.lassocv1se.kernel <- res.uniash <- matrix(NA, nrow = p, ncol = nreps)

for (i in 1:nreps){
  
  ## Simulate Data
  D1 <- genData(n, p, beta=beta, family = "binomial")  
  X <- D1$X
  y <- D1$y
  
  ### Estimate univariate locfdr
  zstat <- tstat <- rep(NA, ncol(X))
  for (j in 1:ncol(X)){
    fit.lm <- glm(y ~ X[,j], family = "binomial")
    zstat[j] <- summary(fit.lm)$coefficients[2,3]
  } 
  
  ## Univariate using locfdr
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)$fdr
  
  ## Univariate using ashr
  univariate.res.ash <- get_lfdr(ash(beta = zstat, se = rep(1, length(zstat))))
  
  ### Fit lasso model
  cv.fit <- cv.ncvreg(X,y,penalty = "lasso", returnX = TRUE, family = "binomial")
  fit <- cv.fit$fit
  
  ## positions of lambda values of interest
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv1se.lam <- min(which(cv.fit$cve - cv.fit$cvse <  min(cv.fit$cve)))
  cv.lam <- which(cv.fit$fit$lambda == cv.fit$lambda.min)
  
  ## local mfdr using model/lambdas from above
  locfdr.cv  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam],  method = "ashr")$mfdr
  locfdr.lam <- ncvreg::local_mfdr(fit, fit$lambda[mfdr.lam],  method = "ashr")$mfdr
  locfdr.cv1se <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam],  method = "ashr")$mfdr
  
  locfdr.cv.k  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "kernel")$mfdr
  locfdr.cv1se.k <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "kernel")$mfdr
  
  ## Store iteration results (local fdrs)
  res.lassomfdr[,i] <- locfdr.lam
  res.lassocv[,i] <- locfdr.cv
  res.lassocv1se[,i] <- locfdr.cv1se
  res.uni[,i] <- univariate.res
  res.uniash[,i] <- univariate.res.ash
  res.lassocv.kernel[,i] <- locfdr.cv.k
  res.lassocv1se.kernel[,i] <- locfdr.cv1se.k
}

## Save simulation results
save(res.lassomfdr,res.lassocv, res.lassocv1se, res.uni, res.uniash, res.lassocv.kernel, res.lassocv1se.kernel,
     id.A, id.B, id.C, file = "simres/logistic_easy.RData")


####################################################
##### Logistic reg "assumptions violated" scenario
####################################################

### Sim Seed
set.seed(12345)

### Setup
n <- 200
p <- 600
bb <- numeric(60)
bb[(0:5)*10+1] <- c(1.2,-1.2,1,-1,.8,-.8)  #### 6 true variables (type A)

## Ids of var types
id.A <- which(bb != 0)
id.B <- 1:t
id.B <- id.B[-id.A]
id.C <- (t+1):p


## Number of simulation replications
nreps <- 200

## Setup objects to store results
res.lassomfdr <- res.lassocv <- res.lassocv1se <- res.uni <- matrix(NA, nrow = p, ncol = nreps)
res.lassocv.kernel <- res.lassocv1se.kernel <- res.uniash <- matrix(NA, nrow = p, ncol = nreps)

for(i in 1:nreps){
  ## Simulate data
  D1 <- genData(n, J=6, J0=6, K=10, K0=1, rho=0, rho.g=0.5, beta=bb, family = "binomial")  
  D2 <- genData(n, p - 60, rho=0.8, beta=0, corr="auto")
  
  X <- cbind(D1$X, D2$X)
  y <- D1$y
  
  #### Univariate locfdr
  zstat <- tstat <- rep(NA, ncol(X))
  for (j in 1:ncol(X)){
    fit.lm <- glm(y ~ X[,j], family = "binomial")
    zstat[j] <- summary(fit.lm)$coefficients[2,3]
  }

  ## Univariate using locfdr
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)$fdr
  
  ## Univariate using ashr
  univariate.res.ash <- get_lfdr(ash(beta = zstat, se = rep(1, length(zstat))))
  
  ### Fit lasso model
  cv.fit <- cv.ncvreg(X,y,penalty = "lasso", returnX = TRUE, family = "binomial")
  fit <- cv.fit$fit
  
  ## positions of lambda values of interest
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv1se.lam <- min(which(cv.fit$cve - cv.fit$cvse <  min(cv.fit$cve)))
  cv.lam <- which(cv.fit$fit$lambda == cv.fit$lambda.min)
  
  ## local mfdr using model/lambdas from above
  locfdr.cv  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "ashr")$mfdr
  locfdr.lam <- ncvreg::local_mfdr(fit, fit$lambda[mfdr.lam], method = "ashr")$mfdr
  locfdr.cv1se <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "ashr")$mfdr
  
  locfdr.cv.k  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "kernel")$mfdr
  locfdr.cv1se.k <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "kernel")$mfdr
  
  ## Store iteration results (local fdrs)
  res.lassomfdr[,i] <- locfdr.lam
  res.lassocv[,i] <- locfdr.cv
  res.lassocv1se[,i] <- locfdr.cv1se
  res.uni[,i] <- univariate.res
  res.uniash[,i] <- univariate.res.ash
  res.lassocv.kernel[,i] <- locfdr.cv.k
  res.lassocv1se.kernel[,i] <- locfdr.cv1se.k
}

## Save simulation results
save(res.lassomfdr,res.lassocv, res.lassocv1se, res.uni, res.uniash, res.lassocv.kernel, res.lassocv1se.kernel,
     id.A, id.B, id.C, file = "simres/logistic_hard.RData")


####################################################
##### Linear reg "assumptions met" scenario
####################################################

### Sim Seed
set.seed(12345)

### Setup
n <- 1000
p <- 600
t <- 60
sig <- sqrt(n)
beta <- c(rep(3,t/2), rep(-3,t/2), rep(0,p-t))

## Ids of var types
id.A <- which(bb != 0)
id.B <- 1:t
id.B <- id.B[-id.A]
id.C <- (t+1):p

## Number of simulation replications
nreps <- 200

## Setup objects to store results
res.lassomfdr <- res.lassocv <- res.lassocv1se <- res.uni <- matrix(NA, nrow = p, ncol = nreps)
res.lassocv.kernel <- res.lassocv1se.kernel <- res.uniash <- matrix(NA, nrow = p, ncol = nreps)

for (i in 1:nreps){
  ### Simulate data
  X <- matrix(rnorm(n*p), ncol = p, nrow = n)
  y <- X %*% beta + rnorm(n,sd = sig)
  
  ### Estimate univariate locfdr
  tstat <- rep(NA, ncol(X))
  for (j in 1:ncol(X)){
    fit.lm <- lm(y ~ X[,j])
    tstat[j] <- summary(fit.lm)$coefficients[2,3]
  }
  zstat <- qnorm(pt(tstat,n - 2))
  zstat[is.infinite(zstat)] <- tstat[is.infinite(zstat)]  ## force infinite z-stats back to their pre-transformed t-stat
  
  ## Univariate using locfdr
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)$fdr
  
  ## Univariate using ashr
  univariate.res.ash <- get_lfdr(ash(beta = zstat, se = rep(1, length(zstat))))
  
  ### Fit lasso model
  cv.fit <- cv.ncvreg(X,y,penalty = "lasso", returnX = TRUE)
  fit <- cv.fit$fit
  
  ## positions of lambda values of interest
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv1se.lam <- min(which(cv.fit$cve - cv.fit$cvse <  min(cv.fit$cve)))
  cv.lam <- which(cv.fit$fit$lambda == cv.fit$lambda.min)
  
  ## local mfdr using model/lambdas from above
  locfdr.cv  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "ashr")$mfdr
  locfdr.lam <- ncvreg::local_mfdr(fit, fit$lambda[mfdr.lam],  method = "ashr")$mfdr
  locfdr.cv1se <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam],  method = "ashr")$mfdr
  
  locfdr.cv.k  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "kernel")$mfdr
  locfdr.cv1se.k <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "kernel")$mfdr
  
  ## Store iteration results (local fdrs)
  res.lassomfdr[,i] <- locfdr.lam
  res.lassocv[,i] <- locfdr.cv
  res.lassocv1se[,i] <- locfdr.cv1se
  res.uni[,i] <- univariate.res
  res.uniash[,i] <- univariate.res.ash
  res.lassocv.kernel[,i] <- locfdr.cv.k
  res.lassocv1se.kernel[,i] <- locfdr.cv1se.k
}

## Save simulation results
save(res.lassomfdr,res.lassocv, res.lassocv1se, res.uni, res.uniash, res.lassocv.kernel, res.lassocv1se.kernel,
     id.A, id.B, id.C, file = "simres/linear_easy.RData")

####################################################
##### Linear reg "assumptions violated" scenario
####################################################

### Sim Seed
set.seed(12345)

### Setup
n <- 200
p <- 600
bb <- numeric(60)
bb[(0:5)*10+1] <- c(6,-6,5,-5,4,-4)  #### 6 true variables (type A)

## Ids of var types
id.A <- which(bb != 0)
id.B <- 1:t
id.B <- id.B[-id.A]
id.C <- (t+1):p

## Number of simulation replications
nreps <- 200

## Setup objects to store results
res.lassomfdr <- res.lassocv <- res.lassocv1se <- res.uni <- matrix(NA, nrow = p, ncol = nreps)
res.lassocv.kernel <- res.lassocv1se.kernel <- res.uniash <- matrix(NA, nrow = p, ncol = nreps)

for (i in 1:nreps){
  ### Simulate correlated data
  D1 <- genData(n, J=6, J0=6, K=10, K0=1, rho=0, rho.g=0.5, beta=bb)  #### 9 correlated vars for each true (type B)
  D2 <- genData(n, p - 60, rho=0.8, beta=0, corr="auto")
  X <- cbind(D1$X, D2$X)
  y <- D1$y
  
  ### Estimate univariate locfdr
  tstat <- rep(NA, ncol(X))
  for (j in 1:ncol(X)){
    fit.lm <- lm(y ~ X[,j])
    tstat[j] <- summary(fit.lm)$coefficients[2,3]
  }
  zstat <- qnorm(pt(tstat,n - 2))
  zstat[is.infinite(zstat)] <- tstat[is.infinite(zstat)]  ## force infinite z-stats back to their pre-transformed t-stat
  
  ## Univariate using locfdr
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)$fdr
  
  ## Univariate using ashr
  univariate.res.ash <- get_lfdr(ash(beta = zstat, se = rep(1, length(zstat))))
  
  ### Fit lasso model
  cv.fit <- cv.ncvreg(X,y,penalty = "lasso", returnX = TRUE)
  fit <- cv.fit$fit
  
  ## positions of lambda values of interest
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv1se.lam <- min(which(cv.fit$cve - cv.fit$cvse <  min(cv.fit$cve)))
  cv.lam <- which(cv.fit$fit$lambda == cv.fit$lambda.min)
  
  ## local mfdr using model/lambdas from above
  locfdr.cv  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "ashr")$mfdr
  locfdr.lam <- ncvreg::local_mfdr(fit, fit$lambda[mfdr.lam], method = "ashr")$mfdr
  locfdr.cv1se <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "ashr")$mfdr
  
  locfdr.cv.k  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "kernel")$mfdr
  locfdr.cv1se.k <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "kernel")$mfdr
  
  ## Store iteration results (local fdrs)
  res.lassomfdr[,i] <- locfdr.lam
  res.lassocv[,i] <- locfdr.cv
  res.lassocv1se[,i] <- locfdr.cv1se
  res.uni[,i] <- univariate.res
  res.uniash[,i] <- univariate.res.ash
  res.lassocv.kernel[,i] <- locfdr.cv.k
  res.lassocv1se.kernel[,i] <- locfdr.cv1se.k
}

## Save simulation results
save(res.lassomfdr,res.lassocv, res.lassocv1se, res.uni, res.uniash, res.lassocv.kernel, res.lassocv1se.kernel,
     id.A, id.B, id.C, file = "simres/linear_hard.RData")


##########################################################################################
## At this point six simulations have run and their results are stored as .RData files
## The code below uses these to generate Fig2 and Table1
##########################################################################################


########################################################################################
#
#           Code to produce Figure 2 (requires simulation results from earlier)
#
########################################################################################

# First general lower row of fig2 (assumptions violated)
load("simres/linear_hard.RData")

## Setup results data.frames
df.cv <- data.frame(Estimate = c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])),
                    Noise = c(rep(0, 60*ncol(res.lassocv)),
                              rep(1, 540*ncol(res.lassocv))))            
  
df.cv1se <- data.frame(Estimate = c(as.vector(res.lassocv1se[1:60,]), as.vector(res.lassocv1se[61:600,])),
                        Noise = c(rep(0, 60*ncol(res.lassocv1se)),
                                  rep(1, 540*ncol(res.lassocv1se))))     

df.mfdr <- data.frame(Estimate = c(as.vector(res.lassomfdr[1:60,]), as.vector(res.lassomfdr[61:600,])),
                      Noise = c(rep(0, 60*ncol(res.lassomfdr)),
                                rep(1, 540*ncol(res.lassomfdr))))    

df.uni <-  data.frame(Estimate = c(as.vector(res.uni[1:60,]), as.vector(res.uni[61:600,])),
                      Noise = c(rep(0, 60*ncol(res.uni)),
                                rep(1, 540*ncol(res.uni))))     

df.uniash <-  data.frame(Estimate = c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])),
                         Noise = c(rep(0, 60*ncol(res.uniash)),
                                   rep(1, 540*ncol(res.uniash))))     

### Loess curves
xx <- seq(0, 1, by = .01)

l.cv <- with(df.cv, loess(Noise ~ Estimate, span = .15))
yy.cv <- predict(l.cv, newdata = xx)
smooth.cv <- data.frame(xx = xx, yy = yy.cv)

l.cv1se <- with(df.cv1se, loess(Noise ~ Estimate, span = .15))
yy.cv1se <- predict(l.cv1se, newdata = xx)
smooth.cv1se <- data.frame(xx = xx, yy = yy.cv1se)

l.mfdr <- with(df.mfdr, loess(Noise ~ Estimate, span = .15))
yy.mfdr <- predict(l.mfdr, newdata = xx)
smooth.mfdr <- data.frame(xx = xx, yy = yy.mfdr)

l.uni <- with(df.uni, loess(Noise ~ Estimate, span = .55))
yy.uni <- predict(l.uni, newdata = xx)
smooth.uni <- data.frame(xx = xx, yy = yy.uni)

l.uniash <- with(df.uniash, loess(Noise ~ Estimate, span = .15))
yy.uniash <- predict(l.uniash, newdata = xx)
smooth.uniash <- data.frame(xx = xx, yy = yy.uniash)


## There are too many points, only plot a sample for the scatterplot
sample.id <- sample(1:nrow(df.cv), size = 8000, replace = FALSE)

## Create scatterplots from random subset of points
p.cv <- ggplot(df.cv[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.cv1se <- ggplot(df.cv1se[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.mfdr <- ggplot(df.mfdr[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.uni <- ggplot(df.uni[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.uniash <- ggplot(df.uniash[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)


## Add loess curves from full sim results
p.cv <- p.cv + geom_path(data = smooth.cv, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "Proportion (noise)")

p.cv1se <- p.cv1se + geom_path(data = smooth.cv1se, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "")

p.mfdr <- p.mfdr + geom_path(data = smooth.mfdr, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "")

p.uni <- p.uni + geom_path(data = smooth.uni, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "")

p.uniash <- p.uniash + geom_path(data = smooth.uniash, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "")


### Now generate upper row of Fig2 (assumptions met)
load("simres/linear_easy.RData")

## Setup results data.frames
df.cv <- data.frame(Estimate = c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])),
                    Noise = c(rep(0, 60*ncol(res.lassocv)),
                              rep(1, 540*ncol(res.lassocv))))            

df.cv1se <- data.frame(Estimate = c(as.vector(res.lassocv1se[1:60,]), as.vector(res.lassocv1se[61:600,])),
                       Noise = c(rep(0, 60*ncol(res.lassocv1se)),
                                 rep(1, 540*ncol(res.lassocv1se))))     

df.mfdr <- data.frame(Estimate = c(as.vector(res.lassomfdr[1:60,]), as.vector(res.lassomfdr[61:600,])),
                      Noise = c(rep(0, 60*ncol(res.lassomfdr)),
                                rep(1, 540*ncol(res.lassomfdr))))    

df.uni <-  data.frame(Estimate = c(as.vector(res.uni[1:60,]), as.vector(res.uni[61:600,])),
                      Noise = c(rep(0, 60*ncol(res.uni)),
                                rep(1, 540*ncol(res.uni))))     

df.uniash <-  data.frame(Estimate = c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])),
                         Noise = c(rep(0, 60*ncol(res.uniash)),
                                   rep(1, 540*ncol(res.uniash))))     

### Loess curves
xx <- seq(0, 1, by = .01)

l.cv <- with(df.cv, loess(Noise ~ Estimate, span = .2))
yy.cv <- predict(l.cv, newdata = xx)
smooth.cv <- data.frame(xx = xx, yy = yy.cv)

l.cv1se <- with(df.cv1se, loess(Noise ~ Estimate, span = .15))
yy.cv1se <- predict(l.cv1se, newdata = xx)
smooth.cv1se <- data.frame(xx = xx, yy = yy.cv1se)

l.mfdr <- with(df.mfdr, loess(Noise ~ Estimate, span = .2))
yy.mfdr <- predict(l.mfdr, newdata = xx)
smooth.mfdr <- data.frame(xx = xx, yy = yy.mfdr)

l.uni <- with(df.uni, loess(Noise ~ Estimate, span = .55))
yy.uni <- predict(l.uni, newdata = xx)
smooth.uni <- data.frame(xx = xx, yy = yy.uni)

l.uniash <- with(df.uniash, loess(Noise ~ Estimate, span = .2))
yy.uniash <- predict(l.uniash, newdata = xx)
smooth.uniash <- data.frame(xx = xx, yy = yy.uniash)

## There are too many points, only plot a sample for the scatterplot
sample.id <- sample(1:nrow(df.cv), size = 8000, replace = FALSE)

## Create scatterplots from random subset of points
p.cv2 <- ggplot(df.cv[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.cv1se2 <- ggplot(df.cv1se[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.mfdr2 <- ggplot(df.mfdr[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.uni2 <- ggplot(df.uni[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.uniash2 <- ggplot(df.uniash[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)


## Add loess curves from full sim results
p.cv2 <- p.cv2 + geom_path(data = smooth.cv, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1) +
  labs(y = "Proportion (noise)")

p.cv1se2 <- p.cv1se2 + geom_path(data = smooth.cv1se, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1) +
  labs(y = "")

p.mfdr2 <- p.mfdr2 + geom_path(data = smooth.mfdr, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1) +
  labs(y = "")

p.uni2 <- p.uni2 + geom_path(data = smooth.uni, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "")

p.uniash2 <- p.uniash2 + geom_path(data = smooth.uniash, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "")


#grid.arrange(p.cv, p.cv1se, p.mfdr, p.uni, p.uniash, 
#             p.cv2, p.cv1se2, p.mfdr2, p.uni2, p.uniash2, nrow = 2)

### Generate and save Figure2 as a pdf
png("figures_tables/Fig2.png", h=6, w=8.5, units = 'in', res = 300)

grid.arrange(arrangeGrob(p.cv, top="lasso mfdr (CV)"), arrangeGrob(p.cv1se, top="lasso mfdr (1SE)"),  
             arrangeGrob(p.uniash, top="Univariate fdr", right = "Assumptions Violated"),
             arrangeGrob(p.cv2), arrangeGrob(p.cv1se2),  
             arrangeGrob(p.uniash2, right = "Assumptions Met"), ncol=3)

dev.off()


########################################################################################
#
#           Code to produce Table 1 (requires simulation results from earlier)
#
########################################################################################


####### Row 1, linear assumptions violated
load(file = "simres/linear_hard.RData")

## Setup results data.frames
df.cv <- data.frame(Estimate = cut(c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])),
                                   breaks = c(0,.2,.4,.6,.8,1)),
                    Noise = c(rep(0, 60*ncol(res.lassocv)),
                              rep(1, 540*ncol(res.lassocv))))            

df.cv1se <- data.frame(Estimate = cut(c(as.vector(res.lassocv1se[1:60,]), as.vector(res.lassocv1se[61:600,])),
                                      breaks = c(0,.2,.4,.6,.8,1)),
                       Noise = c(rep(0, 60*ncol(res.lassocv1se)),
                                 rep(1, 540*ncol(res.lassocv1se))))     

df.mfdr <- data.frame(Estimate = cut(c(as.vector(res.lassomfdr[1:60,]), as.vector(res.lassomfdr[61:600,])),
                                     breaks = c(0,.2,.4,.6,.8,1.1)),
                      Noise = c(rep(0, 60*ncol(res.lassomfdr)),
                                rep(1, 540*ncol(res.lassomfdr))))    

df.uni <-  data.frame(Estimate = cut(c(as.vector(res.uni[1:60,]), as.vector(res.uni[61:600,])),
                                     breaks = c(0,.2,.4,.6,.8,1)),
                      Noise = c(rep(0, 60*ncol(res.uni)),
                                rep(1, 540*ncol(res.uni))))     

df.uniash <-  data.frame(Estimate = cut(c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])),
                                        breaks = c(0,.2,.4,.6,.8,1)),
                         Noise = c(rep(0, 60*ncol(res.uniash)),
                                   rep(1, 540*ncol(res.uniash))))     

props.cv <- df.cv %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.cv1se <- df.cv1se %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.mfdr <- df.mfdr %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.uni <- df.uni %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.uniash <- df.uniash %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))

tab.res <- data.frame(cv = props.cv$prop_noise,
                      cv1se = props.cv1se$prop_noise,
                      mfdr = props.mfdr$prop_noise,
                      uni = props.uni$prop_noise,
                      uniash = props.uniash$prop_noise)

rownames(tab.res) <- levels(props.cv$Estimate)
lin_hard_tab <- round(t(tab.res),3)
rownames(lin_hard_tab) <- paste("Linear", rownames(lin_hard_tab))

#### Row 2, logistic
load(file = "simres/logistic_hard.RData")

## Setup results data.frames
df.cv <- data.frame(Estimate = cut(c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])),
                                   breaks = c(0,.2,.4,.6,.8,1)),
                    Noise = c(rep(0, 60*ncol(res.lassocv)),
                              rep(1, 540*ncol(res.lassocv))))            

df.cv1se <- data.frame(Estimate = cut(c(as.vector(res.lassocv1se[1:60,]), as.vector(res.lassocv1se[61:600,])),
                                      breaks = c(0,.2,.4,.6,.8,1)),
                       Noise = c(rep(0, 60*ncol(res.lassocv1se)),
                                 rep(1, 540*ncol(res.lassocv1se))))     

df.mfdr <- data.frame(Estimate = cut(c(as.vector(res.lassomfdr[1:60,]), as.vector(res.lassomfdr[61:600,])),
                                     breaks = c(0,.2,.4,.6,.8,1.1)),
                      Noise = c(rep(0, 60*ncol(res.lassomfdr)),
                                rep(1, 540*ncol(res.lassomfdr))))    

df.uni <-  data.frame(Estimate = cut(c(as.vector(res.uni[1:60,]), as.vector(res.uni[61:600,])),
                                     breaks = c(0,.2,.4,.6,.8,1)),
                      Noise = c(rep(0, 60*ncol(res.uni)),
                                rep(1, 540*ncol(res.uni))))     

df.uniash <-  data.frame(Estimate = cut(c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])),
                                        breaks = c(0,.2,.4,.6,.8,1)),
                         Noise = c(rep(0, 60*ncol(res.uniash)),
                                   rep(1, 540*ncol(res.uniash))))     

props.cv <- df.cv %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.cv1se <- df.cv1se %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.mfdr <- df.mfdr %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.uni <- df.uni %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.uniash <- df.uniash %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))

tab.res <- data.frame(cv = props.cv$prop_noise,
                      cv1se = props.cv1se$prop_noise,
                      mfdr = props.mfdr$prop_noise,
                      uni = props.uni$prop_noise,
                      uniash = props.uniash$prop_noise)

rownames(tab.res) <- levels(props.cv$Estimate)
log_hard_tab <- round(t(tab.res),3)
rownames(log_hard_tab) <- paste("Logistic", rownames(log_hard_tab))


######## Row 3 - Cox 
load(file = "simres/cox_hard.RData")

## Setup results data.frames
df.cv <- data.frame(Estimate = cut(c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])),
                                   breaks = c(0,.2,.4,.6,.8,1)),
                    Noise = c(rep(0, 60*ncol(res.lassocv)),
                              rep(1, 540*ncol(res.lassocv))))            

df.cv1se <- data.frame(Estimate = cut(c(as.vector(res.lassocv1se[1:60,]), as.vector(res.lassocv1se[61:600,])),
                                      breaks = c(0,.2,.4,.6,.8,1)),
                       Noise = c(rep(0, 60*ncol(res.lassocv1se)),
                                 rep(1, 540*ncol(res.lassocv1se))))     

df.mfdr <- data.frame(Estimate = cut(c(as.vector(res.lassomfdr[1:60,]), as.vector(res.lassomfdr[61:600,])),
                                     breaks = c(0,.2,.4,.6,.8,1.1)),
                      Noise = c(rep(0, 60*ncol(res.lassomfdr)),
                                rep(1, 540*ncol(res.lassomfdr))))    

df.uni <-  data.frame(Estimate = cut(c(as.vector(res.uni[1:60,]), as.vector(res.uni[61:600,])),
                                     breaks = c(0,.2,.4,.6,.8,1)),
                      Noise = c(rep(0, 60*ncol(res.uni)),
                                rep(1, 540*ncol(res.uni))))     

df.uniash <-  data.frame(Estimate = cut(c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])),
                                        breaks = c(0,.2,.4,.6,.8,1)),
                         Noise = c(rep(0, 60*ncol(res.uniash)),
                                   rep(1, 540*ncol(res.uniash))))     

props.cv <- df.cv %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.cv1se <- df.cv1se %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.mfdr <- df.mfdr %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.uni <- df.uni %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))
props.uniash <- df.uniash %>% group_by(Estimate) %>% summarize(prop_noise = mean(Noise))

tab.res <- data.frame(cv = props.cv$prop_noise,
                      cv1se =props.cv1se$prop_noise, 
                      mfdr = props.mfdr$prop_noise,
                      uni = props.uni$prop_noise,
                      uniash = props.uniash$prop_noise)



rownames(tab.res) <- levels(props.cv$Estimate)
cox_hard_tab <- round(t(tab.res),3)
rownames(cox_hard_tab) <- paste("Cox", rownames(cox_hard_tab))

ttab1 <- rbind(lin_hard_tab, log_hard_tab, cox_hard_tab)
ttab1 <- ttab1[which(rownames(ttab1) %in% c("Linear cv", "Linear cv1se", "Linear uniash",
                                            "Logistic cv", "Logistic cv1se", "Logistic uniash",
                                            "Cox cv", "Cox cv1se", "Cox uniash")), ]

png("figures_tables/Table1.png", h=5, w=7, units = 'in', res = 300)
grid.table(ttab1)
dev.off()


########################################################################################
#
#           Code to produce Figure 3 (requires simulation results from earlier)
#
########################################################################################

load("simres/linear_hard.RData")

## Ids for A/B vars within first 60 (last 540 are all "C")
AB.ids <- rep("B", 60)
AB.ids[id.A] <- "A"

## Setup results data.frames
df.cv <- data.frame(Estimate = (c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])) < .1),
                    Type = c(rep(AB.ids, ncol(res.lassocv)),
                              rep("C", 540*ncol(res.lassocv))))            

df.cv.n <- group_by(df.cv, Type) %>% summarise(N = sum(Estimate)/ncol(res.lassocv))

df.cv1se <-  data.frame(Estimate = (c(as.vector(res.lassocv1se[1:60,]), as.vector(res.lassocv1se[61:600,])) < .1),
                        Type = c(rep(AB.ids, ncol(res.lassocv1se)),
                                 rep("C", 540*ncol(res.lassocv1se))))            

df.cv1se.n <- group_by(df.cv1se, Type) %>% summarise(N = sum(Estimate)/ncol(res.lassocv1se))


df.mfdr <- data.frame(Estimate = (c(as.vector(res.lassomfdr[1:60,]), as.vector(res.lassomfdr[61:600,])) < .1),
                     Type = c(rep(AB.ids, ncol(res.lassomfdr)),
                              rep("C", 540*ncol(res.lassomfdr))))            

df.mfdr.n <- group_by(df.mfdr, Type) %>% summarise(N = sum(Estimate)/ncol(res.lassomfdr))


df.uni <-  data.frame(Estimate = (c(as.vector(res.uni[1:60,]), as.vector(res.uni[61:600,])) < .1),
                      Type = c(rep(AB.ids, ncol(res.uni)),
                               rep("C", 540*ncol(res.uni))))            

df.uni.n <- group_by(df.uni, Type) %>% summarise(N = sum(Estimate)/ncol(res.uni))


df.uniash <-  data.frame(Estimate = (c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])) < .1),
                         Type = c(rep(AB.ids, ncol(res.uniash)),
                                  rep("C", 540*ncol(res.uniash))))            

df.uniash.n <- group_by(df.uni, Type) %>% summarise(N = sum(Estimate)/ncol(res.uniash))

power.res.hard <- data.frame(Val = c(df.cv.n$N, df.cv1se.n$N, df.mfdr.n$N, df.uni.n$N, df.uniash.n$N),
                        Method = c(rep("lasso mfdr (CV)", 3), rep("lasso mfdr (1SE)", 3), rep("Lasso (mFDR)", 3),
                                   rep("Univariate (locfdr)", 3), rep("Univariate fdr", 3)),
                        Var = c("Causal (A)","Correlated (B)","Noise (C)"), Scenario = "Assumptions Violated")


### AUC (assumptions Violated)
cvs <- c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])) 
est.cv.v <- cvs[which(df.cv$Type %in% c('A','C'))]
tr.cv.v <- ifelse(df.cv[which(df.cv$Type %in% c('A','C')),2] == 'C', 1, 0)
auc.cv.v <- pROC::auc(tr.cv.v, est.cv.v)

us <- c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])) 
est.uni.v <- us[which(df.uniash$Type %in% c('A','C'))]
tr.uni.v <- ifelse(df.uniash[which(df.uniash$Type %in% c('A','C')),2] == 'C', 1, 0)
auc.uni.v <- pROC::auc(tr.uni.v, est.uni.v)


load("simres/linear_easy.RData")
## Ids for A/B vars within first 60 (last 540 are all "C")
AB.ids <- rep("A", 60)

## Setup results data.frames
df.cv <- data.frame(Estimate = (c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])) < .1),
                    Type = c(rep(AB.ids, ncol(res.lassocv)),
                             rep("C", 540*ncol(res.lassocv))))            

df.cv.n <- group_by(df.cv, Type) %>% summarise(N = sum(Estimate)/ncol(res.lassocv))

df.cv1se <-  data.frame(Estimate = (c(as.vector(res.lassocv1se[1:60,]), as.vector(res.lassocv1se[61:600,])) < .1),
                        Type = c(rep(AB.ids, ncol(res.lassocv1se)),
                                 rep("C", 540*ncol(res.lassocv1se))))            

df.cv1se.n <- group_by(df.cv1se, Type) %>% summarise(N = sum(Estimate)/ncol(res.lassocv1se))


df.mfdr <- data.frame(Estimate = (c(as.vector(res.lassomfdr[1:60,]), as.vector(res.lassomfdr[61:600,])) < .1),
                      Type = c(rep(AB.ids, ncol(res.lassomfdr)),
                               rep("C", 540*ncol(res.lassomfdr))))            

df.mfdr.n <- group_by(df.mfdr, Type) %>% summarise(N = sum(Estimate)/ncol(res.lassomfdr))


df.uni <-  data.frame(Estimate = (c(as.vector(res.uni[1:60,]), as.vector(res.uni[61:600,])) < .1),
                      Type = c(rep(AB.ids, ncol(res.uni)),
                               rep("C", 540*ncol(res.uni))))            

df.uni.n <- group_by(df.uni, Type) %>% summarise(N = sum(Estimate)/ncol(res.uni))


df.uniash <-  data.frame(Estimate = (c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])) < .1),
                         Type = c(rep(AB.ids, ncol(res.uniash)),
                                  rep("C", 540*ncol(res.uniash))))            

df.uniash.n <- group_by(df.uni, Type) %>% summarise(N = sum(Estimate)/ncol(res.uniash))

power.res.easy <- data.frame(Val = c(df.cv.n$N, df.cv1se.n$N, df.mfdr.n$N, df.uni.n$N, df.uniash.n$N),
                             Method = c(rep("lasso mfdr (CV)", 2), rep("lasso mfdr (1SE)", 2), rep("Lasso (mFDR)", 2),
                                        rep("Univariate (locfdr)", 2), rep("Univariate fdr", 2)),
                             Var = c("Causal (A)","Noise (C)"), Scenario = "Assumptions Met")

### AUC (assumptions Met)
cvs <- c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])) 
est.cv.m <- cvs[which(df.cv$Type %in% c('A','C'))]
tr.cv.m <- ifelse(df.cv[which(df.cv$Type %in% c('A','C')),2] == 'C', 1, 0)
auc.cv.m <- pROC::auc(tr.cv.m, est.cv.m)

us <- c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])) 
est.uni.m <- us[which(df.uniash$Type %in% c('A','C'))]
tr.uni.m <- ifelse(df.uniash[which(df.uniash$Type %in% c('A','C')),2] == 'C', 1, 0)
auc.uni.m <- pROC::auc(tr.uni.m, est.uni.m)

### Power plot (Fig 3)

power.combined <- rbind(power.res.hard, power.res.easy)
power.combined <- power.combined[which(power.combined$Method %in% c("lasso mfdr (CV)", "lasso mfdr (1SE)", "Univariate fdr")),]

png("figures_tables/Fig3.png", h=4, w=6, units = 'in', res = 300)
p <- ggplot(power.combined, aes(x = relevel(Method, ref = "lasso mfdr (CV)"), y = Val, fill = Var)) + ggtitle("") +
  geom_bar(stat = "identity", position=position_stack(reverse=TRUE)) + scale_fill_manual(name = "Feature Type", values = pal(3)[c(2,3,1)]) +
  coord_flip() + scale_y_continuous(expression(paste("Features with ", widehat(mfdr), " or ", widehat(fdr), " < 0.10"))) + scale_x_discrete(" ")
p + facet_grid(.~ Scenario, scales = "free") + theme(axis.text.y = element_text(angle = 45, hjust = 1))
dev.off()

## Numbers appearing in the text accompanying Fig 3
power.combined

### Assumptions Violated
3.625/3.448 ## 5% more A vars found
3.760/0.105 ## Univariate declares significance for 35.8 times as many B vars
0.465/0.070 ## Univariate declares significance for 6.6 times as many C vars

### Assumptions Met
17.975/11.180  ## 61% more A vars in

### AUC results
auc.cv.m
auc.uni.m
auc.cv.v
auc.uni.v

########################################################################################
#
#           Simulations discussed in Sec 4.3 (comparisons to Selective Inf and SS)
#
########################################################################################

### Seed
set.seed(12345)

######## Assumptions met steup
n <- 1000
p <- 600

sig <- sqrt(n)
t <- 60
beta <- c(rep(4,t/2),rep(-4,t/2),rep(0,p-t))
nreps <- 200


fdr.true <- fdr.true.cv <- fdr.true.uni <- matrix(NA, nrow = nreps, ncol = p)
tres.pv <- fres.pv <- tres.spacing <- fres.spacing <- tres.modspacing <- tres.knockoff <- NULL
fres.modspacing <- tres.covtest <- fres.covtest <- tres.mfdr <- fres.mfdr <- fres.knockoff <- NULL
tres.ss <- cres.ss <- frate.ss <- numeric(nreps)
true.cv <- false.cv <- numeric(nreps)

for (i in 1:nreps){
  ### Generate Data
  X <- matrix(rnorm(n*p), ncol = p, nrow = n)
  y <- X %*% beta + rnorm(n,sd = sig)
  
  ### Use Sample Spliting
  ss.res <- hdi(X,y)
  ss.pvs <- fdr.adjust(ss.res$pval.corr)
  tres.ss[i] <- sum(ss.pvs[1:t] < .1)
  frate.ss[i] <- sum(ss.pvs[(t+1):p] < .1)/sum(ss.pvs < .1)
  
  
  ### Use selective inf
  fit.lar <- tryCatch(lar(X,y, maxsteps = 100), error=function(e) NULL)
  sig.est <- tryCatch(estimateSigma(X, y), error=function(e) NULL)
  res <- tryCatch(larInf(fit.lar, sigma = sig.est$sigmahat, k = 100), error=function(e) NULL)
  
  
  if (!is.null(res)){
    step.pv <- forwardStop(res$pv)
    step.spacing <- forwardStop(res$pv.spacing)
    step.modspacing <- forwardStop(res$pv.modspac)
    step.covtest <- forwardStop(res$pv.covtest)
    
    ### Results for each
    if (step.pv == 0){
      true.pv <- 0
    } else {
      true.pv <- sum(res$vars[1:step.pv] <= t)
    }
    false.pv <- sum(res$vars[1:step.pv] > t)
    frate.pv <- false.pv/step.pv
    ### Spacing test 
    if (step.spacing == 0){
      true.spacing <- 0
    } else {
      true.spacing <- sum(res$vars[1:step.spacing] <= t)
    }
    false.spacing <- sum(res$vars[1:step.spacing] > t)
    frate.spacing <- false.spacing/step.spacing 
    ### mod spacing test
    if (step.modspacing == 0){
      true.modspacing <- 0
    } else {
      true.modspacing <- sum(res$vars[1:step.modspacing] <= t)
    }
    false.modspacing <- sum(res$vars[1:step.modspacing] > t)
    frate.modspacing <- false.modspacing/step.modspacing 
    ### Covtest  
    if (step.covtest == 0){
      true.covtest <- 0
    } else {
      true.covtest <- sum(res$vars[1:step.covtest] <= t)
    }
    false.covtest <- sum(res$vars[1:step.covtest] > t)
    frate.covtest <- false.covtest/step.covtest
  } else {
    true.covtest <- frate.covtest <- frate.modspacing <- true.modspacing <- true.spacing <- frate.spacing <- true.pv <- frate.pv <- NA
  }  
  
  ### Knock-off filter
  suppressWarnings(knres <- knockoff.filter(X = X, y = y, fdr = .1))
  true.kn <- sum(knres$selected <= t) # true selections
  false.kn <-sum(knres$selected > t) # false selections
  frate.kn <- false.kn/length(knres$selected) # false rate
  
  ### store iteration results
  tres.pv <- c(tres.pv, true.pv)
  fres.pv <- c(fres.pv, frate.pv)
  
  tres.spacing <- c(tres.spacing, true.spacing)
  fres.spacing <- c(fres.spacing, frate.spacing)
  
  tres.modspacing <- c(tres.modspacing, true.modspacing)
  fres.modspacing <- c(fres.modspacing, frate.modspacing)
  
  tres.covtest <- c(tres.covtest, true.covtest)
  fres.covtest <- c(fres.covtest, frate.covtest)
  
  tres.knockoff <- c(tres.knockoff, true.kn)
  fres.knockoff <- c(fres.knockoff ,frate.kn)
}

final.true <- data.frame(pv = tres.pv, spacing = tres.spacing, modspacing = tres.modspacing, covtest = tres.covtest, ss = tres.ss, knock = tres.knockoff)
final.false <- data.frame(pv = fres.pv, spacing = fres.spacing, modspacing = fres.modspacing, covtest = fres.covtest, ss = frate.ss, knock = fres.knockoff)

save(final.true, final.false, file = "simres/model_easy.RData")


################### Assumptions Violated ###############################

### Seed
set.seed(12345)

### Setup
n <- 200
p <- 600
bb <- numeric(60)
bb[(0:5)*10+1] <- c(6,-6,5,-5,4,-4)  #### 6 true variables (type A)

id.A <- which(bb != 0)
id.B <- 1:60
id.B <- id.B[-id.A]
id.C <- 61:p

t <- 60
nreps <- 200

fdr.true <- fdr.true.cv <- fdr.true.uni <- matrix(NA, nrow = nreps, ncol = p)
tres.pv <- fres.pv <- tres.spacing <- fres.spacing <- tres.modspacing <- tres.knockoff <- NULL
cres.pv <- cres.pv <- cres.spacing <- cres.covtest <- cres.modspacing <- cres.knockoff <- NULL
fres.modspacing <- tres.covtest <- fres.covtest <- tres.mfdr <- fres.mfdr <- fres.knockoff <- NULL
tres.ss <- cres.ss <- frate.ss <- numeric(nreps)
true.cv <- cor.cv <- false.cv <- numeric(nreps)

for (i in 1:nreps){
  ### Corr data
  D1 <- genData(n, J=6, J0=6, K=10, K0=1, rho=0, rho.g=0.5, beta=bb)  #### 9 correlated vars for each true (type B)
  D2 <- genData(n, p - 60, rho=0.8, beta=0, corr="auto")
  
  X <- cbind(D1$X, D2$X)
  y <- D1$y
  
  ### Use Sample Spliting
  ss.res <- hdi(X,y)
  ss.pvs <- fdr.adjust(ss.res$pval.corr)
  tres.ss[i] <- sum(ss.pvs[id.A] < .1)
  cres.ss[i] <- sum(ss.pvs[id.B] < .1)
  frate.ss[i] <- sum(ss.pvs[id.C] < .1)/sum(ss.pvs < .1)
  
  ### Use selective inf
  fit.lar <- tryCatch(lar(X,y, maxsteps = 105), error=function(e) NULL)
  sig.est <- tryCatch(estimateSigma(X, y), error=function(e) NULL)
  res <- tryCatch(larInf(fit.lar, sigma = sig.est$sigmahat, k = 100), error=function(e) NULL)
  
  
  if (!is.null(res)){
    step.pv <- forwardStop(res$pv)
    step.spacing <- forwardStop(res$pv.spacing)
    step.modspacing <- forwardStop(res$pv.modspac)
    step.covtest <- forwardStop(res$pv.covtest)
    
    ### Results for each
    if (step.pv == 0){
      true.pv <- 0
      cor.pv <- 0
    } else {
      true.pv <- sum(res$vars[1:step.pv] %in% id.A)
      cor.pv <- sum(res$vars[1:step.pv] %in% id.B)
    }
    false.pv <- sum(res$vars[1:step.pv] > t)
    frate.pv <- false.pv/step.pv
    
    ### Spacing test 
    if (step.spacing == 0){
      true.spacing <- 0
      cor.spacing <- 0
    } else {
      true.spacing <- sum(res$vars[1:step.spacing] %in% id.A)
      cor.spacing <- sum(res$vars[1:step.spacing] %in% id.B)
    }
    false.spacing <- sum(res$vars[1:step.spacing] > t)
    frate.spacing <- false.spacing/step.spacing 
    
    ### mod spacing test
    if (step.modspacing == 0){
      true.modspacing <- 0
      cor.modspacing <- 0
    } else {
      true.modspacing <- sum(res$vars[1:step.modspacing] %in% id.A)
      cor.modspacing <- sum(res$vars[1:step.modspacing] %in% id.B)
    }
    false.modspacing <- sum(res$vars[1:step.modspacing] > t)
    frate.modspacing <- false.modspacing/step.modspacing 
    ### Covtest  
    if (step.covtest == 0){
      true.covtest <- 0
      cor.covtest <- 0
    } else {
      true.covtest <- sum(res$vars[1:step.covtest] %in% id.A)
      cor.covtest <- sum(res$vars[1:step.covtest] %in% id.B)
    }
    false.covtest <- sum(res$vars[1:step.covtest] > t)
    frate.covtest <- false.covtest/step.covtest
  } else {
    true.covtest <- cor.covtest <- frate.covtest <- NA
    cor.pv <- cor.modspacing <- frate.modspacing <- true.modspacing <- true.spacing <- NA
    cor.spacing <- frate.spacing <- true.pv <- frate.pv <- NA
  }  
  
  ### Knock-off Filter
  suppressWarnings(knres <- knockoff.filter(X = X, y = y, fdr = .1))
  true.kn <- sum(knres$selected %in% id.A) # var A selections
  cor.kn <- sum(knres$selected %in% id.B) # var B selections
  false.kn <- sum(knres$selected > t) # var C selections
  frate.kn <- false.kn/length(knres$selected) # false rate
  
  ### store iteration results
  tres.pv <- c(tres.pv, true.pv)
  cres.pv <- c(cres.pv, cor.pv)
  fres.pv <- c(fres.pv, frate.pv)
  
  tres.spacing <- c(tres.spacing, true.spacing)
  cres.spacing <- c(cres.spacing, cor.spacing)
  fres.spacing <- c(fres.spacing, frate.spacing)
  
  tres.modspacing <- c(tres.modspacing, true.modspacing)
  cres.modspacing <- c(cres.modspacing, cor.modspacing)
  fres.modspacing <- c(fres.modspacing, frate.modspacing)
  
  tres.covtest <- c(tres.covtest, true.covtest)
  cres.covtest <- c(cres.covtest, cor.covtest)
  fres.covtest <- c(fres.covtest, frate.covtest)
  
  tres.knockoff <- c(tres.knockoff, true.kn)
  cres.knockoff <- c(cres.knockoff, cor.kn)
  fres.knockoff <- c(fres.knockoff, frate.kn)
  
}

final.true <- data.frame(pv = tres.pv, spacing = tres.spacing, modspacing = tres.modspacing, covtest = tres.covtest, ss = tres.ss, knock = tres.knockoff)
final.cor <- data.frame(pv = cres.pv, spacing = cres.spacing, modspacing = cres.modspacing, covtest = cres.covtest, ss = cres.ss, knock = cres.knockoff)
final.false <- data.frame(pv = fres.pv, spacing = fres.spacing, modspacing = fres.modspacing, covtest = fres.covtest, ss = frate.ss, knock = fres.knockoff)

save(id.A, id.B, id.C, final.true, final.cor, final.false, file = "simres/model_hard.RData")

########################################################################################
#
#           Code to produce Table 2 (requires simulation results from earlier)
#
########################################################################################



load(file = "simres/model_easy.RData")
p <- ncol(fdr.true)

selinf.A <- apply(final.true, 2, mean)
final.false[is.na(final.false)] <- 0  ### Fdr of 0/0 is defined as 0
selinf.nC <- apply(final.false, 2, mean)
selinf.perC <- colSums(final.false)/(colSums(final.true) + colSums(final.false))

load(file = "simres/linear_easy.RData")
locmfdr.A <- sum(rowMeans(res.lassocv[1:60,] < .1))
locmfdr.nC <- sum(rowMeans(res.lassocv[61:600,] < .1))
locmfdr.perC <- sum(rowMeans(res.lassocv[61:600,] < .1))/sum(rowMeans(res.lassocv < .1))

mod.res.easy <- data.frame(rbind(c(locmfdr.A, selinf.A),
                       NA,
                       c(locmfdr.nC, selinf.nC),
                       c(locmfdr.perC, selinf.perC)))

rownames(mod.res.easy) <- c("Avg 'A'", "Avg 'B'", "Avg C", "Avg FDR")
colnames(mod.res.easy) <- c("loc-mfdr", colnames(mod.res.easy)[-1])


load(file = "simres/model_hard.RData")
selinf.A <- apply(final.true, 2, mean)
selinf.B <- apply(final.cor, 2, mean)
final.false[is.na(final.false)] <- 0 ### Fdr of 0/0 is defined as 0
selinf.nC <- apply(final.false, 2, mean)
selinf.perC <- colSums(final.false)/(colSums(final.true) + colSums(final.false))

load(file = "simres/linear_hard.RData")
locmfdr.A <- sum(rowMeans(res.lassocv[id.A,] < .1))
locmfdr.B <- sum(rowMeans(res.lassocv[id.B,] < .1))
locmfdr.nC <- sum(rowMeans(res.lassocv[61:600,] < .1))
locmfdr.perC <- sum(rowMeans(res.lassocv[61:600,] < .1))/sum(rowMeans(res.lassocv< .1))

mod.res.hard <- data.frame(rbind(c(locmfdr.A,selinf.A),
                        c(locmfdr.B,selinf.B),
                        c(locmfdr.nC,selinf.nC),
                        c(locmfdr.perC,selinf.perC)))

rownames(mod.res.hard) <- c("Avg 'A' - Violated", "Avg 'B' - Violated", "Avg C - Violated", "Avg FDR - Violated")
colnames(mod.res.hard) <- c("loc-mfdr", colnames(mod.res.hard)[-1])


png("figures_tables/Table2.png", h=7, w=10, units = 'in', res = 300)
grid.table(round(t(rbind(mod.res.hard, mod.res.easy))[,-6], 2))
dev.off()


########################################################################################
#
#        Code used to produce case study results (Table 3, Figure 4, Table 4)
#
########################################################################################


### Case study #1 (table 3, figure 4)

## Load Shedden dataset
load("case_study_data/Shedden2008.RData")

### Set up model matrix
ZZ <- model.matrix(~ factor(Sex) + factor(Race) + factor(AdjChemo) + factor(SmHist) + factor(Margin) + factor(Grade), Z)
XX <- cbind(X, ZZ)
w <- rep(0:1, c(ncol(ZZ), ncol(X))) ## Don't the penalize clinical covariates

### Fit model
set.seed(12345)
cv.fit <- cv.ncvsurv(XX, S, penalty = 'lasso', returnX = TRUE, penalty.factor = w)
fit <- cv.fit$fit

## positions of lambda values of interest
cv1se.lam <- min(which(cv.fit$cve - cv.fit$cvse <  min(cv.fit$cve)))
cv.lam <- which(cv.fit$fit$lambda == cv.fit$lambda.min)

## local mfdr using model/lambdas from above
#locfdr.cv  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "ashr", nullweight = 1, mixcompdist = "halfnormal")$pen.vars
#locfdr.cv1se <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "ashr", nullweight = 1, mixcompdist = "halfnormal")$pen.vars

locfdr.cv  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "kernel")$pen.vars
locfdr.cv1se <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "kernel")$pen.vars


sum(locfdr.cv$Selected == "*")  # 43 selections at lambda CV
sum(locfdr.cv1se$Selected == "*")  # 1 selections at lambda CV1se

## Results
shedden.cv.res <- data.frame(name = row.names(locfdr.cv[order(locfdr.cv$mfdr),])[1:10],
                          z = locfdr.cv[order(locfdr.cv$mfdr),][1:10,2],
                          fdr = locfdr.cv[order(locfdr.cv$mfdr),][1:10,3])


shedden.cv1se.res <- data.frame(name = row.names(locfdr.cv1se[order(locfdr.cv1se$mfdr),])[1:10],
                             z = locfdr.cv1se[order(locfdr.cv1se$mfdr),][1:10,2],
                             fdr = locfdr.cv1se[order(locfdr.cv1se$mfdr),][1:10,3])

#### Univariate locfdr (for comparison)
zstat <- rep(NA, ncol(X))
for (j in 1:ncol(X)){
  fit.uni <- coxph(S ~ X[,j] + ZZ[,-1])
  zstat[j] <- summary(fit.uni)$coefficients[1,4]
}

## Univariate using ashr
univariate.res.ash <- get_lfdr(ash(beta = zstat, se = rep(1, length(zstat))))

shedden.uni.res <- data.frame(name = colnames(X[,order(univariate.res.ash)[1:10]]),
                         z = zstat[order(univariate.res.ash)][1:10],
                         fdr = univariate.res.ash[order(univariate.res.ash)][1:10])

#mean(abs(cor(X)))

#############
### Table 3
#############

u.names <- fData[row.names(fData) %in% shedden.uni.res$name,]
u.names2 <- u.names[match(shedden.uni.res$name, rownames(u.names)),1]

cv.names <- fData[row.names(fData) %in% shedden.cv.res$name,]
cv.names2 <- cv.names[match(shedden.cv.res$name, rownames(cv.names)),1]

cv1se.names <- fData[row.names(fData) %in% shedden.cv1se.res$name,]
cv1se.names2 <- cv1se.names[match(shedden.cv1se.res$name, rownames(cv1se.names)),1]

## note: the <NA> names are missing in the fData file, they are replaced 
## with the probe_id in the table appearing in the manuscript

png("figures_tables/Table3.png", h=6, w=8.5, units = 'in', res = 300)
grid.table(data.frame(feature.uni = u.names2, fdr.uni = round(shedden.uni.res$fdr,4),
                       feature.cv1se = cv1se.names2, mfdr.cv1se = round(shedden.cv1se.res$fdr,4),
                       feature.cv = cv.names2, mfdr.cv = round(shedden.cv.res$fdr,4)))
dev.off()


#############
### Figure 4
#############

set.seed(12345)
dfd <- data.frame(z = c(locfdr.cv$z, locfdr.cv1se$z, zstat, rnorm(length(zstat))),
                  ID = c(rep("mfdr (CV)", length(locfdr.cv$z)),
                         rep("mfdr (1SE)", length(locfdr.cv1se$z)), 
                         rep("Univariate fdr", length(zstat)),
                         rep("Theoretical N(0,1)", length(zstat))),
                  lty = c(rep("1", length(locfdr.cv$z)),
                          rep("1", length(locfdr.cv1se$z)), 
                          rep("1", length(zstat)),
                          rep("2", length(zstat))))

dfd$ID <- relevel(dfd$ID, ref = "Theoretical N(0,1)")
myColors <- c("black", 'red', 'blue','green')
names(myColors) <- levels(dfd$ID)
colScale <- scale_colour_manual(name = "",values = myColors)

pd <- ggplot(data = dfd, aes(x = z, color = ID, lty = lty)) + geom_line(stat = "density", lwd = 1.1) + 
  labs(y = "Density", colour = "") + scale_linetype(guide = FALSE) + scale_x_continuous(breaks = c(-5:5),
                                                                                        limits = c(-5.5,5.5)) + colScale

png("figures_tables/Fig4.png", h=4.5, w=6, units = 'in', res = 300)
pd
dev.off()


### Case study #2 (table 4)

## Load TCGA BRCA1 dataset
TCGA <- readRDS("case_study_data//bcTCGA.rds")

## fit model
X <- TCGA$X
y <- TCGA$y
set.seed(12345)
cv.fit <- cv.ncvreg(X, y, penalty = "lasso")
fit <- cv.fit$fit

## positions of lambda values of interest
cv1se.lam <- min(which(cv.fit$cve - cv.fit$cvse <  min(cv.fit$cve)))
cv.lam <- which(cv.fit$fit$lambda == cv.fit$lambda.min)

## local mfdr using model/lambdas from above
locfdr.cv  <- ncvreg::local_mfdr(fit, fit$lambda[cv.lam], method = "ashr", mixcompdist = "halfuniform")
locfdr.cv1se <- ncvreg::local_mfdr(fit, fit$lambda[cv1se.lam], method = "ashr", mixcompdist = "halfuniform")

#summary(fit, fit$lambda[cv.lam])
#summary(fit, fit$lambda[cv1se.lam])
#cv.fit$cve[cv1se.lam]/cv.fit$cve[cv.lam]

### Univariate testing
tstat <- numeric(ncol(X))
for (j in 1:ncol(X)){
  fit.lm <- lm(y ~ X[,j])
  tstat[j] <- summary(fit.lm)$coefficients[2,3]
}
zstat <- qnorm(pt(tstat,length(y)- 2))
zstat[is.infinite(zstat)] <- tstat[is.infinite(zstat)] ## Make infinite z-stats into their pre-transformed t-stat

#univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)$fdr
univariate.res <- get_lfdr(ash(beta = zstat, se = rep(1, length(zstat))))

id <- order(univariate.res)[1:10]  ## Top 10 univariate selections
uni.output <- data.frame(name = colnames(X[,order(univariate.res)[1:10]]),
                         z = zstat[order(univariate.res)][1:10],
                         fdr = univariate.res[order(univariate.res)][1:10])


#############
### Table 4
#############

# Note: Chromosome locations are not contained in these data and were manually added to the table appearing in the manuscript

png("figures_tables/Table4.png", h=7, w=8, units = 'in', res = 300)
grid.table(data.frame(Gene = uni.output$name, uni_fdr = uni.output$fdr,
                       CV1se_mfdr = locfdr.cv1se[id,3],
                       CV_mfdr = locfdr.cv[id,3]))
dev.off()


#############
### Figure 5
#############

### Sim Seed
set.seed(145)
#set.seed(138)

### Setup
n <- 200
p <- 600
t <- 60
bb <- numeric(60)
bb[(0:5)*10+1] <- c(6,-6,5,-5,4,-4)  #### 6 true variables (type A)

## Ids of var types
id.A <- which(bb != 0)
id.B <- 1:t
id.B <- id.B[-id.A]
id.C <- (t+1):p

### Simulate correlated data
D1 <- genData(n, J=6, J0=6, K=10, K0=1, rho=0, rho.g=0.5, beta=bb)  #### 9 correlated vars for each true (type B)
D2 <- genData(n, p - 60, rho=0.8, beta=0, corr="auto")
X <- cbind(D1$X, D2$X)
y <- D1$y

### Fit model
cv.fit <- cv.ncvreg(X, y, penalty = "lasso")
fit <- cv.fit$fit

locfdr.cv.ashr  <- ncvreg::local_mfdr(fit, cv.fit$lambda.min, method = "ashr")
locfdr.cv.kernel  <- ncvreg::local_mfdr(fit, cv.fit$lambda.min, method = "kernel")


fun <- approxfun(density(locfdr.cv.kernel$z))

p1 <- ggplot(locfdr.cv.ashr, aes(x = z)) + geom_histogram(aes(y = ..density..), bins = 50, col = "black", fill = "grey")+ 
  stat_function(fun=dnorm,args=list(mean=0, sd=1), color = "black", size = 1, linetype = "dashed") + 
  stat_function(fun=fun,args=list(), color = "black", size = 1) + geom_rug() + 
  labs(x = "z-value", y = "Density", title = " ")

df <- rbind(data.frame(locfdr.cv.kernel, Method = "density"),
            data.frame(locfdr.cv.ashr, Method = "ashr"))

#p2 <- ggplot(locfdr.cv.kernel, aes(x = z, y = mfdr)) + geom_point() + labs(x = "z-value", y = "mfdr Estimate", title = "Kernel")
#p3 <- ggplot(locfdr.cv.ashr, aes(x = z, y = mfdr)) + geom_point() + labs(x = "z-value", y = "mfdr Estimate", title = "Ashr")

p4 <- ggplot() + geom_jitter(data = df, aes(x = z, y = mfdr, col = Method), width = .1, height = .01, alpha = .9)  + 
  geom_rug(data = df, aes(x = z)) + labs(x = "z-value", y = expression(widehat(mfdr)), title = " ") + theme(legend.position="top")

png("figures_tables/Fig5.png", h=5, w=8, units = 'in', res = 300)
grid.arrange(p1,p4, nrow = 1)
dev.off()



#ashr.fit <- ash.workhorse(locfdr.cv.kernel$z, rep(1, nrow(locfdr.cv.kernel)), mixcompdist = "normal")

     