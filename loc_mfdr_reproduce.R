####################################################################################################
#
#                                                                                     
#   Filename    :	loc_mfdr_reproduce.R										  
#   Input data files  :    ---                                                        
#   Output data files :    fig1.pdf, fig2.pdf, fig3.pdf, fig4.pdf, table1.pdf, 
#			                     table2.pdf, table3.pdf, table4.pdf 
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
library(Rccp)
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
  temp <- locmfdr.new(fit$fit, lambda = lseq[i], number = ncol(X))
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
  temp <- locmfdr.new(fit$fit, lambda = lseq[i], number = ncol(X))
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


## Assemble the 4 plots into a single figure and store as pdf
pdf("figures_tables/Fig1.pdf", height=7, width=8.5)
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

### Setup
n <- 1000
p <- 600
bb <- numeric(60)
bb[(0:5)*10+1] <- c(5.5,-5.5,5,-5,4.5,-4.5)

id.A <- which(bb != 0)
id.B <- 1:60
id.B <- id.B[-id.A]
id.C <- 61:p

t <- 60
beta <- c(rep(.15,t/2), rep(-.15,t/2), rep(0,p-t))
nreps <- 200
pi0 = 1

trueA <- trueB <- false <- trueA2 <- trueB2 <- false2 <- matrix(NA, nrow = nreps, ncol = 5)
trueA3 <- trueB3 <- false3 <-  matrix(NA, nrow = nreps, ncol = 5)
fdr.true <- fdr.true.cv <- fdr.true.uni <- matrix(NA, nrow = nreps, ncol = t)
multi.A <- multi.A2 <- multi.B <- multi.B2 <- multi.C <- multi.C2 <- rep(NA, nreps)
uni.A <- uni.B <- uni.C <- rep(NA,nreps)
res.lassomfdr <- res.lassocv <- res.uni <- matrix(NA, nrow = p, ncol = nreps)
for (i in 1:nreps){
  
  D1 <- genData(n, p, beta=beta, family = "survival")
  
  X <- D1$X
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
  # Uses locfdr package
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)
  
  
  ### Fit model
  X <- std(X)
  fit <- ncvsurv(X,y, penalty = "lasso", returnX = TRUE)
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv.fit <- cv.ncvsurv(X,y,penalty = "lasso", returnX = TRUE)
  cv.lam <- which(cv.fit$fit$lambda == cv.fit$lambda.min)
  
  locfdr.lam <- locmfdr(fit, fit$lambda[mfdr.lam])[,2]
  locfdr.cv <- locmfdr(cv.fit$fit, cv.fit$lambda.min)[,2]
  
  #### store fdr for true vars
  fdr.true[i,] <- locfdr.lam[1:t]
  fdr.true.cv[i,] <- locfdr.cv[1:t]
  fdr.true.uni[i,] <- univariate.res$fdr[1:t]
  
  
  ### Seperate true/false vars into fdr bins 
  bins <- cut(locfdr.lam, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA[i,] <- table(bins[1:t])
  trueB[i,] <- rep(0,5)
  false[i,] <- table(bins[(t+1):p])
  
  bins2 <- cut(locfdr.cv, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA2[i,] <- table(bins2[1:t])
  trueB2[i,] <- rep(0,5)
  false2[i,] <- table(bins2[(t+1):p])
  
  bins3 <- cut(univariate.res$fdr, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA3[i,] <- table(bins3[1:t])
  trueB3[i,] <- rep(0,5)
  false3[i,] <- table(bins3[(t+1):p])
  
  
  ##### Power using 10% threshold
  uni.A[i] <- sum(univariate.res$fdr[1:t] < .1)
  uni.C[i] <- sum(univariate.res$fdr[(t+1):p] < .1)
  
  multi.A[i] <- sum(locfdr.lam[1:t] < .1)
  multi.C[i] <- sum(locfdr.lam[(t+1):p] < .1)
  
  multi.A2[i] <- sum(locfdr.cv[1:t] < .1)
  multi.C2[i] <- sum(locfdr.cv[(t+1):p] < .1)
  
  res.lassomfdr[,i] <- locfdr.lam
  res.lassocv[,i] <- locfdr.cv
  res.uni[,i] <- univariate.res$fdr
}

save(trueA, trueA2, trueA3, trueB, trueB2, trueB3, false, false2, false3, bins, 
     uni.A, uni.B, uni.C, multi.A, multi.B, multi.C, multi.A2, multi.B2, multi.C2, res.lassomfdr,
     res.lassocv, res.uni, file = "simres/cox_easy.RData")


############################################
##### Cox reg "assumptions violated" scenario
############################################ 

### Setup
n <- 200
p <- 600
bb <- numeric(60)
bb[(0:5)*10+1] <- c(1, -1, .75, -.75, .5, -.5)  #### 6 true variables (type A)

id.A <- which(bb != 0)
id.B <- 1:60
id.B <- id.B[-id.A]
id.C <- 61:p

nsim <- 200

trueA <- trueB <- false <- trueA2 <- trueB2 <- false2 <- matrix(NA, nrow = nsim, ncol = 5)
trueA3 <- trueB3 <- false3 <-  matrix(NA, nrow = nsim, ncol = 5)
fdr.true <- fdr.true.cv <- fdr.true.uni <- matrix(NA, nrow = nsim, ncol = 6)
multi.A <- multi.A2 <- multi.B <- multi.B2 <- multi.C <- multi.C2 <- rep(NA, nsim)
uni.A <- uni.B <- uni.C <- rep(NA,nsim)
res.lassomfdr <- res.lassocv <- res.uni <- matrix(NA, nrow = p, ncol = nsim)

for(i in 1:nsim){
  D1 <- genData(n, J=6, J0=6, K=10, K0=1, rho=0, rho.g=0.5, beta=bb, family = "survival") 
  D2 <- genData(n, p - 60, beta=0)
  
  X <- cbind(D1$X, D2$X)
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
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)
  
  
  ### Fit model
  X <- std(X)
  fit <- ncvsurv(X,y, penalty = "lasso", returnX = TRUE)
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv.fit <- cv.ncvsurv(X,y,penalty = "lasso", returnX = TRUE)
  cv.lam <- which(cv.fit$fit$lambda == cv.fit$lambda.min)
  
  locfdr.lam <- locmfdr(fit, fit$lambda[mfdr.lam])[,2]
  locfdr.cv <- locmfdr(cv.fit$fit, cv.fit$lambda.min)[,2]
  
  #### store fdr for true vars
  fdr.true[i,] <- locfdr.lam[(0:5)*10+1]
  fdr.true.cv[i,] <- locfdr.cv[(0:5)*10+1]
  fdr.true.uni[i,] <- univariate.res$fdr[(0:5)*10+1]
  
  
  ### Seperate true/false vars into fdr bins 
  bins <- cut(locfdr.lam, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA[i,] <- table(bins[(0:5)*10+1])
  trueB[i,] <- table(bins[(1:60)[-((0:5)*10+1)]])
  false[i,] <- table(bins[61:p])
  
  bins2 <- cut(locfdr.cv, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA2[i,] <- table(bins2[(0:5)*10+1])
  trueB2[i,] <- table(bins2[(1:60)[-((0:5)*10+1)]])
  false2[i,] <- table(bins2[61:p])
  
  bins3 <- cut(univariate.res$fdr, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA3[i,] <- table(bins3[(0:5)*10+1])
  trueB3[i,] <- table(bins3[(1:60)[-((0:5)*10+1)]])
  false3[i,] <- table(bins3[61:p])
  
  
  ##### Power using 10% threshold
  uni.A[i] <- sum(univariate.res$fdr[id.A] < .1)
  uni.B[i] <- sum(univariate.res$fdr[id.B] < .1)
  uni.C[i] <- sum(univariate.res$fdr[id.C] < .1)
  
  multi.A[i] <- sum(locfdr.lam[id.A] < .1)
  multi.B[i] <- sum(locfdr.lam[id.B] < .1)
  multi.C[i] <- sum(locfdr.lam[id.C] < .1)
  
  multi.A2[i] <- sum(locfdr.cv[id.A] < .1)
  multi.B2[i] <- sum(locfdr.cv[id.B] < .1)
  multi.C2[i] <- sum(locfdr.cv[id.C] < .1)
  
  res.lassomfdr[,i] <- locfdr.lam
  res.lassocv[,i] <- locfdr.cv
  res.uni[,i] <- univariate.res$fdr
}

## Save results
save(trueA, trueA2, trueA3, trueB, trueB2, trueB3, false, false2, false3, bins, uni.A, uni.B, uni.C, multi.A, multi.B, multi.C, multi.A2, multi.B2, multi.C2, 
     res.lassomfdr, res.lassocv, res.uni, file = "simres/cox_hard.RData")




####################################################
##### Logistic reg "assumptions met" scenario
####################################################


### Setup
n <- 1000
p <- 600
bb <- numeric(60)
bb[(0:5)*10+1] <- c(5.5,-5.5,5,-5,4.5,-4.5)

id.A <- which(bb != 0)
id.B <- 1:60
id.B <- id.B[-id.A]
id.C <- 61:p

sig <- 1
t <- 60
beta <- c(rep(.10,t/2), rep(-.10,t/2), rep(0,p-t))
nreps <- 200
pi0 = 1

trueA <- trueB <- false <- trueA2 <- trueB2 <- false2 <- matrix(NA, nrow = nreps, ncol = 5)
trueA3 <- trueB3 <- false3 <-  matrix(NA, nrow = nreps, ncol = 5)
fdr.true <- fdr.true.cv <- fdr.true.uni <- matrix(NA, nrow = nreps, ncol = t)
multi.A <- multi.A2 <- multi.B <- multi.B2 <- multi.C <- multi.C2 <- rep(NA, nreps)
uni.A <- uni.B <- uni.C <- rep(NA,nreps)
res.lassomfdr <- res.lassocv <- res.uni <- matrix(NA, nrow = p, ncol = nreps)

for (i in 1:nreps){
  D1 <- genData(n, p, beta=beta, family = "binomial")  
  
  X <- D1$X
  y <- D1$y
  
  ### Fit model
  X <- std(X)
  fit <- ncvreg(X,y, penalty = "lasso", family = "binomial", returnX = TRUE)
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv.fit <- cv.ncvreg(X,y,penalty = "lasso", family = "binomial")
  cv.lam <- which(fit$lambda == cv.fit$lambda.min)
  
  locfdr.lam <- locmfdr(fit, fit$lambda[mfdr.lam])[,2]
  locfdr.cv <- locmfdr(fit, cv.fit$lambda.min)[,2]
  
  
  ### Estimate univariate locfdr
  zstat <- tstat <- rep(NA, ncol(X))
  for (j in 1:ncol(X)){
    fit.lm <- glm(y ~ X[,j], family = "binomial")
    zstat[j] <- summary(fit.lm)$coefficients[2,3]
  } 
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)
  
  #### store fdr for true vars
  fdr.true[i,] <- locfdr.lam[1:t]
  fdr.true.cv[i,] <- locfdr.cv[1:t]
  fdr.true.uni[i,] <- univariate.res$fdr[1:t]
  
  
  ### Seperate true/false vars into fdr bins 
  bins <- cut(locfdr.lam, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA[i,] <- table(bins[1:t])
  trueB[i,] <- rep(0,5)
  false[i,] <- table(bins[(t+1):p])
  
  bins2 <- cut(locfdr.cv, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA2[i,] <- table(bins2[1:t])
  trueB2[i,] <- rep(0,5)
  false2[i,] <- table(bins2[(t+1):p])
  
  bins3 <- cut(univariate.res$fdr, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA3[i,] <- table(bins3[1:t])
  trueB3[i,] <- rep(0,5)
  false3[i,] <- table(bins3[(t+1):p])
  
  
  ##### Power using 10% threshold
  uni.A[i] <- sum(univariate.res$fdr[1:t] < .1)
  #uni.B[i] <- sum(univariate.res$fdr[id.B] < .1)
  uni.C[i] <- sum(univariate.res$fdr[(t+1):p] < .1)
  
  multi.A[i] <- sum(locfdr.lam[1:t] < .1)
  #multi.B[i] <- sum(locfdr.lam[id.B] < .1)
  multi.C[i] <- sum(locfdr.lam[(t+1):p] < .1)
  
  multi.A2[i] <- sum(locfdr.cv[1:t] < .1)
  #multi.B2[i] <- sum(locfdr.cv[id.B] < .1)
  multi.C2[i] <- sum(locfdr.cv[(t+1):p] < .1)
  
  res.lassomfdr[,i] <- locfdr.lam
  res.lassocv[,i] <- locfdr.cv
  res.uni[,i] <- univariate.res$fdr
}

save(trueA, trueA2, trueA3, trueB, trueB2, trueB3, false, false2, false3, bins, 
     uni.A, uni.B, uni.C, multi.A, multi.B, multi.C, multi.A2, multi.B2, multi.C2, 
     res.lassomfdr, res.lassocv, res.uni, file = "simres/logistic_easy.RData")


####################################################
##### Logistic reg "assumptions violated" scenario
####################################################


### Setup
n <- 200
p <- 600
bb <- numeric(60)
bb[(0:5)*10+1] <- c(1.25,-1.25,1,-1,.75,-.75)  #### 6 true variables (type A)

id.A <- which(bb != 0)
id.B <- 1:60
id.B <- id.B[-id.A]
id.C <- 61:p

nsim <- 200

trueA <- trueB <- false <- trueA2 <- trueB2 <- false2 <- matrix(NA, nrow = nsim, ncol = 5)
trueA3 <- trueB3 <- false3 <-  matrix(NA, nrow = nsim, ncol = 5)
fdr.true <- fdr.true.cv <- fdr.true.uni <- matrix(NA, nrow = nsim, ncol = 6)
multi.A <- multi.A2 <- multi.B <- multi.B2 <- multi.C <- multi.C2 <- rep(NA, nsim)
uni.A <- uni.B <- uni.C <- rep(NA,nsim)
res.lassomfdr <- res.lassocv <- res.uni <- matrix(NA, nrow = p, ncol = nsim)

for(i in 1:nsim){
  D1 <- genData(n, J=6, J0=6, K=10, K0=1, rho=0, rho.g=0.5, beta=bb, family = "binomial")  #### 9 correlated vars for each true (type B)
  D2 <- genData(n, p - 60, rho=0.8, beta=0, corr="auto")
  
  X <- cbind(D1$X, D2$X)
  y <- D1$y
  
  #### Univariate locfdr
  zstat <- tstat <- rep(NA, ncol(X))
  for (j in 1:ncol(X)){
    fit.lm <- glm(y ~ X[,j], family = "binomial")
    zstat[j] <- summary(fit.lm)$coefficients[2,3]
  }
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)
  
  
  ### Fit model
  X <- std(X)
  fit <- ncvreg(X,y, penalty = "lasso", family = "binomial", returnX = TRUE)
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv.fit <- cv.ncvreg(X,y,penalty = "lasso", family = "binomial")
  cv.lam <- which(fit$lambda == cv.fit$lambda.min)
  
  locfdr.lam <- locmfdr(fit, fit$lambda[mfdr.lam])[,2]
  locfdr.cv <- locmfdr(fit, cv.fit$lambda.min)[,2]
  
  #### store fdr for true vars
  fdr.true[i,] <- locfdr.lam[(0:5)*10+1]
  fdr.true.cv[i,] <- locfdr.cv[(0:5)*10+1]
  fdr.true.uni[i,] <- univariate.res$fdr[(0:5)*10+1]
  
  
  ### Seperate true/false vars into fdr bins 
  bins <- cut(locfdr.lam, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA[i,] <- table(bins[(0:5)*10+1])
  trueB[i,] <- table(bins[(1:60)[-((0:5)*10+1)]])
  false[i,] <- table(bins[61:p])
  
  bins2 <- cut(locfdr.cv, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA2[i,] <- table(bins2[(0:5)*10+1])
  trueB2[i,] <- table(bins2[(1:60)[-((0:5)*10+1)]])
  false2[i,] <- table(bins2[61:p])
  
  bins3 <- cut(univariate.res$fdr, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA3[i,] <- table(bins3[(0:5)*10+1])
  trueB3[i,] <- table(bins3[(1:60)[-((0:5)*10+1)]])
  false3[i,] <- table(bins3[61:p])
  
  
  ##### Power using 10% threshold
  uni.A[i] <- sum(univariate.res$fdr[id.A] < .1)
  uni.B[i] <- sum(univariate.res$fdr[id.B] < .1)
  uni.C[i] <- sum(univariate.res$fdr[id.C] < .1)
  
  multi.A[i] <- sum(locfdr.lam[id.A] < .1)
  multi.B[i] <- sum(locfdr.lam[id.B] < .1)
  multi.C[i] <- sum(locfdr.lam[id.C] < .1)
  
  multi.A2[i] <- sum(locfdr.cv[id.A] < .1)
  multi.B2[i] <- sum(locfdr.cv[id.B] < .1)
  multi.C2[i] <- sum(locfdr.cv[id.C] < .1)
  
  res.lassomfdr[,i] <- locfdr.lam
  res.lassocv[,i] <- locfdr.cv
  res.uni[,i] <- univariate.res$fdr
}

save(trueA, trueA2, trueA3, trueB, trueB2, trueB3, false, false2, false3, bins, 
     uni.A, uni.B, uni.C, multi.A, multi.B, multi.C, multi.A2, multi.B2, multi.C2,
     res.lassomfdr, res.lassocv, res.uni, file = "simres/logistic_hard.RData")


####################################################
##### Linear reg "assumptions met" scenario
####################################################

### Setup
n <- 1000
p <- 600

sig <- sqrt(n)
t <- 60
beta <- c(rnorm(t/2,4,1),-rnorm(t/2,4,1),rep(0,p-t))
nreps <- 200
pi0 = 1

trueA <- trueB <- false <- trueA2 <- trueB2 <- false2 <- matrix(NA, nrow = nreps, ncol = 5)
trueA3 <- trueB3 <- false3 <-  matrix(NA, nrow = nreps, ncol = 5)
res.lassomfdr <- res.lassocv <- res.uni <- matrix(NA, nrow = nreps, ncol = p)
multi.A <- multi.A2 <- multi.B <- multi.B2 <- multi.C <- multi.C2 <- rep(NA, nreps)
uni.A <- uni.B <- uni.C <- rep(NA,nreps)

for (i in 1:nreps){
  X <- matrix(rnorm(n*p), ncol = p, nrow = n)
  y <- X %*% beta + rnorm(n,sd = sig)
  
  ### Fit model
  X <- std(X)
  y <- scale(y, scale = FALSE)
  fit <- ncvreg(X,y, penalty = "lasso", returnX = TRUE)
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv.fit <- cv.ncvreg(X,y,penalty = "lasso")
  cv.lam <- which(fit$lambda == cv.fit$lambda.min)
  
  ### Calculate zj for each lambda (mfdr)
  R <- y - predict(fit, X, lambda = fit$lambda[mfdr.lam], type = "response")
  z <- (1/n)*t(X) %*% R + fit$beta[-1,mfdr.lam]
  ### Calculate zj for each lambda (CV)
  R.cv <- y - predict(fit, X, lambda = fit$lambda[cv.lam], type = "response")
  z.cv <- (1/n)*t(X) %*% R.cv + fit$beta[-1,cv.lam]
  
  
  ### Estimate locfdr for each lambda
  f <- density(z)
  ff <- approxfun(f$x, f$y)
  f.cv <- density(z.cv)
  ff.cv <- approxfun(f.cv$x, f.cv$y)
  S <- predict(fit, type = "nvars")
  
  sig.lam <- sqrt(fit$loss[mfdr.lam]/(n - S[mfdr.lam] + 1))
  sig.cv <- sqrt(fit$loss[cv.lam]/(n - S[cv.lam] + 1))
  locfdr.lam <- pmin(pi0*dnorm(z, 0, sig.lam/sqrt(n))/ff(z), 1)
  locfdr.cv <- pmin(pi0*dnorm(z.cv, 0, sig.cv/sqrt(n))/ff.cv(z.cv), 1)
  
  
  
  ### Estimate univariate locfdr
  tstat <- rep(NA, ncol(X))
  for (j in 1:ncol(X)){
    fit.lm <- lm(y ~ X[,j])
    tstat[j] <- summary(fit.lm)$coefficients[2,3]
  }
  zstat <- qnorm(pt(tstat,n - 2))
  zstat[is.infinite(zstat)] <- tstat[is.infinite(zstat)]  #### Make infinite z-stats into their pre-transformed t-stat
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)
  
  #### store fdr for true vars
  res.lassomfdr[i,] <- locfdr.lam
  res.lassocv[i,] <- locfdr.cv
  res.uni[i,] <- univariate.res$fdr
  
  
  ### Seperate true/false vars into fdr bins 
  bins <- cut(locfdr.lam, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA[i,] <- table(bins[1:t])
  trueB[i,] <- rep(0,5)
  false[i,] <- table(bins[(t+1):p])
  
  bins2 <- cut(locfdr.cv, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA2[i,] <- table(bins2[1:t])
  trueB2[i,] <- rep(0,5)
  false2[i,] <- table(bins2[(t+1):p])
  
  bins3 <- cut(univariate.res$fdr, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA3[i,] <- table(bins3[1:t])
  trueB3[i,] <- rep(0,5)
  false3[i,] <- table(bins3[(t+1):p])
  
  
  ##### Power using 10% threshold
  uni.A[i] <- sum(univariate.res$fdr[1:t] < .1)
  #uni.B[i] <- sum(univariate.res$fdr[id.B] < .1)
  uni.C[i] <- sum(univariate.res$fdr[(t+1):p] < .1)
  
  multi.A[i] <- sum(locfdr.lam[1:t] < .1)
  #multi.B[i] <- sum(locfdr.lam[id.B] < .1)
  multi.C[i] <- sum(locfdr.lam[(t+1):p] < .1)
  
  multi.A2[i] <- sum(locfdr.cv[1:t] < .1)
  #multi.B2[i] <- sum(locfdr.cv[id.B] < .1)
  multi.C2[i] <- sum(locfdr.cv[(t+1):p] < .1)
}

save(trueA, trueA2, trueA3, trueB, trueB2, trueB3, false, false2, false3, bins, 
     uni.A, uni.B, uni.C, multi.A, multi.B, multi.C, multi.A2, multi.B2, multi.C2, 
     res.lassomfdr, res.lassocv, res.uni, file = "simres/linear_easy.RData")

####################################################
##### Linear reg "assumptions violated" scenario
####################################################

### Setup
n <- 200
p <- 600
bb <- numeric(60)
bb[(0:5)*10+1] <- c(6,-6,5,-5,4,-4)

id.A <- which(bb != 0)
id.B <- 1:60
id.B <- id.B[-id.A]
id.C <- 61:p

sig <- sqrt(n)
nreps <- 200

trueA <- trueB <- false <- trueA2 <- trueB2 <- false2 <- matrix(NA, nrow = nreps, ncol = 5)
trueA3 <- trueB3 <- false3 <-  matrix(NA, nrow = nreps, ncol = 5)
res.lassomfdr <- res.lassocv <- res.uni <- matrix(NA, nrow = nreps, ncol = p)
tstat <- rep(NA, p)
multi.A <- multi.A2 <- multi.B <- multi.B2 <- multi.C <- multi.C2 <- rep(NA, nreps)
uni.A <- uni.B <- uni.C <- rep(NA,nreps)

for (i in 1:nreps){

  ### Corr data
  D1 <- genData(n, J=6, J0=6, K=10, K0=1, rho=0, rho.g=0.5, beta=bb)  #### 9 correlated vars for each true (type B)
  D2 <- genData(n, p - 60, rho=0.8, beta=0, corr="auto")

  X <- cbind(D1$X, D2$X)
  y <- D1$y
  
  
  ### Fit model
  X <- std(X)
  y <- scale(y, scale = FALSE)
  fit <- ncvreg(X,y, penalty = "lasso", returnX = TRUE)
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv.fit <- cv.ncvreg(X,y,penalty = "lasso")
  cv.lam <- which(fit$lambda == cv.fit$lambda.min)
  
  ### Calculate zj for each lambda (mfdr)
  R <- y - predict(fit, X, lambda = fit$lambda[mfdr.lam], type = "response")
  z <- (1/n)*t(X) %*% R + fit$beta[-1,mfdr.lam]
  ### Calculate zj for each lambda (CV)
  R.cv <- y - predict(fit, X, lambda = fit$lambda[cv.lam], type = "response")
  z.cv <- (1/n)*t(X) %*% R.cv + fit$beta[-1,cv.lam]
  
  
  ### Estimate locfdr for each lambda
  f <- density(z)
  ff <- approxfun(f$x, f$y)
  f.cv <- density(z.cv)
  ff.cv <- approxfun(f.cv$x, f.cv$y)
  S <- predict(fit, type = "nvars")
  
  sig.lam <- sqrt(fit$loss[mfdr.lam]/(n - S[mfdr.lam] + 1))
  sig.cv <- sqrt(fit$loss[cv.lam]/(n - S[cv.lam] + 1))
  locfdr.lam <- pmin(dnorm(z, 0, sig.lam/sqrt(n))/ff(z), 1)
  locfdr.cv <- pmin(dnorm(z.cv, 0, sig.cv/sqrt(n))/ff.cv(z.cv), 1)
  
  
  
  ### Estimate univariate locfdr
  tstat <- rep(NA, ncol(X))
  for (j in 1:ncol(X)){
    fit.lm <- lm(y ~ X[,j])
    tstat[j] <- summary(fit.lm)$coefficients[2,3]
  }
  zstat <- qnorm(pt(tstat,n - 2))
  zstat[is.infinite(zstat)] <- tstat[is.infinite(zstat)]  #### Make infinite z-stats into their pre-transformed t-stat
  univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)
  
  #### store fdr for true vars
  res.lassomfdr[i,] <- locfdr.lam
  res.lassocv[i,] <- locfdr.cv
  res.uni[i,] <- univariate.res$fdr
  
  
  ### Seperate true/false vars into fdr bins 
  bins <- cut(locfdr.lam, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA[i,] <- table(bins[(0:5)*10+1])
  trueB[i,] <- table(bins[(1:60)[-((0:5)*10+1)]])
  false[i,] <- table(bins[61:p])
  
  bins2 <- cut(locfdr.cv, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA2[i,] <- table(bins2[(0:5)*10+1])
  trueB2[i,] <- table(bins2[(1:60)[-((0:5)*10+1)]])
  false2[i,] <- table(bins2[61:p])
  
  bins3 <- cut(univariate.res$fdr, breaks  = c(0,.2,.4,.6,.8,1.1))
  trueA3[i,] <- table(bins3[(0:5)*10+1])
  trueB3[i,] <- table(bins3[(1:60)[-((0:5)*10+1)]])
  false3[i,] <- table(bins3[61:p])
  
  
  ##### Power using 10% threshold
  uni.A[i] <- sum(univariate.res$fdr[id.A] < .1)
  uni.B[i] <- sum(univariate.res$fdr[id.B] < .1)
  uni.C[i] <- sum(univariate.res$fdr[id.C] < .1)
  
  multi.A[i] <- sum(locfdr.lam[id.A] < .1)
  multi.B[i] <- sum(locfdr.lam[id.B] < .1)
  multi.C[i] <- sum(locfdr.lam[id.C] < .1)
  
  multi.A2[i] <- sum(locfdr.cv[id.A] < .1)
  multi.B2[i] <- sum(locfdr.cv[id.B] < .1)
  multi.C2[i] <- sum(locfdr.cv[id.C] < .1)
}

save(trueA, trueA2, trueA3, trueB, trueB2, trueB3, false, false2, false3, bins, 
     uni.A, uni.B, uni.C, multi.A, multi.B, multi.C, multi.A2, multi.B2, multi.C2, 
     res.lassomfdr, res.lassocv, res.uni, file = "simres/linear_hard.RData")

##########################################################################################
## At this point all six simulations have run and their results are stored as .RData files
## The code below uses these to generate Fig2 and Table1
##########################################################################################


########################################################################################
#
#           Code to produce Table 2 (requires simulation results from earlier)
#
########################################################################################

# First general lower row of fig2 (assumptions violated)
load("simres/linear_hard.RData")

### Below will treat 'B' variables as valid discoveries, comment out to exclude B
AC <- rep(1,600)
false.id <- c(rep(0,60),rep(1,540)) ### This counts 'B' variables as valid discoveries

df <- data.frame(est = as.vector(t(res.uni)), id = rep(false.id, nrow(res.uni)))
df2 <- data.frame(est = as.vector(t(res.lassocv)), id = rep(false.id, nrow(res.uni)))
df3 <- data.frame(est = as.vector(t(res.lassomfdr)), id = rep(false.id, nrow(res.uni)))

set.seed(999)
sample.id <- sample(1:nrow(df), size = 4000)  ### if too many points to plot (plotting them all sometimes crashes R)
#sample.id <- 1:nrow(df) ### if few enough sim reps to plot all pts

df.tot <- rbind(data.frame(df2[sample.id,], Method = "mfdr (CV)"),
                data.frame(df3[sample.id,], Method = "mfdr (mFdr)"),
                data.frame(df[sample.id,], Method = "Univariate fdr"))
line <- rbind(data.frame(Method = "Univariate fdr", a = lm.fit$coefficients[1], b = lm.fit$coefficients[2]),
              data.frame(Method = "mfdr (CV)", a = lm.fit2$coefficients[1], b = lm.fit2$coefficients[2]),
              data.frame(Method = "mfdr (mFdr)", a = lm.fit3$coefficients[1], b = lm.fit3$coefficients[2]))
ref <- data.frame(a = 0, b = 1)

p1 <- ggplot(df.tot, aes(est,id, col = Method)) + geom_jitter(width = 0, height = .03, alpha = .15) +
  xlab("Estimated local false discovery rate") + ylab("Noise Variable Indicator") + ggtitle("Linear - Violated Scenario") +
  geom_smooth(se = FALSE, method = "loess", span = 1.6) + 
  geom_abline(data = ref, aes(intercept = a,slope = b), linetype = 2, size = 1) +
  facet_grid(. ~ Method)


### Now generate upper row of Fig2 (assumptions met)
load("simres/linear_easy.RData")
false.id <- c(rep(0,60),rep(1,540))

df <- data.frame(est = as.vector(t(res.uni)), id = rep(false.id, nrow(res.uni)))
df2 <- data.frame(est = as.vector(t(res.lassocv)), id = rep(false.id, nrow(res.uni)))
df3 <- data.frame(est = as.vector(t(res.lassomfdr)), id = rep(false.id, nrow(res.uni)))

set.seed(999)
sample.id <- sample(1:nrow(df), size = 4000) ### if too many points to plot
#sample.id <- 1:nrow(df) 

df.tot <- rbind(data.frame(df2[sample.id,], Method = "mfdr (CV)"),
                data.frame(df3[sample.id,], Method = "mfdr (mFdr)"),
                data.frame(df[sample.id,], Method = "Univariate fdr"))

p2 <- ggplot(df.tot, aes(est,id, col = Method)) + geom_jitter(width = 0, height = .03, alpha = .15) +
  xlab("Estimated local false discovery rate") + ylab("Noise Variable Indicator") + ggtitle("Linear - Met Scenario") +
  geom_smooth(method = "loess", se = FALSE, span = 1.5) +
  geom_abline(data = ref, aes(intercept = a,slope = b), linetype = 2, size = 1) +
  facet_grid(. ~ Method)

### Generate and save Figure2 as a pdf
pdf("figures_tables/Fig2.pdf", height=5.5, width=8)
multiplot(p2, p1, cols=1)
dev.off()


########################################################################################
#
#           Code to produce Table 1 (requires simulation results from earlier)
#
########################################################################################

load(file = "simres/linear_hard.RData")
#### If only C are false discoveries
res.lin.mfdr <- colSums(false)/(colSums(trueA) + colSums(trueB) + colSums(false))
names(res.lin.mfdr) <- levels(bins)
res.lin.cv <- colSums(false2)/(colSums(trueA2) + colSums(trueB2) + colSums(false2))
names(res.lin.cv) <- levels(bins)
res.lin.uni <- colSums(false3)/(colSums(trueA3) + colSums(trueB3) + colSums(false3))
names(res.lin.uni) <- levels(bins)


load(file = "simres/logistic_hard.RData")
#### If only C are false discoveries
res.log.mfdr <- colSums(false)/(colSums(trueA) + colSums(trueB) + colSums(false))
names(res.log.mfdr) <- levels(bins)
res.log.cv <- colSums(false2)/(colSums(trueA2) + colSums(trueB2) + colSums(false2))
names(res.log.cv) <- levels(bins)
res.log.uni <- colSums(false3)/(colSums(trueA3) + colSums(trueB3) + colSums(false3))
names(res.log.uni) <- levels(bins)


load(file = "simres/cox_hard.RData")
#### If only C are false discoveries
res.cox.mfdr <- colSums(false)/(colSums(trueA) + colSums(trueB) + colSums(false))
names(res.cox.mfdr) <- levels(bins)
res.cox.cv <- colSums(false2)/(colSums(trueA2) + colSums(trueB2) + colSums(false2))
names(res.cox.cv) <- levels(bins)
res.cox.uni <- colSums(false3)/(colSums(trueA3) + colSums(trueB3) + colSums(false3))
names(res.cox.uni) <- levels(bins)

ttab1 <- round(rbind(res.lin.uni, res.lin.mfdr, res.lin.cv,
                     res.log.uni, res.log.mfdr, res.log.cv,
                     res.cox.uni, res.cox.mfdr, res.cox.cv), 2)

rownames(ttab1) <- c("Linear-Univariate", "Linear-mFDR-lam", "Linear-CV-lam",
                     "Logistic-Univariate", "Logistic-mFDR-lam", "Logistic-CV-lam",
                     "Cox-Univariate", "Cox-mFDR-lam", "Cox-CV-lam")

pdf("Table1.pdf", height=7, width=8.5)
grid.table(ttab1)
dev.off()


########################################################################################
#
#           Code to produce Figure 3 (requires simulation results from earlier)
#
########################################################################################

load("simres/linear_hard.RData")
results <- rbind(c(mean(uni.A), mean(multi.A), mean(multi.A2)),
                 c(mean(uni.B), mean(multi.B), mean(multi.B2)),
                 c(mean(uni.C), mean(multi.C), mean(multi.C2)))
colnames(results) <- c("Univariate fdr", "mfdr (mFdr)", "mfdr (CV)")
rownames(results) <- c("Causal (A)", "Correlated (B)", "Noise (C)")

#### Plot
dd <- melt(results)
dd1 <- data.frame(Var2 = dd$Var2, Feature = dd$Var1, value = dd$value, scenario = rep("Assumptions Violated", nrow(dd)))

load("simres/linear_easy.RData")
results <- rbind(c(mean(uni.A), mean(multi.A), mean(multi.A2)),
                 c(mean(uni.B), mean(multi.B), mean(multi.B2)),
                 c(mean(uni.C), mean(multi.C), mean(multi.C2)))
colnames(results) <- c("Univariate fdr", "mfdr (mFdr)", "mfdr (CV)")
rownames(results) <- c("Causal (A)", "Correlated (B)", "Noise (C)")

dd <- melt(results)
dd2 <- data.frame(Var2 = dd$Var2, Feature = dd$Var1, value = dd$value, scenario = rep("Assumptions Met", nrow(dd)))
ddd <- rbind(dd1, dd2)

pdf("Fig3.pdf", height=5, width=7.5)
p <- ggplot(ddd, aes(x = Var2, y = value, fill = Feature)) + ggtitle("") +
  geom_bar(stat = "identity", position=position_stack(reverse=TRUE)) + scale_fill_manual(name = "Feature Type", values = pal(3)[c(2,3,1)]) +
  coord_flip() + scale_y_continuous(expression(paste("Features with ", widehat(mfdr), " or ", widehat(fdr), " < 0.10"))) + scale_x_discrete("Method")
p + facet_grid(. ~ scenario, scales = "free")
dev.off()

########################################################################################
#
#           Simulations discussed in Sec 4.3 (comparisons to Selective Inf and SS)
#
########################################################################################

### Assumptions met
n <- 1000
p <- 600

sig <- sqrt(n)
t <- 60
beta <- c(rep(4,t/2),rep(-4,t/2),rep(0,p-t))
nreps <- 200
pi0 = 1

fdr.true <- fdr.true.cv <- fdr.true.uni <- matrix(NA, nrow = nreps, ncol = p)
tres.pv <- fres.pv <- tres.spacing <- fres.spacing <- tres.modspacing <- tres.knockoff <- NULL
fres.modspacing <- tres.covtest <- fres.covtest <- tres.mfdr <- fres.mfdr <- fres.knockoff <- NULL
tres.ss <- cres.ss <- frate.ss <- numeric(nreps)
true.cv <- false.cv <- numeric(nreps)

for (i in 1:nreps){
  X <- matrix(rnorm(n*p), ncol = p, nrow = n)
  y <- X %*% beta + rnorm(n,sd = sig)
  
  ### Fit model
  X <- std(X)
  y <- scale(y, scale = FALSE)
  fit <- ncvreg(X,y, penalty = "lasso", returnX = TRUE)
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv.fit <- cv.ncvreg(X,y,penalty = "lasso")
  cv.lam <- which(fit$lambda == cv.fit$lambda.min)
  
  true.cv[i] <- sum(cv.fit$fit$beta[2:61,cv.lam] != 0)
  false.cv[i] <- sum(cv.fit$fit$beta[62:601,cv.lam] != 0)
  
  
  ### Calculate zj for each lambda (mfdr)
  R <- y - predict(fit, X, lambda = fit$lambda[mfdr.lam], type = "response")
  z <- (1/n)*t(X) %*% R + fit$beta[-1,mfdr.lam]
  ### Calculate zj for each lambda (CV)
  R.cv <- y - predict(fit, X, lambda = fit$lambda[cv.lam], type = "response")
  z.cv <- (1/n)*t(X) %*% R.cv + fit$beta[-1,cv.lam]
  
  
  ### Estimate locfdr for each lambda
  f <- density(z)
  ff <- approxfun(f$x, f$y)
  f.cv <- density(z.cv)
  ff.cv <- approxfun(f.cv$x, f.cv$y)
  S <- predict(fit, type = "nvars")
  
  sig.lam <- sqrt(fit$loss[mfdr.lam]/(n - S[mfdr.lam] + 1))
  sig.cv <- sqrt(fit$loss[cv.lam]/(n - S[cv.lam] + 1))
  locfdr.lam <- pmin(pi0*dnorm(z, 0, sig.lam/sqrt(n))/ff(z), 1)
  locfdr.cv <- pmin(pi0*dnorm(z.cv, 0, sig.cv/sqrt(n))/ff.cv(z.cv), 1)
  
  
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
  
  #### store fdr for all vars
  fdr.true[i,] <- locfdr.lam
  fdr.true.cv[i,] <- locfdr.cv
  
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

save(final.true, final.false, t, fdr.true, fdr.true.cv,  true.cv, false.cv, p, 
     file = "simres/model_easy.RData")


### Assumptions Violated
### Setup
n <- 200
p <- 600
bb <- numeric(60)
bb[(0:5)*10+1] <- c(5.5,-5.5,5,-5,4.5,-4.5)

id.A <- which(bb != 0)
id.B <- 1:60
id.B <- id.B[-id.A]
id.C <- 61:p

t <- 60
nreps <- 200
pi0 = 1

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
  
  
  ### Fit model
  X <- std(X)
  y <- scale(y, scale = FALSE)
  fit <- ncvreg(X,y, penalty = "lasso", returnX = TRUE)
  mfdr.lam <- max(which(mfdr(fit)$mFDR < .1))
  cv.fit <- cv.ncvreg(X,y,penalty = "lasso")
  cv.lam <- which(fit$lambda == cv.fit$lambda.min)
  
  true.cv[i] <- sum(cv.fit$fit$beta[id.A + 1,cv.lam] != 0)
  cor.cv[i] <- sum(cv.fit$fit$beta[id.B + 1,cv.lam] != 0)
  false.cv[i] <- sum(cv.fit$fit$beta[id.C + 1,cv.lam] != 0)
  
  ### Calculate zj for each lambda (mfdr)
  R <- y - predict(fit, X, lambda = fit$lambda[mfdr.lam], type = "response")
  z <- (1/n)*t(X) %*% R + fit$beta[-1,mfdr.lam]
  ### Calculate zj for each lambda (CV)
  R.cv <- y - predict(fit, X, lambda = fit$lambda[cv.lam], type = "response")
  z.cv <- (1/n)*t(X) %*% R.cv + fit$beta[-1,cv.lam]
  
  
  ### Estimate locfdr for each lambda
  f <- density(z)
  ff <- approxfun(f$x, f$y)
  f.cv <- density(z.cv)
  ff.cv <- approxfun(f.cv$x, f.cv$y)
  S <- predict(fit, type = "nvars")
  
  sig.lam <- sqrt(fit$loss[mfdr.lam]/(n - S[mfdr.lam] + 1))
  sig.cv <- sqrt(fit$loss[cv.lam]/(n - S[cv.lam] + 1))
  locfdr.lam <- pmin(pi0*dnorm(z, 0, sig.lam/sqrt(n))/ff(z), 1)
  locfdr.cv <- pmin(pi0*dnorm(z.cv, 0, sig.cv/sqrt(n))/ff.cv(z.cv), 1)
  
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
  
  #### store fdr for all vars
  fdr.true[i,] <- locfdr.lam
  fdr.true.cv[i,] <- locfdr.cv
  
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

save(id.A, id.B, id.C, final.true, final.cor, final.false, fdr.true, fdr.true.cv,
     true.cv, cor.cv, false.cv, p, file = "simres/model_hard.RData")

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
selinf.C <- colSums(final.false)/(colSums(final.true) + colSums(final.false))

locmfdr.A <- sum(colMeans(fdr.true[,1:t] < .10))
locmfdr.nC <- sum(colMeans(fdr.true[,(t+1):p] < .10))
locmfdr.C <- sum(colSums(fdr.true[,(t+1):p] < .10))/sum(colSums(fdr.true < .10))

RR <- data.frame(rbind(c(locmfdr.A, selinf.A),
                       NA,
                       c(locmfdr.nC, selinf.nC),
                       c(locmfdr.C, selinf.C)))

rownames(RR) <- c("Avg 'A'", "Avg 'B'", "Avg C", "Avg FDR")
colnames(RR) <- c("loc-mfdr", colnames(RR)[-1])


load(file = "simres/model_hard.RData")

selinf.A <- apply(final.true, 2, mean)
selinf.B <- apply(final.cor, 2, mean)
final.false[is.na(final.false)] <- 0 ### Fdr of 0/0 is defined as 0
selinf.nC <- apply(final.false, 2, mean)
selinf.C <- colSums(final.false)/(colSums(final.true) + colSums(final.false))

locmfdr.A <- sum(colMeans(fdr.true[,id.A] < .10))
locmfdr.B <- sum(colMeans(fdr.true[,id.B] < .10))
locmfdr.nC <- sum(colMeans(fdr.true[,id.C] < .10))
locmfdr.C <- sum(colSums(fdr.true[,id.C] < .10))/sum(colSums(fdr.true < .10))

RR2 <- data.frame(rbind(c(locmfdr.A,selinf.A),
                        c(locmfdr.B,selinf.B),
                        c(locmfdr.nC,selinf.nC),
                        c(locmfdr.C,selinf.C)))

rownames(RR2) <- c("Avg 'A'", "Avg 'B'", "Avg C", "Avg FDR")
colnames(RR2) <- c("loc-mfdr", colnames(RR2)[-1])

pdf("Table2.pdf", height=7, width=8.5)
grid.table(round(t(rbind(RR2, RR)), 2))
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
XX <- std(X)
ZZ <- model.matrix(~ factor(Sex) + factor(Race) + factor(AdjChemo) + factor(SmHist) + factor(Margin) + factor(Grade), Z)
w <- rep(0:1, c(ncol(ZZ)-1, ncol(X))) ## Don't the penalize clinical covariates
fit2 <- ncvsurv(cbind(std(ZZ[,-1]), XX), S, penalty = 'lasso', returnX = TRUE, penalty.factor = w)

### local mfdr at MFdr lambda
step <- max(which(mfdr(fit2)$mFDR < .1))
shedden.locfdr <- case.locmfdr(fit2, fit2$lambda[step])

id <- order(shedden.locfdr[,2])[1:10] - 17  ## Subtract 17 because of 17 unpenalized vars
mfdr.output <- data.frame(name = colnames(X[,id]),
                          z = shedden.locfdr[order(shedden.locfdr[,2]),][1:10,1],
                          fdr = shedden.locfdr[order(shedden.locfdr[,2]),][1:10,2])


## local mfdr at CV lambda
set.seed(12345)
cv.fit2 <- cv.ncvsurv(cbind(std(ZZ[,-1]), XX), S, penalty = 'lasso', penalty.factor = w)
cv.lambda <- cv.fit2$lambda.min

shedden.locfdr2 <- case.locmfdr(fit2, cv.lambda)

idd <- order(shedden.locfdr2[,2])[1:10] - 17  ## Subtract 17 because of 17 unpenalized vars
cv.output <- data.frame(name = colnames(X[,idd]),
                        z = shedden.locfdr2[order(shedden.locfdr2[,2]),][1:10,1],
                        fdr = shedden.locfdr2[order(shedden.locfdr2[,2]),][1:10,2])


#### Univariate locfdr (for comparison)
zstat <- rep(NA, ncol(XX))
for (j in 1:ncol(XX)){
  fit.uni <- coxph(S ~ XX[,j] + ZZ[,-1])
  zstat[j] <- summary(fit.uni)$coefficients[1,4]
}
univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)

uni.output <- data.frame(name = colnames(X[,order(univariate.res$fdr)[1:10]]),
                         z = zstat[order(univariate.res$fdr)][1:10],
                         fdr = univariate.res$fdr[order(univariate.res$fdr)][1:10])

#############
### Table 3
#############

u.names <- fData[row.names(fData) %in% uni.output$name,]
u.names2 <- u.names[match(uni.output$name, rownames(u.names)),1]

mf.names <- fData[row.names(fData) %in% mfdr.output$name,]
mf.names2 <- mf.names[match(mfdr.output$name, rownames(mf.names)),1]

cv.names <- fData[row.names(fData) %in% cv.output$name,]
cv.names2 <- cv.names[match(cv.output$name, rownames(cv.names)),1]

## note: the <NA> names are missing in the fData file, they are replaced 
## with the probe_id in the table appearing in the manuscript

pdf("Table3.pdf", height=7, width=8.5)
grid.table((data.frame(feature.uni = u.names2, fdr.uni = uni.output$fdr,
                       feature.mfdr = mf.names2, mfdr.mfdr = mfdr.output$fdr,
                       feature.cv = cv.names2, mfdr.cv = cv.output$fdr)))
dev.off()



#############
### Figure 4
#############

pdf("Fig4.pdf", height=4, width=5.5)
plot(density(shedden.locfdr2$z), lwd = 2, col = 2, xlab = "z", main = "Density Comparisons")
lines(seq(-5,5,by = .01), dnorm(seq(-5,5,by = .01)), col = 1, lwd = 2, lty = 2)
lines(density(shedden.locfdr$z), lwd = 2, col = 3)
lines(density(zstat), lwd = 2, col = 4)
legend("topright", legend = c("lasso (CV)", "lasso (mfdr)", "Univariate", "Theoretical Null"), col = c(2,3,4,1), lty = c(1,1,1,2), lwd = 2)
dev.off()




### Case study #2 (table 4)

## Load TCGA BRCA1 dataset
TCGA <- readRDS("case_study_data//bcTCGA.rds")

## Standardize and fit penalized regression models
X <- TCGA$X
XX <- std(X)
y <- TCGA$y
fit <- ncvreg(XX, y, penalty = "lasso", returnX = TRUE)
cv.fit <- cv.ncvreg(XX, y, penalty = "lasso")
step <- max(which(mfdr(fit)[,3] < .1))

### local mfdr at two different lambda values
lassofdr.res1 <- locmfdr(fit, fit$lambda[step])
lassofdr.res2 <- locmfdr(fit, cv.fit$lambda.min)

### Univariate testing
tstat <- numeric(ncol(X))
for (j in 1:ncol(X)){
  fit.lm <- lm(y ~ X[,j])
  tstat[j] <- summary(fit.lm)$coefficients[2,3]
}
zstat <- qnorm(pt(tstat,n - 2))
zstat[is.infinite(zstat)] <- tstat[is.infinite(zstat)] ## Make infinite z-stats into their pre-transformed t-stat

univariate.res <- locfdr(zstat, nulltype = 0, plot = 0)

id <- order(univariate.res$fdr)[1:10]  ## Top 10 univariate selections
uni.output <- data.frame(name = colnames(X[,order(univariate.res$fdr)[1:10]]),
                         z = zstat[order(univariate.res$fdr)][1:10],
                         fdr = univariate.res$fdr[order(univariate.res$fdr)][1:10])


#############
### Table 4
#############

# Note: Chromosome locations are not contained in these data and were manually added to the table appearing in the manuscript

pdf("Table4.pdf", height=7, width=8.5)
grid.table((data.frame(Gene = uni.output$name, uni_loc = uni.output$fdr,
                       mfdr_locmfdr = lassofdr.res1[id,2],
                       cv_locmfdr= lassofdr.res1[id,2])))
dev.off()

