library(hdi)
library(ncvreg)
library(ashr)
library(Matrix)

## Sources functions used in the simulations and the construction of figures
source("C:\\Users\\millerry\\Documents\\loc-mfdr-paper\\functions.R")

### Setup parameters
n <- 200
t <- 60
bb <- numeric(t)
bb[(0:5)*10+1] <- c(6,-6,5,-5,4,-4)  #### 6 true variables (type A)

## Ids of var types
id.A <- which(bb != 0)
id.B <- 1:t
id.B <- id.B[-id.A]

## For n = 200, vary p (100, 200, 400, 800), rho (0, .3, .6, .9), corr (exh, auto)
## Store mfdr's and proj_fdr's for each across 100 iterations

nreps <- 100
rhos <- c(0,.3,.6,.9)
cors <-  c("auto", "exch")


ps <- c(100,200,400,800,1600)


for(q in 1:length(ps)) {
p <- ps[q]
id.C <- (t+1):p

proj_res <- array(NA, dim = c(nreps, p, length(rhos), length(cors)))
mfdr_res <- array(NA, dim = c(nreps, p, length(rhos), length(cors)))


for(k in 1:length(cors)){
  for(j in 1:length(rhos)){
      for(i in 1:nreps){
### Generate data using funs sourced
D1 <- genData(n, J=6, J0=6, K=10, K0=1, rho=0, rho.g=0.5, beta=bb)  #### 9 correlated vars for each true (type B)
D2 <- genData(n, p - 60, rho= rhos[j], beta=0, corr = cors[k])
X <- cbind(D1$X, D2$X)
X <- ncvreg::std(X)
y <- D1$y

lpj <- lasso.proj(x = X, y = y,standardize = FALSE, return.Z = TRUE)
zpj <- lpj$bhat/lpj$se
pfdr <- get_lfdr(ash(as.vector(zpj),rep(1, length(zpj))))

fit <- cv.ncvreg(X, y)
res <- local_mfdr(fit$fit, lambda = fit$lambda.min)
mfdr <- get_lfdr(ash(betahat = res$z, se = rep(1, length(res$z))))

## dim1 = rep, dim2 = p, dim3 = rho, dim4 = corr
proj_res[i, ,j,k] <- pfdr
mfdr_res[i, ,j,k] <- mfdr
      
}}}
save(proj_res, mfdr_res, p, n, t, bb, id.A, id.B, id.C, rhos, cors,
     file = paste0("H:\\sim1p",p,".RData"))
}


#############################
#### rho = 0, cor = auto
#############################
dfp <- data.frame(Estimate = c(as.vector(proj_res[,1:t,1,1]), 
                               as.vector(proj_res[,(t+1):nrow(proj_res),1,1])),
                  Indicator = c(rep(0, t*nreps),
                                rep(1, (p-t)*nreps)))

dfm <- data.frame(Estimate = c(as.vector(mfdr_res[,1:t,1,1]),
                               as.vector(mfdr_res[,(t+1):nrow(mfdr_res),1,1])),
                  Indicator = c(rep(0, t*nreps),
                                rep(1, (p-t)*nreps)))
                                                 
### Loess curves
xx <- seq(0, 1, by = .01)
                     
sm.dfp <- with(dfp, loess(Indicator ~ Estimate, span = .15))
y.dfp <- predict(sm.dfp, newdata = xx)
smooth.dfp <- data.frame(xx = xx, yy = y.dfp)

p.dfp1 <- ggplot(dfp) + geom_jitter(aes(x = Estimate, y = Indicator), height = .1, alpha = .1)
p.dfp1 <- p.dfp1 + geom_path(data = smooth.dfp, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "Proportion (noise)")
p.dfp1

sm.dfm <- with(dfm, loess(Indicator ~ Estimate, span = .15))
y.dfm <- predict(sm.dfm, newdata = xx)
smooth.dfm <- data.frame(xx = xx, yy = y.dfm)

p.dfm1 <- ggplot(dfm) + geom_jitter(aes(x = Estimate, y = Indicator), height = .1, alpha = .1)
p.dfm1 <- p.dfm1 + geom_path(data = smooth.dfm, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "Proportion (noise)")
p.dfm1


#############################
#### rho = 0.5, cor = auto
#############################
dfp <- data.frame(Estimate = c(as.vector(proj_res[,1:t,2,1]), 
                               as.vector(proj_res[,(t+1):nrow(proj_res),2,1])),
                  Indicator = c(rep(0, t*nreps),
                                rep(1, (p-t)*nreps)))

dfm <- data.frame(Estimate = c(as.vector(mfdr_res[,1:t,2,1]),
                               as.vector(mfdr_res[,(t+1):nrow(mfdr_res),2,1])),
                  Indicator = c(rep(0, t*nreps),
                                rep(1, (p-t)*nreps)))

### Loess curves
xx <- seq(0, 1, by = .01)

sm.dfp <- with(dfp, loess(Indicator ~ Estimate, span = .15))
y.dfp <- predict(sm.dfp, newdata = xx)
smooth.dfp <- data.frame(xx = xx, yy = y.dfp)

p.dfp2 <- ggplot(dfp) + geom_jitter(aes(x = Estimate, y = Indicator), height = .1, alpha = .1)
p.dfp2 <- p.dfp2 + geom_path(data = smooth.dfp, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "Proportion (noise)")
p.dfp2

sm.dfm <- with(dfm, loess(Indicator ~ Estimate, span = .15))
y.dfm <- predict(sm.dfm, newdata = xx)
smooth.dfm <- data.frame(xx = xx, yy = y.dfm)

p.dfm2 <- ggplot(dfm) + geom_jitter(aes(x = Estimate, y = Indicator), height = .1, alpha = .1)
p.dfm2 <- p.dfm2 + geom_path(data = smooth.dfm, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "Proportion (noise)")
p.dfm2


#############################
#### rho = 0, cor = exch
#############################
dfp <- data.frame(Estimate = c(as.vector(proj_res[,1:t,1,2]), 
                               as.vector(proj_res[,(t+1):nrow(proj_res),1,2])),
                  Indicator = c(rep(0, t*nreps),
                                rep(1, (p-t)*nreps)))

dfm <- data.frame(Estimate = c(as.vector(mfdr_res[,1:t,1,2]),
                               as.vector(mfdr_res[,(t+1):nrow(mfdr_res),1,2])),
                  Indicator = c(rep(0, t*nreps),
                                rep(1, (p-t)*nreps)))

### Loess curves
xx <- seq(0, 1, by = .01)

sm.dfp <- with(dfp, loess(Indicator ~ Estimate, span = .15))
y.dfp <- predict(sm.dfp, newdata = xx)
smooth.dfp <- data.frame(xx = xx, yy = y.dfp)

p.dfp3 <- ggplot(dfp) + geom_jitter(aes(x = Estimate, y = Indicator), height = .1, alpha = .1)
p.dfp3 <- p.dfp3 + geom_path(data = smooth.dfp, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "Proportion (noise)")
p.dfp3

sm.dfm <- with(dfm, loess(Indicator ~ Estimate, span = .15))
y.dfm <- predict(sm.dfm, newdata = xx)
smooth.dfm <- data.frame(xx = xx, yy = y.dfm)

p.dfm3 <- ggplot(dfm) + geom_jitter(aes(x = Estimate, y = Indicator), height = .1, alpha = .1)
p.dfm3 <- p.dfm3 + geom_path(data = smooth.dfm, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "Proportion (noise)")
p.dfm3


#############################
#### rho = 0.5, cor = exch
#############################
dfp <- data.frame(Estimate = c(as.vector(proj_res[,1:t,2,2]), 
                               as.vector(proj_res[,(t+1):nrow(proj_res),2,2])),
                  Indicator = c(rep(0, t*nreps),
                                rep(1, (p-t)*nreps)))

dfm <- data.frame(Estimate = c(as.vector(mfdr_res[,1:t,2,2]),
                               as.vector(mfdr_res[,(t+1):nrow(mfdr_res),2,2])),
                  Indicator = c(rep(0, t*nreps),
                                rep(1, (p-t)*nreps)))

### Loess curves
xx <- seq(0, 1, by = .01)

sm.dfp <- with(dfp, loess(Indicator ~ Estimate, span = .15))
y.dfp <- predict(sm.dfp, newdata = xx)
smooth.dfp <- data.frame(xx = xx, yy = y.dfp)

p.dfp4 <- ggplot(dfp) + geom_jitter(aes(x = Estimate, y = Indicator), height = .1, alpha = .1)
p.dfp4 <- p.dfp4 + geom_path(data = smooth.dfp, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "Proportion (noise)")
p.dfp4

sm.dfm <- with(dfm, loess(Indicator ~ Estimate, span = .15))
y.dfm <- predict(sm.dfm, newdata = xx)
smooth.dfm <- data.frame(xx = xx, yy = y.dfm)

p.dfm4 <- ggplot(dfm) + geom_jitter(aes(x = Estimate, y = Indicator), height = .1, alpha = .1)
p.dfm4 <- p.dfm4 + geom_path(data = smooth.dfm, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "Proportion (noise)")
p.dfm4

library(gridExtra)
grid.arrange(arrangeGrob(p.dfp1, top="rho = 0"),  
             arrangeGrob(p.dfp2, top="rho = 0.5", right = "De-biased z-statistics"),
             arrangeGrob(p.dfm1), 
             arrangeGrob(p.dfm2, right = "mfdr z-statistics"), ncol=2)

grid.arrange(arrangeGrob(p.dfp3, top="rho = 0"),  
             arrangeGrob(p.dfp4, top="rho = 0.5", right = "De-biased z-statistics"),
             arrangeGrob(p.dfm3), 
             arrangeGrob(p.dfm4, right = "mfdr z-statistics"), ncol=2)






#############################
#### Power
#############################
thr <- 0.1

# Exch, rgho = 0
PA3 <- sum(proj_res[,id.A,1,2] < thr)/nreps 
PB3 <- sum(proj_res[,id.B,1,2] < thr)/nreps 
PC3 <- sum(proj_res[,id.C,1,2] < thr)/nreps 

MA3 <- sum(mfdr_res[,id.A,1,2] < thr)/nreps 
MB3 <- sum(mfdr_res[,id.B,1,2] < thr)/nreps 
MC3 <- sum(mfdr_res[,id.C,1,2] < thr)/nreps 

# Exch, rho = .5
PA4 <- sum(proj_res[,id.A,2,2] < thr)/nreps 
PB4 <- sum(proj_res[,id.B,2,2] < thr)/nreps 
PC4 <- sum(proj_res[,id.C,2,2] < thr)/nreps 

MA4 <- sum(mfdr_res[,id.A,2,2] < thr)/nreps 
MB4 <- sum(mfdr_res[,id.B,2,2] < thr)/nreps 
MC4 <- sum(mfdr_res[,id.C,2,2] < thr)/nreps 

pp <- data.frame(Val = c(PA3, PB3, PC3, MA3, MB3, MC3, PA4, PB4, PC4, MA4, MB4, MC4),
                 Type = c(rep(c("A","B","C"), 4)),
                 Stat = c("unbiased", "unbiased", "unbiased", "mfdr", "mfdr", "mfdr", "unbiased", "unbiased", "unbiased", "mfdr", "mfdr", "mfdr"),
                 Rho = c(rep("Rho = 0", 6), rep("Rho = 0.5", 6)))

pp$Type <- fct_reorder(pp$Type, desc(pp$Type))

# Create the barplot
library(forcats)
ggplot(data=pp, aes(x=Rho, y=Val, fill=Type)) +
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Paired")+
  theme_minimal() + facet_wrap(~Stat)

  