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

### Install latest version of ncvreg
#devtools::install_github("pbreheny/ncvreg")
library(ncvreg)

## Sets the working directory to current location of this folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Sources functions used in the simulations and the construction of figures
source("functions.R")

#########################################################
###
###  Supplemental Figure #1 (Logistic Reg calibration)
###
#########################################################

# First general lower row of fig2 (assumptions violated)
load("simres/logistic_hard.RData")

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


l.cv <- with(df.cv, loess(Noise ~ Estimate, span = .02))
yy.cv <- predict(l.cv, newdata = xx)
smooth.cv <- data.frame(xx = xx, yy = yy.cv)

l.cv1se <- with(df.cv1se, loess(Noise ~ Estimate, span = .15))
yy.cv1se <- predict(l.cv1se, newdata = xx)
smooth.cv1se <- data.frame(xx = xx, yy = yy.cv1se)

l.mfdr <- with(df.mfdr, loess(Noise ~ Estimate, span = .15))
yy.mfdr <- predict(l.mfdr, newdata = xx)
smooth.mfdr <- data.frame(xx = xx, yy = yy.mfdr)

l.uni <- with(df.uni, loess(Noise ~ Estimate, span = .5))
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
p.cv <- p.cv + geom_path(data = smooth.cv, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda CV) - Violated")

p.cv1se <- p.cv1se + geom_path(data = smooth.cv1se, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda CV-1se) - Violated")

p.mfdr <- p.mfdr + geom_path(data = smooth.mfdr, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda 10% mFDR) - Violated")

p.uni <- p.uni + geom_path(data = smooth.uni, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "fdr (locfdr) - Violated")

p.uniash <- p.uniash + geom_path(data = smooth.uniash, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "fdr (ashr) - Violated")


### Now generate upper row of Fig2 (assumptions met)
load("simres/logistic_easy.RData")

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

l.uni <- with(df.uni, loess(Noise ~ Estimate, span = .6))
yy.uni <- predict(l.uni, newdata = xx)
smooth.uni <- data.frame(xx = xx, yy = yy.uni)

l.uniash <- with(df.uniash, loess(Noise ~ Estimate, span = .15))
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
p.cv2 <- p.cv2 + geom_path(data = smooth.cv, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda CV) - Met")

p.cv1se2 <- p.cv1se2 + geom_path(data = smooth.cv1se, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda CV-1se) - Met")

p.mfdr2 <- p.mfdr2 + geom_path(data = smooth.mfdr, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda 10% mFDR) - Met")

p.uni2 <- p.uni2 + geom_path(data = smooth.uni, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "fdr (locfdr) - Met")

p.uniash2 <- p.uniash2 + geom_path(data = smooth.uniash, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "fdr (ashr) - Met")


### Supplemental figure 1 as a pdf
pdf("supplemental_figures_tables/Supp1.pdf", height=7, width=10)
grid.arrange(p.cv, p.cv1se, p.mfdr, p.uni, p.uniash, 
             p.cv2, p.cv1se2, p.mfdr2, p.uni2, p.uniash2, nrow = 2)
dev.off()

#########################################################
###
###  Supplemental Figure #2 (Cox Reg calibration)
###
#########################################################


# First general lower row of fig2 (assumptions violated)
load("simres/cox_hard.RData")

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

l.mfdr <- with(df.mfdr, loess(Noise ~ Estimate, span = .25))
yy.mfdr <- predict(l.mfdr, newdata = xx)
smooth.mfdr <- data.frame(xx = xx, yy = yy.mfdr)

l.uni <- with(df.uni, loess(Noise ~ Estimate, span = .5))
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
p.cv <- p.cv + geom_path(data = smooth.cv, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda CV) - Violated")

p.cv1se <- p.cv1se + geom_path(data = smooth.cv1se, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda CV-1se) - Violated")

p.mfdr <- p.mfdr + geom_path(data = smooth.mfdr, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda 10% mFDR) - Violated")

p.uni <- p.uni + geom_path(data = smooth.uni, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "fdr (locfdr) - Violated")

p.uniash <- p.uniash + geom_path(data = smooth.uniash, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "fdr (ashr) - Violated")


### Now generate upper row of Fig2 (assumptions met)
load("simres/cox_easy.RData")

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

l.uni <- with(df.uni, loess(Noise ~ Estimate, span = .6))
yy.uni <- predict(l.uni, newdata = xx)
smooth.uni <- data.frame(xx = xx, yy = yy.uni)

l.uniash <- with(df.uniash, loess(Noise ~ Estimate, span = .15))
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
p.cv2 <- p.cv2 + geom_path(data = smooth.cv, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda CV) - Met")

p.cv1se2 <- p.cv1se2 + geom_path(data = smooth.cv1se, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda CV-1se) - Met")

p.mfdr2 <- p.mfdr2 + geom_path(data = smooth.mfdr, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mfdr (lambda 10% mFDR) - Met")

p.uni2 <- p.uni2 + geom_path(data = smooth.uni, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "fdr (locfdr) - Met")

p.uniash2 <- p.uniash2 + geom_path(data = smooth.uniash, aes(x = xx, y = yy)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "fdr (ashr) - Met")

### Save figure
pdf("supplemental_figures_tables/Supp2.pdf", height=7, width=10)
grid.arrange(p.cv, p.cv1se, p.mfdr, p.uni, p.uniash, 
             p.cv2, p.cv1se2, p.mfdr2, p.uni2, p.uniash2, nrow = 2)
dev.off()

