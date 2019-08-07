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



##############################################################
###
###  Supplemental Figure #1 (Logistic Reg calibration, ashr)
###
##############################################################

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

l.cv <- with(df.cv, loess(Noise ~ Estimate, span = .05))
yy.cv <- predict(l.cv, newdata = xx)
smooth.cv <- data.frame(xx = xx, yy = yy.cv)

l.cv1se <- with(df.cv1se, loess(Noise ~ Estimate, span = .05))
yy.cv1se <- predict(l.cv1se, newdata = xx)
smooth.cv1se <- data.frame(xx = xx, yy = yy.cv1se)

l.mfdr <- with(df.mfdr, loess(Noise ~ Estimate, span = .05))
yy.mfdr <- predict(l.mfdr, newdata = xx)
smooth.mfdr <- data.frame(xx = xx, yy = yy.mfdr)

l.uni <- with(df.uni, loess(Noise ~ Estimate, span = .6))
yy.uni <- predict(l.uni, newdata = xx)
smooth.uni <- data.frame(xx = xx, yy = yy.uni)

l.uniash <- with(df.uniash, loess(Noise ~ Estimate, span = .1))
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


### Now generate upper row
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

l.cv <- with(df.cv, loess(Noise ~ Estimate, span = .2))
yy.cv <- predict(l.cv, newdata = xx)
smooth.cv <- data.frame(xx = xx, yy = yy.cv)

l.cv1se <- with(df.cv1se, loess(Noise ~ Estimate, span = .15))
yy.cv1se <- predict(l.cv1se, newdata = xx)
smooth.cv1se <- data.frame(xx = xx, yy = yy.cv1se)

l.mfdr <- with(df.mfdr, loess(Noise ~ Estimate, span = .2))
yy.mfdr <- predict(l.mfdr, newdata = xx)
smooth.mfdr <- data.frame(xx = xx, yy = yy.mfdr)

l.uni <- with(df.uni, loess(Noise ~ Estimate, span = .65))
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


### Supplemental figure 1 as a pdf
png("supplemental_figures_tables/Supp1.png", h=6, w=8.5, units = 'in', res = 300)
grid.arrange(arrangeGrob(p.cv, top="lasso mfdr (CV)"), arrangeGrob(p.cv1se, top="lasso mfdr (CV)"),  
             arrangeGrob(p.uniash, top="Univariate fdr", right = "Assumptions Violated"),
             arrangeGrob(p.cv2), arrangeGrob(p.cv1se2),  
             arrangeGrob(p.uniash2, right = "Assumptions Met"), ncol=3)
dev.off()

#########################################################
###
###  Supplemental Figure #2 (Cox Reg calibration)
###
#########################################################

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

l.cv <- with(df.cv, loess(Noise ~ Estimate, span = .05))
yy.cv <- predict(l.cv, newdata = xx)
smooth.cv <- data.frame(xx = xx, yy = yy.cv)

l.cv1se <- with(df.cv1se, loess(Noise ~ Estimate, span = .05))
yy.cv1se <- predict(l.cv1se, newdata = xx)
smooth.cv1se <- data.frame(xx = xx, yy = yy.cv1se)

l.mfdr <- with(df.mfdr, loess(Noise ~ Estimate, span = .05))
yy.mfdr <- predict(l.mfdr, newdata = xx)
smooth.mfdr <- data.frame(xx = xx, yy = yy.mfdr)

l.uni <- with(df.uni, loess(Noise ~ Estimate, span = .55))
yy.uni <- predict(l.uni, newdata = xx)
smooth.uni <- data.frame(xx = xx, yy = yy.uni)

l.uniash <- with(df.uniash, loess(Noise ~ Estimate, span = .1))
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


### Now generate upper row
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


### Supplemental figure 2 
png("supplemental_figures_tables/Supp2.png", h=6, w=8.5, units = 'in', res = 300)
grid.arrange(arrangeGrob(p.cv, top="lasso mfdr (CV)"), arrangeGrob(p.cv1se, top="lasso mfdr (CV)"),  
             arrangeGrob(p.uniash, top="Univariate fdr", right = "Assumptions Violated"),
             arrangeGrob(p.cv2), arrangeGrob(p.cv1se2),  
             arrangeGrob(p.uniash2, right = "Assumptions Met"), ncol=3)
dev.off()

##############################################################
###
###  Supplemental Figure #3 (Power Logistic)
###
##############################################################

load("simres/logistic_hard.RData")

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
                             Method = c(rep("lasso mfdr (CV)", 3), rep("lasso mfdr (CV1se)", 3), rep("Lasso (mFDR)", 3),
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


load("simres/logistic_easy.RData")
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
                             Method = c(rep("lasso mfdr (CV)", 2), rep("lasso mfdr (CV1se)", 2), rep("Lasso (mFDR)", 2),
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
power.combined <- power.combined[which(power.combined$Method %in% c("lasso mfdr (CV)", "lasso mfdr (CV1se)", "Univariate (locfdr)")),]

png("supplemental_figures_tables/Supp3.png", h=4, w=6, units = 'in', res = 300)
p <- ggplot(power.combined, aes(x = relevel(Method, ref = "lasso mfdr (CV)"), y = Val, fill = Var)) + ggtitle("Logistic Regression") +
  geom_bar(stat = "identity", position=position_stack(reverse=TRUE)) + scale_fill_manual(name = "Feature Type", values = pal(3)[c(2,3,1)]) +
  coord_flip() + scale_y_continuous(expression(paste("Features with ", widehat(mfdr), " or ", widehat(fdr), " < 0.10"))) + scale_x_discrete(" ")
p + facet_grid(.~ Scenario, scales = "free") + theme(axis.text.y = element_text(angle = 45, hjust = 1))
dev.off()

##############################################################
###
###  Supplemental Figure #4 (Power Cox)
###
##############################################################

load("simres/cox_hard.RData")

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
                             Method = c(rep("lasso mfdr (CV)", 3), rep("lasso mfdr (CV1se)", 3), rep("Lasso (mFDR)", 3),
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


load("simres/cox_easy.RData")
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
                             Method = c(rep("lasso mfdr (CV)", 2), rep("lasso mfdr (CV1se)", 2), rep("Lasso (mFDR)", 2),
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
power.combined <- power.combined[which(power.combined$Method %in% c("lasso mfdr (CV)", "lasso mfdr (CV1se)", "Univariate fdr")),]

png("supplemental_figures_tables/Supp4.png", h=4, w=6, units = 'in', res = 300)
p <- ggplot(power.combined, aes(x = relevel(Method, ref = "lasso mfdr (CV)"), y = Val, fill = Var)) + ggtitle("Cox Regression") +
  geom_bar(stat = "identity", position=position_stack(reverse=TRUE)) + scale_fill_manual(name = "Feature Type", values = pal(3)[c(2,3,1)]) +
  coord_flip() + scale_y_continuous(expression(paste("Features with ", widehat(mfdr), " or ", widehat(fdr), " < 0.10"))) + scale_x_discrete(" ")
p + facet_grid(.~ Scenario, scales = "free") + theme(axis.text.y = element_text(angle = 45, hjust = 1))
dev.off()


##############################################################
###
###  Supplemental Figure #5 (ashr vs. kernel - Met)
###
##############################################################

load("simres/linear_easy.RData")

## Setup results data.frames
df.cv <- data.frame(Estimate = c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])),
                    Noise = c(rep(0, 60*ncol(res.lassocv)),
                              rep(1, 540*ncol(res.lassocv))))   

df.cvkern <- data.frame(Estimate = c(as.vector(res.lassocv.kernel[1:60,]), as.vector(res.lassocv.kernel[61:600,])),
                        Noise = c(rep(0, 60*ncol(res.lassocv.kernel)),
                                  rep(1, 540*ncol(res.lassocv.kernel))))   
  
df.uniash <-  data.frame(Estimate = c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])),
                           Noise = c(rep(0, 60*ncol(res.uniash)),
                                     rep(1, 540*ncol(res.uniash))))     
df.unikern <-  data.frame(Estimate = c(as.vector(res.uni[1:60,]), as.vector(res.uni[61:600,])),
                      Noise = c(rep(0, 60*ncol(res.uni)),
                                rep(1, 540*ncol(res.uni))))     

### Loess curves
xx <- seq(0, 1, by = .01)

l.cv <- with(df.cv, loess(Noise ~ Estimate, span = .15))
yy.cv <- predict(l.cv, newdata = xx)
smooth.cv <- data.frame(xx = xx, yy = yy.cv)

l.cvk<- with(df.cvkern, loess(Noise ~ Estimate, span = .9))
yy.cv1k<- predict(l.cvk, newdata = xx)
smooth.cvk <- data.frame(xx = xx, yy = yy.cv1se)

l.uni <- with(df.unikern, loess(Noise ~ Estimate, span = .75))
yy.uni <- predict(l.uni, newdata = xx)
smooth.uni <- data.frame(xx = xx, yy = yy.uni)

l.uniash <- with(df.uniash, loess(Noise ~ Estimate, span = .25))
yy.uniash <- predict(l.uniash, newdata = xx)
smooth.uniash <- data.frame(xx = xx, yy = yy.uniash)

## There are too many points, only plot a sample for the scatterplot
sample.id <- sample(1:nrow(df.cv), size = 8000, replace = FALSE)

## Create scatterplots from random subset of points
p.cv2 <- ggplot(df.cv[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.cvk2 <- ggplot(df.cvkern[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.uni2 <- ggplot(df.unikern[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.uniash2 <- ggplot(df.uniash[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)


## Add loess curves from full sim results
p.cv2 <- p.cv2 + geom_path(data = smooth.cv, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1) +
  labs(y = "Proportion (noise)")

p.cvk2 <- p.cvk2 + geom_path(data = smooth.cvk, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1) +
  labs(y = "")

p.uni2 <- p.uni2 + geom_path(data = smooth.uni, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "")

p.uniash2 <- p.uniash2 + geom_path(data = smooth.uniash, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y ="Proportion (noise)")

png("supplemental_figures_tables/Supp5.png", h=6, w=8.5, units = 'in', res = 300)
grid.arrange(arrangeGrob(p.cv2, top="ashr"), arrangeGrob(p.cvk2, top="kernel", right = "lasso mfdr (CV)"),
             arrangeGrob(p.uniash2), arrangeGrob(p.uni2, right = "Univariate fdr"), ncol=2, top = "Assumptions Met Scenario")
dev.off()


##############################################################
###
###  Supplemental Figure #6 (ashr vs. kernel - Violated)
###
##############################################################


load("simres/linear_hard.RData")

## Setup results data.frames
df.cv <- data.frame(Estimate = c(as.vector(res.lassocv[1:60,]), as.vector(res.lassocv[61:600,])),
                    Noise = c(rep(0, 60*ncol(res.lassocv)),
                              rep(1, 540*ncol(res.lassocv))))   

df.cvkern <- data.frame(Estimate = c(as.vector(res.lassocv.kernel[1:60,]), as.vector(res.lassocv.kernel[61:600,])),
                        Noise = c(rep(0, 60*ncol(res.lassocv.kernel)),
                                  rep(1, 540*ncol(res.lassocv.kernel))))   

df.uniash <-  data.frame(Estimate = c(as.vector(res.uniash[1:60,]), as.vector(res.uniash[61:600,])),
                         Noise = c(rep(0, 60*ncol(res.uniash)),
                                   rep(1, 540*ncol(res.uniash))))     
df.unikern <-  data.frame(Estimate = c(as.vector(res.uni[1:60,]), as.vector(res.uni[61:600,])),
                          Noise = c(rep(0, 60*ncol(res.uni)),
                                    rep(1, 540*ncol(res.uni))))     

### Loess curves
xx <- seq(0, 1, by = .01)

l.cv <- with(df.cv, loess(Noise ~ Estimate, span = .15))
yy.cv <- predict(l.cv, newdata = xx)
smooth.cv <- data.frame(xx = xx, yy = yy.cv)

l.cvk<- with(df.cvkern, loess(Noise ~ Estimate, span = .45))
yy.cv1k<- predict(l.cvk, newdata = xx)
smooth.cvk <- data.frame(xx = xx, yy = yy.cv1se)

l.uni <- with(df.unikern, loess(Noise ~ Estimate, span = .55))
yy.uni <- predict(l.uni, newdata = xx)
smooth.uni <- data.frame(xx = xx, yy = yy.uni)

l.uniash <- with(df.uniash, loess(Noise ~ Estimate, span = .10))
yy.uniash <- predict(l.uniash, newdata = xx)
smooth.uniash <- data.frame(xx = xx, yy = yy.uniash)

## There are too many points, only plot a sample for the scatterplot
sample.id <- sample(1:nrow(df.cv), size = 8000, replace = FALSE)

## Create scatterplots from random subset of points
p.cv2 <- ggplot(df.cv[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.cvk2 <- ggplot(df.cvkern[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.uni2 <- ggplot(df.unikern[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)
p.uniash2 <- ggplot(df.uniash[sample.id,]) + geom_jitter(aes(x = Estimate, y = Noise), height = .1, alpha = .1)


## Add loess curves from full sim results
p.cv2 <- p.cv2 + geom_path(data = smooth.cv, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1) +
  labs(y = "Proportion (noise)")

p.cvk2 <- p.cvk2 + geom_path(data = smooth.cvk, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1) +
  labs(y = "")

p.uni2 <- p.uni2 + geom_path(data = smooth.uni, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y = "")

p.uniash2 <- p.uniash2 + geom_path(data = smooth.uniash, aes(x = xx, y = yy), lwd = 1.1) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", lwd = 1)  +
  labs(y ="Proportion (noise)")

png("supplemental_figures_tables/Supp6.png", h=6, w=8.5, units = 'in', res = 300)
grid.arrange(arrangeGrob(p.cv2, top="ashr"), arrangeGrob(p.cvk2, top="kernel", right = "lasso mfdr (CV)"),
             arrangeGrob(p.uniash2), arrangeGrob(p.uni2, right = "Univariate fdr"), ncol=2, top = "Assumptions Violated Scenario")
dev.off()
