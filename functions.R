

### Function to calculate the local mfdr from a ncvreg object at a specified lambda value
locmfdr <- function(fit, lambda){
  ### Setup
  X <- fit$X
  y <- fit$y
  n <- nrow(X)
  p <- ncol(X)
  S <- predict(fit, type = "nvars")
  step <- which(fit$lambda == lambda)
  
  ### Linear regression
  if(class(fit)[1] == "ncvreg"){
    if(fit$family == "gaussian"){
      ### Calculate z_j
      R <- y - predict(fit, X, lambda = lambda, type = "response")
      z <- (1/n)*t(X) %*% R + fit$beta[-1,step]
      sig.est <- sqrt(fit$loss[step]/(n - S[step] + 1))
      z <- z/(sig.est/sqrt(n))
      
      ### Calculate locfdrs
      f <- density(z)
      ff <- approxfun(f$x, f$y)
      est.gam <- pmin(dnorm(z, 0, 1)/ff(z), 1) ## rows = variables, cols = lambdas
    }
    
    ### Logistic regression 
    else if (fit$family == "binomial"){
      P <-  predict(fit, X, lambda = lambda, type = 'response')
      U <- y - P   ### Calculating the score
      W <- diag(as.vector(P*(1 - P)))  
      
      ### Calculate zj
      z <- rep(NA, p)
      for (j in 1:p){
        vj <- t(X[,j]) %*% W %*% X[,j]
        z[j] <- (X[,j] %*% U + vj * fit$beta[(j+1),step])/(sqrt(vj))  ### j+1 bc intercept
      }
      
      ### Calculate locfdrs
      f <- density(z)
      ff <- approxfun(f$x, f$y)
      est.gam <- pmin(dnorm(z, 0, 1)/ff(z), 1)
    }
  }
  
  ### Cox regression
  if(class(fit)[1] == "ncvsurv"){
    d <- fit$fail
    rsk <- rev(cumsum(rev(exp(fit$Eta[,step]))))
    P <- outer(exp(fit$Eta[,step]), rsk, '/')
    P[upper.tri(P)] <- 0
    U <- d - P%*%d  
    
    ### Calculate W (maybe diag for speed?)
    W <- -P %*% diag(d) %*% t(P)
    #W <- matrix(0, n,n)
    diag(W) <- diag(P %*% diag(d) %*% t(1-P))
    
    ### Calculate v_j and then z_j
    z <- rep(NA, p)
    for (j in 1:p){
      vj <- t(X[,j]) %*% W %*% X[,j]
      z[j] <- (X[,j] %*% U + vj * fit$beta[j,step])/(sqrt(vj))
    }
    
    ### Calculate locfdr
    f <- density(z)
    ff <- approxfun(f$x, f$y)
    est.gam <- pmin(dnorm(z, 0, 1)/ff(z), 1)
  }
  
  results <- data.frame(z = z, locfdr = est.gam)
  return(results)
}


### Function used to generate correlated data following structural shown in the A-B-C causal diagram
genData <- function(n, J, K=1, beta, family=c("gaussian","binomial","survival"), J0=ceiling(J/2), K0=K, SNR=1, sig = c("homogeneous","heterogeneous"), sig.g = c("homogeneous","heterogeneous"), rho = 0, rho.g = rho, corr=c("exchangeable", "autoregressive")) {
  family <- match.arg(family)
  sig <- match.arg(sig)
  sig.g <- match.arg(sig.g)
  corr <- match.arg(corr)
  
  ## Gen X, S
  if (corr=="exchangeable") {
    X <- genX(n=n, J=J, K=K, rho=rho, rho.g=rho.g)
  } else {
    RHO <- matrix(rho^(0:(J-1)), J, J, byrow=TRUE)
    S <- bandSparse(J, k=0:(J-1), diagonals=RHO, symmetric=TRUE)
    R <- chol(S)
    X <- as.matrix(matrix(rnorm(n*J), n, J) %*% R)
  }
  
  j <- rep(1:J,rep(K,J))
  
  ## Gen beta
  if (missing(beta) || length(beta)==1) {
    k <- rep(1:K,J)
    b <- (j <= J0) * (k <= K0)
    s <- c(1,-1)[1+j%%2] * c(1,-1)[1+k%%2]
    if (missing(beta)) {
      S <- matrix(rho, nrow=J*K, ncol=J*K)
      for (i in 1:J) S[(i-1)*K+1:K,(i-1)*K+1:K] <- rho.g
      diag(S) <- rep(1,J*K)
      if (sig=="heterogeneous") b <- b*j
      if (sig.g=="heterogeneous") b <- b*k
      b <- b*s
      beta <- b*sqrt(SNR)/sqrt(crossprod(b,S)%*%b)
    } else beta <- b*s*beta
  }
  
  ## Gen y
  y <- genY(X%*%beta, family=family, sigma=sqrt(n))
  return(list(X=X,y=y,beta=beta,family=family,group=j))
}

## rho  : correlation across all explanatory variables
## rho.g: correlation within group (must be at least rho)
genX <- function(n, J, K=1, rho=0, rho.g=rho, corr=corr) {
  a <- sqrt(rho/(1-rho.g))
  b <- sqrt((rho.g-rho)/(1-rho.g))
  Z <- rnorm(n)
  ZZ <- t(matrix(rep(rnorm(n*J), rep(K,n*J)), ncol=n))
  ZZZ <- matrix(rnorm(n*J*K),nrow=n)
  return(matrix(as.numeric(a*Z + b*ZZ + ZZZ),nrow=n)/sqrt(1+a^2+b^2))
}

genY <- function(eta,family=c("gaussian","binomial","survival"),sigma=1,lam.0=1) {
  family=match.arg(family)
  n <- length(eta)
  if (family=="gaussian") y <- rnorm(n,mean=eta,sd=sigma)
  else if (family=="binomial")
  {
    pi. <- exp(eta)/(1+exp(eta))
    pi.[eta > log(.9999/.0001)] <- 1
    pi.[eta < log(.0001/.9999)] <- 0
    y <- rbinom(n,1,pi.)
  } else if (family == "survival")
  {
    haz <- lam.0*exp(eta)
    y <- rexp(n, haz)
  }
  return(y)
}


### A version of locmfdr with more features (sorting, in-house standardization checks, etc)
locmfdr.new <- function(fit, lambda, X = NULL, y = NULL, number = NULL, alpha = NULL){
  ### Check for valid args
  if(!is.null(number)){
    if(number < 1) stop("'number' should be a positive integer")
  }
  if(!is.null(alpha)){
    if(alpha > 1 | alpha <= 0) stop("'alpha' should be in the interval (0,1]")
  }
  
  ### Setup standardized X
  if(is.null(X)){
    if(is.null(fit$X)){
      stop("This procedure requires X and y, either supply X and y, or fit the model using the option 'returnX = TRUE'")
    } else {XX <- fit$X
    y <- fit$y} 
  } else {
    if(class(fit)[1] == "ncvsurv"){XX <- std(X[fit$order,])}
    else {XX <- std(X)}
  }
  
  
  ### Setup general
  lid <- which(fit$lambda == lambda)
  ns <- attr(XX, "nonsingular")
  sc <- attr(XX, "scale")
  n <- nrow(XX)
  p <- ncol(XX)
  S <- predict(fit, type = "nvars", lambda = lambda)
  pen.idx <- fit$penalty.factor > 0
  
  ##### Linear Regression
  if(class(fit)[1] == "ncvreg"){
    if(fit$family == "gaussian"){
      ### Setup standardized beta and centered y
      yy <- y - mean(y)
      bb <- c(mean(y), fit$beta[ns+1,lid]*sc)
      
      ### Calculate standardized z_j's
      R <- yy - cbind(1, XX) %*% bb
      z <- (1/n)*t(XX) %*% R + bb[-1]
      sig.est <- sqrt(fit$loss[lid]/(n - S + 1))
      z <- z/(sig.est/sqrt(n))
    }
    
    ### Logistic regression 
    else if (fit$family == "binomial"){
      ## Setup standardized beta
      if(is.factor(y)){
        y <- as.numeric(y)-1
      }
      bb <- c(mean(y), fit$beta[ns+1,lid]*sc)
      
      ### Setup the score vector and W matrix
      P <-  predict(fit, X, lambda = lambda, type = 'response')
      U <- y - P   
      W <- diag(as.vector(P*(1 - P)))  
      
      ### Calculate v_j and z_j
      z <- numeric(p)
      for (j in 1:p){
        vj <- t(XX[,j]) %*% W %*% XX[,j]
        z[j] <- (XX[,j] %*% U + vj * bb[j+1])/(sqrt(vj))  ### j+1 bc intercept
      }
    }
  }
  
  
  ### Cox regression
  if(class(fit)[1] == "ncvsurv"){
    ## Setup standardized beta
    bb <- fit$beta[ns,lid]*sc
    
    ### Calculate score vector and W (maybe diagonalize it for speed?)
    d <- fit$fail
    rsk <- rev(cumsum(rev(exp(fit$Eta[,lid]))))
    P <- outer(exp(fit$Eta[,lid]), rsk, '/')
    P[upper.tri(P)] <- 0
    U <- d - P%*%d  
    W <- -P %*% diag(d) %*% t(P)
    #W <- matrix(0,n,n)
    diag(W) <- diag(P %*% diag(d) %*% t(1-P))
    
    ### Calculate v_j and z_j
    z <- numeric(p)
    for (j in 1:p){
      vj <- t(XX[,j]) %*% W %*% XX[,j]
      z[j] <- (XX[,j] %*% U + vj * bb[j])/(sqrt(vj))
    }
  }
  
  ### Calculate locfdr
  f <- density(z[pen.idx])
  ff <- approxfun(f$x, f$y)
  est.gam <- pmin(dnorm(z, 0, 1)/ff(z), 1)
  est.gam[!pen.idx] <- NA
  
  
  #### Calculate Fdr (using both est cdf and empirical cdf)
  est.Fdr <- emp.Fdr <- numeric(p)
  est.Fdr[!pen.idx] <- emp.Fdr[!pen.idx] <- NA
  
  ### setup results
  
  if(class(fit)[1] == "ncvreg"){
    results <- data.frame(beta = fit$beta[-1,lid], z = z, locfdr = est.gam, est.mFdr = est.Fdr, emp.mFdr = emp.Fdr)
    results$selected <- ifelse(fit$beta[-1,lid] != 0, "*"," ")
    rownames(results) <- rownames(fit$beta)[-1]
  } else {
    results <- data.frame(beta = fit$beta[,lid], z = z, locfdr = est.gam, est.mFdr = est.Fdr, emp.mFdr = emp.Fdr)
    results$selected <- ifelse(fit$beta[,lid] != 0, "*"," ")
    rownames(results) <- rownames(fit$beta)
  }
  
  return(results[pen.idx,])
}



## Case Study mfdr function
#### DATA MUST BE STANDARDIZED BEFORE MODEL FITTING WITH THIS FUCTION
case.locmfdr <- function(fit, lambda){
  ### Setup
  X <- fit$X
  y <- fit$y
  n <- nrow(X)
  p <- ncol(X)
  S <- predict(fit, type = "nvars")
  step <- which(fit$lambda == lambda)
  pen.idx <- fit$penalty.factor > 0
  
  ### Linear regression
  if(class(fit)[1] == "ncvreg"){
    if(fit$family == "gaussian"){
      ### Calculate z_j
      R <- y - predict(fit, X, lambda = lambda, type = "response")
      z <- (1/n)*t(X) %*% R + fit$beta[-1,step]
      sig.est <- sqrt(fit$loss[step]/(n - S[step] + 1))
      z <- z/(sig.est/sqrt(n))
      
    }
    
    ### Logistic regression 
    else if (fit$family == "binomial"){
      P <-  predict(fit, X, lambda = lambda, type = 'response')
      U <- y - P   ### Calculating the score
      W <- diag(as.vector(P*(1 - P)))  
      
      ### Calculate zj
      z <- rep(NA, p)
      for (j in 1:p){
        vj <- t(X[,j]) %*% W %*% X[,j]
        z[j] <- (X[,j] %*% U + vj * fit$beta[(j+1),step])/(sqrt(vj))  ### j+1 bc intercept
      }
      
    }
  }
  
  ### Cox regression
  if(class(fit)[1] == "ncvsurv"){
    d <- fit$fail
    rsk <- rev(cumsum(rev(exp(fit$Eta[,step]))))
    P <- outer(exp(fit$Eta[,step]), rsk, '/')
    P[upper.tri(P)] <- 0
    U <- d - P%*%d  
    
    ### Calculate W (maybe diag for speed?)
    W <- -P %*% diag(d) %*% t(P)
    #W <- matrix(0, n,n)
    diag(W) <- diag(P %*% diag(d) %*% t(1-P))
    
    ### Calculate v_j and then z_j
    z <- rep(NA, p)
    for (j in 1:p){
      vj <- t(X[,j]) %*% W %*% X[,j]
      z[j] <- (X[,j] %*% U + vj * fit$beta[j,step])/(sqrt(vj))
    }
    
  }
  
  ### Calculate locfdr
  
  f <- density(z[pen.idx])
  ff <- approxfun(f$x, f$y)
  est.gam1 <- pmin(dnorm(z, 0, 1)/ff(z), 1)
  est.gam1[!pen.idx] <- NA
  
  
  results <- data.frame(z = z, locfdr_def = est.gam1)
  #return(list(results, f, aa, bkde.xy)) ### Use to generate comparison plots
  return(results)
}



## Color palatte used for the figures
pal <- function(n, alpha=1) {
  if (n==2) {
    val <- hcl(seq(15,375,len=4), l=60, c=150, alpha=alpha)[c(1,3)]
  } else {
    val <- hcl(seq(15,375,len=n+1), l=60, c=150, alpha=alpha)[1:n]
  }
  val
}

### Function to stack multiple plots
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
