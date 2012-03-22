###################
# Content Listing #
###################

# Generic & internal functions:
#  - crit.fun()
#  - bwmc() [C wrapper]

# Preprocessing functions:
#  - makeMyD
#  - initializeHP (contains devHelper, checkMyData and getMclustHPests)

# Run functions:
#  - ebCoexpressZeroStep   [Contains C wrapper]
#  - ebCoexpressOneStep    [Contains C wrapper]
#  - ebCoexpressFullTCAECM [Contains C wrapper]
#  - ebCoexpressMeta       [Contains C wrapper]

# Show functions:
#  - showPair
#  - showNetwork
#  - rankMyGenes

####################
# Loading packages #
####################

library(EBarrays)
library(mclust)
library(minqa)

##################################
# Generic and internal functions #
##################################

crit.fun <- function (ecPostProbs, targetFDR=0.05) 
{
  x <- ecPostProbs ; cc <- targetFDR
  yy <- cumsum(sort(x))/(1:length(x))
  mm <- yy < cc ; index <- sum(mm)
  if(index > 0)   out <- 1 - sort(x)[index]
  if(index == 0)   out <- 1
  names(out) <- NULL
  return(out)
}

bwmc <- function(X)
{
  cn <- dim(X)[2] ; rn <- dim(X)[1]
  meds <- apply(X,2,FUN=median)
  mads <- apply(X,2,FUN=mad,constant=1)
  return(as.matrix(.Call("bwmcCworker",Xx=X,rnx=rn, cnx=cn,
                         medsx=meds, madsx=mads, PACKAGE="EBcoexpress"),rn,cn))
}

##################
# User functions #
##################

makeMyD <- function(X, conditions, useBWMC=FALSE, gpsep="~")
{
  k <- length(unique(conditions))
  m <- dim(X)[1] ; labels <- 1:m
  if(!is.null(dimnames(X)[[1]]))   labels <- dimnames(X)[[1]]
  p <- choose(m,2)
  
  D <- matrix(NaN, p, k)
  ii <- rep.int(seq_len(m-1), times=rev(seq_len(m-1)))
  foo <- function(x) {return(x:m)}
  jj <- unlist(sapply(2:m, FUN=foo))
  dimnames(D)[[1]] <- mapply(labels[ii],gpsep,labels[jj],FUN=paste,sep="")
  dimnames(D)[[2]] <- paste("Condition",1:k,sep="")

  for(j in 1:k) # Running through conditions
  {
    cases <- t(X[,conditions==j]) # unns[j] by m
    if(useBWMC)   casx <- bwmc(cases)
    if(!useBWMC)  casx <- cor(cases)
    D[,j] <- casx[lower.tri(casx)]
    if(useBWMC)
    {
      if(sum(D[,j]==Inf) > 0)
              stop(paste("Infinite BWMC values detected in condition ",j,
                   "!\nPlease check for large quantities of non-unique ",
                   "values, or set useBWMC= to FALSE",sep=""))
      if(sum(is.na(D[,j])) > 0)
              stop(paste("0/0 BWMC values detected in condition ",j,
                   "!\nPlease check for large quantities of non-unique ",
                   "values, or set useBWMC= to FALSE",sep=""))
    }
  }
  return(D)
}

initializeHP <- function(D, conditions, seed=NULL, plottingOn=FALSE,
                         colx="red", applyTransform=TRUE,
                         verbose=0, subsize=NULL, ...)
{

  getMclustHPests <- function(D, conditions,
                          seed=NULL, adjustForUncertainty=TRUE,
                          plottingOn=FALSE, colx="red",
                          applyTransform=TRUE, withDev=FALSE,
                          verbose=0, ...)
  {

    devHelper <- function(dat, est)
    {
      dens <- density(dat)
      fundDev <- function(x)
      {
        eG <- est[[1]] ; eMus <- est[[2]] ; eTaus <- est[[3]]
        eWs <- est[[4]] ; eWs <- eWs/sum(eWs)
        valE <- 0
        for(i in 1:eG)
          valE <- valE + eWs[i]*dnorm(x, eMus[i], eTaus[i])
        return(valE)
      }
      fit <- sapply(dens$x, FUN=fundDev)
      dev <- sum((dens$y - fit)^2)
      return(dev)
    }

    k <- length(unique(conditions)) # No. of conditions
    unns <- array(NaN, k)
    for(j in 1:k)
      unns[j] <- sum(conditions==j) # No. of samples per condition
    wns <- floor(mean(unns))

    fZ <- function(rho) { return(0.5*log((1+rho)/(1-rho))) }
    if(applyTransform)  D <- fZ(D)

    dat <- c(D)
    if(!is.null(seed)) set.seed(seed)
    out <- suppressWarnings(Mclust(dat, G=1:3, modelNames="V"))
    eG <- out$G
    eMus <- out$parameters$mean
    eVars <- out$parameters$variance$sigmasq
    eBaus <- sqrt(eVars)
    eTaus <- suppressWarnings(sqrt(eVars - 1/(wns-3)))
    if(any(!is.finite(eTaus)))
    {
      if(verbose > 0)
      {
        cat("Negative tau estimates ...\n")
        cat("Changing correction scheme ...\n")
      }
      eTaus <- eBaus/2
    }
    clas <- out$classification
    eWs <- rep(0,eG)
    if(!adjustForUncertainty)
      for(i in 1:eG)
        eWs[i] <- sum(clas==i)
    if(adjustForUncertainty)
      for(i in 1:eG)
        eWs[i] <- sum((1-out$uncert)[clas==i])
    eWs <- eWs/sum(eWs)

    ord <- order(eMus, decreasing=FALSE)
    eMus <- eMus[ord]
    eBaus <- eBaus[ord]
    eTaus <- eTaus[ord]
    eWs <- eWs[ord]

    checkMyData <- function(dat, G, mus, taus, weights=NULL, colx=NULL, ...)
    {
      if(is.null(weights)) weights <- rep(1,G)
      if(is.null(colx)) colx <- "blue"
      p <- length(dat)
      weights <- weights/sum(weights)
      xs <- seq(min(dat),max(dat),length.out=5000)
      funa <- function(x)
      {
        tot <- 0
        for(i in 1:G)
          tot <- tot + weights[i]*dnorm(x, mus[i], taus[i])
        return(tot)
      }
      ys <- sapply(xs, FUN=funa)
      yse <- density(dat)$y
      ymaxa <- max(max(ys),max(yse))
      plot(density(dat), ylim=c(0,ymaxa), ...)
      lines(xs, ys, lty=2, col=colx)
      return(invisible(NULL))
    }

    estObj <- list(G=eG, MUS=eMus, TAUS=eTaus, WEIGHTS=eWs)
    if(plottingOn)
      checkMyData(dat, eG, eMus, eBaus, eWs, colx=colx, ...)
    if(!withDev)
      return(estObj)
    biaObj <- list(eG, eMus, eBaus, eWs)
    dev <- devHelper(dat, biaObj)
    return(list(estObj,dev))
  } # End of getMclustHPests

  if(!is.null(subsize))
  {
    p <- dim(D)[[1]]
    if(p < subsize)    stop("Unusable size for subset of D given!")
    if(!is.null(seed)) set.seed(seed)
    D <- D[sample(1:p, subsize, replace=FALSE),]
  }

  if(plottingOn) par(mfrow=c(2,1))
  ost1 <- getMclustHPests(D, conditions, seed, plottingOn, colx,
                          applyTransform, adjustForUn=FALSE, withDev=TRUE,
                          main="Unweighted Mclust Fit", verbose=verbose)
  ost2 <- getMclustHPests(D, conditions, seed, plottingOn, colx,
                          applyTransform, adjustForUn=TRUE, withDev=TRUE,
                          main="Weighted Mclust Fit", verbose=verbose)
  if(plottingOn) par(mfrow=c(1,1))
  est1 <- ost1[[1]]
  dev1 <- ost1[[2]]
  est2 <- ost2[[1]]
  dev2 <- ost2[[2]]

  if(verbose > 0)
  {
    cat(paste("U-Deviance:",round(dev1,5),"\n"))
    cat(paste("W-Deviance:",round(dev2,5),"\n"))
  }
  if(dev1 < dev2)
  {
    if(verbose > 0) cat("Returning Unweighted Estimate ...\n")
    return(est1)
  }
  if(dev1 >= dev2)
  {
    if(verbose > 0) cat("Returning Weighted Estimate ...\n")
    return(est2)
  }
}

priorDiagnostic <- function(D, conditions, ebOutObj, focusCond,
                         seed=NULL, colx="red", applyTransform=TRUE,
                         subsize=NULL, ...)
{

  if(!is.null(subsize))
  {
    p <- dim(D)[[1]]
    if(p < subsize)    stop("Unusable size for subset of D given!")
    if(!is.null(seed)) set.seed(seed)
    D <- D[sample(1:p, subsize, replace=FALSE),]
  }

  k <- length(unique(conditions)) # No. of conditions
  if(!(focusCond %in% 1:k))
    stop("Error! The parameter focusCond must be a value in 1:K!")
  unns <- array(NaN, k)
  for(j in 1:k)
    unns[j] <- sum(conditions==j) # No. of samples per condition
  nf <- unns[focusCond]

  fZ <- function(rho) { return(0.5*log((1+rho)/(1-rho))) }
  if(applyTransform)  D <- fZ(D)

  dat <- D[,focusCond]
  if(!is.null(seed)) set.seed(seed)
  params <- ebOutObj$MODEL$HPS
  eG <- params$G
  eMUS <- params$MUS
  eTAUS <- params$TAUS
  eWEIS <- params$WEIGHTS

  p <- length(dat)
  xs <- seq(min(dat),max(dat),length.out=5000)
  funa <- function(x)
  {
    tot <- 0
    for(i in 1:eG)
      tot <- tot + eWEIS[i]*dnorm(x, eMUS[i], sqrt(eTAUS[i]^2+(1/(nf-3))))
    return(tot)
  }
  ys <- sapply(xs, FUN=funa)
  yse <- density(dat)$y
  ymaxa <- max(max(ys),max(yse))
  plot(density(dat), ylim=c(0,ymaxa), lwd=2,
    main=paste("Diagnostic using Condition",focusCond),
    xlab="Transformed Correlation Space", ylab="Density", ...)
  lines(xs, ys, lty=2, col=colx, lwd=2)
  return(invisible(NULL))
}

############

ebCoexpressZeroStep <- function(D, conditions, pattern, hpEsts,
                        controlOptions=list())
{

bool <- (sort(unique(conditions)) == sort(unique(unlist(pattern@patterns))))
if(sum(!bool) != 0)
  stop("Conditions and pattern are inconsistent!")
controlList <- list(applyTransform=TRUE, verbose=1)
cOnames <- names(controlOptions)
if(!is.null(cOnames))
{
  if("verbose" %in% cOnames)
    controlList$verbose <- controlOptions$verbose
  if("applyTransform" %in% cOnames)
    controlList$applyTransform <- controlOptions$applyTransform
}

n <- length(conditions)         # No. of samples
k <- length(unique(conditions)) # No. of conditions
unns <- array(NaN, k)
for(j in 1:k)
  unns[j] <- sum(conditions==j) # No. of samples per condition
L <- length(pattern@patterns)   # No. of patterns

fZ <- function(rho) { return(0.5*log((1+rho)/(1-rho))) }
p <- dim(D)[1]
if(dim(D)[2] != k)
  stop("2nd dimension of D does not match number of unique conditions!")
if(controlList$applyTransform)
  D <- fZ(D)

verbose <- controlList$verbose

### End of Setup ###

# Hyperparam setup

G <- hpEsts[[1]]
eMus <- hpEsts[[2]]
eTaus <- hpEsts[[3]]
eWs <- hpEsts[[4]]
eWs <- eWs/sum(eWs)

theta <- NULL
if(G==1)
  { theta <- c(eMus, eTaus) ; tl <- 2 } # Length 2
if(G==2)
  { theta <- c(eMus, eTaus, eWs[1]) ; tl <- 5 } # Length 5
if(G==3)
  { theta <- c(eMus, eTaus, eWs[1:2]) ; tl <- 8 } # Length 8
if(!(G %in% c(1,2,3)))
  stop("Invalid hpEsts object!")

# EM functions

f0toC <- function(zvals, sizes, theta)
{
  return(as.numeric(.Call("f0worker",zvalsx=zvals,sizes3x=(sizes-3),
                          thetax=theta, dx=length(zvals), PACKAGE="EBcoexpress")))
}

theGManager <- function(ind, theta, pat)
{
  if(G==1) return(theGWorker1(ind, theta, pat))
  if(G==2) return(theGWorker2(ind, theta, pat))
  if(G==3) return(theGWorker3(ind, theta, pat))
  stop("Unreachable line reached in G Manager!")
}

theGWorker1 <- function(ind, theta, pat)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      zs <- D[ind,subInd]
      szs <- unns[subInd]
      df0 <- df0 + f0toC(zvals=zs,sizes=szs,theta)
    }
  return(df0)
}

theGWorker2 <- function(ind, theta, pat)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      w2 <- (1-theta[5])
      zs <- D[ind,subInd]
      szs <- unns[subInd]
      l1 <- theta[5]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[1],theta[3])))
      l2 <- w2*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[2],theta[4])))
      df0 <- df0 + log(l1+l2)
    }
  return(df0)
}

theGWorker3 <- function(ind, theta, pat)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      w3 <- (1-theta[7]-theta[8])
      zs <- D[ind,subInd]
      szs <- unns[subInd]
      l1 <- theta[7]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[1],theta[4])))
      l2 <- theta[8]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[2],theta[5])))
      l3 <- w3*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[3],theta[6])))
      df0 <- df0 + log(l1+l2+l3)
    }
  return(df0)
}

# Direct zero-step EM ppbs calculation using predefined mix

mix <- c(0.8,rep((1-0.8)/(L-1), L-1)) # Floor mix
ppbs <- matrix(0,p,L) # ppbs is #pairs by #patterns

time0 <- proc.time()
  for(j in 1:L)
    ppbs[,j] <- exp(log(mix[j]) +
                  sapply(1:p, FUN=theGManager, theta=theta, pat=j))
  ppbs <- ppbs/(apply(ppbs,1,FUN=sum))
  dimnames(ppbs)[[1]] <- dimnames(D)[[1]]
  dimnames(ppbs)[[2]] <- c("EC",paste("DC",(1:(L-1)),sep=""))
time1 <- proc.time()

if(verbose > 0)
  cat(paste("Zero-Stepper Time:",round((time1-time0)[3],5),"\n"))

out <- list(MODEL=list(MIX=mix, HPS=hpEsts), POSTPROBS=ppbs)
return(out)
}

###
###
###

ebCoexpressOneStep <- function(D, conditions, pattern, hpEsts,
                        controlOptions=list())
{

bool <- (sort(unique(conditions)) == sort(unique(unlist(pattern@patterns))))
if(sum(!bool) != 0)
  stop("Conditions and pattern are inconsistent!")
controlList <- list(convtol=0.0005, enforceFloor=TRUE, applyTransform=TRUE,
                    verbose=1, subsize=NULL, m2MaxIter=100)
cOnames <- names(controlOptions)
if(!is.null(cOnames))
{
  if("verbose" %in% cOnames)
    controlList$verbose <- controlOptions$verbose
  if("convtol" %in% cOnames)
    controlList$convtol <- controlOptions$convtol
  if("enforceFloor" %in% cOnames)
    controlList$enforceFloor <- controlOptions$enforceFloor
  if("applyTransform" %in% cOnames)
    controlList$applyTransform <- controlOptions$applyTransform
  if("subsize" %in% cOnames)
    controlList$subsize <- controlOptions$subsize
  if("m2MaxIter" %in% cOnames)
    controlList$m2MaxIter <- controlOptions$m2MaxIter
}

n <- length(conditions)         # No. of samples
k <- length(unique(conditions)) # No. of conditions
unns <- array(NaN, k)
for(j in 1:k)
  unns[j] <- sum(conditions==j) # No. of samples per condition
L <- length(pattern@patterns)   # No. of patterns

fZ <- function(rho) { return(0.5*log((1+rho)/(1-rho))) }
p <- dim(D)[1]
if(dim(D)[2] != k)
  stop("2nd dimension of D does not match number of unique conditions!")
if(controlList$applyTransform)
  D <- fZ(D)

convtol <- controlList$convtol
enforceFloor <- controlList$enforceFloor
verbose <- controlList$verbose
m2MaxIter <- controlList$m2MaxIter
subp <- controlList$subsize
subFlag <- FALSE
if(!is.null(subp))
  subFlag <- TRUE
if(subFlag)
  if(subp > p)
    stop("The provided subsize is greater than the number of pairs!")

### End of Setup ###

# Hyperparam setup

G <- hpEsts[[1]]     # It is critically vital that this be visible
eMus <- hpEsts[[2]]
eTaus <- hpEsts[[3]]
eWs <- hpEsts[[4]]
eWs <- eWs/sum(eWs)

theta <- NULL
if(G==1)
  { theta <- c(eMus, eTaus) ; tl <- 2 } # Length 2
if(G==2)
  { theta <- c(eMus, eTaus, eWs[1]) ; tl <- 5 } # Length 5
if(G==3)
  { theta <- c(eMus, eTaus, eWs[1:2]) ; tl <- 8 } # Length 8
if(!(G %in% c(1,2,3)))
  stop("Invalid hpEsts object!")

mix <- c(0.9,rep((1-0.9)/(L-1), L-1)) # Realistic mix
Z <- matrix(0,p,L) # Z is #pairs by #patterns
Zra <- NULL
randa <- NULL
if(subFlag)
  randa <- sort(sample(1:p, subp, replace=FALSE))

# EM functions

f0toC <- function(zvals, sizes, theta)
{
  return(as.numeric(.Call("f0worker",zvalsx=zvals,sizes3x=(sizes-3),
                          thetax=theta, dx=length(zvals), PACKAGE="EBcoexpress")))
}

theMiniLik <- function(theta)
{
  if(G==3)
    if(theta[7]+theta[8] > 1)   # This is the constraint for G=3
      return(Inf)
  tot <- 0
  # Note: randa (length=subp), a subset of 1:p, must be visible,
  #  as must Zra <- Z[randa,]
  pred.dens2 <- Zra # Something subp-by-L
  for(ind in 1:length(randa))
    for(pat in 1:L)
      pred.dens2[ind,pat] <- theGManager(randa[ind], theta, pat)
  tot <- sum(Zra * pred.dens2)
  return(-tot)
}

theGLik <- function(theta)
{
  if(G==3)
    if(theta[7]+theta[8] > 1)   # This is the constraint for G=3
      return(Inf)
  tot <- 0
  pred.dens <- Z # Something p-by-L
  for(ind in 1:p)
    for(pat in 1:L)
      pred.dens[ind,pat] <- theGManager(ind, theta, pat)
  tot <- sum(Z * pred.dens)
  return(-tot)
}

theGManager <- function(ind, theta, pat)
{
  if(G==1) return(theGWorker1(ind, theta, pat))
  if(G==2) return(theGWorker2(ind, theta, pat))
  if(G==3) return(theGWorker3(ind, theta, pat))
  stop("Unreachable line reached in G Manager!")
}

theGWorker1 <- function(ind, theta, pat)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      zs <- D[ind,subInd]
      szs <- unns[subInd]
      df0 <- df0 + f0toC(zvals=zs,sizes=szs,theta)
    }
  return(df0)
}

theGWorker2 <- function(ind, theta, pat)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      w2 <- (1-theta[5])
      zs <- D[ind,subInd]
      szs <- unns[subInd]
      l1 <- theta[5]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[1],theta[3])))
      l2 <- w2*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[2],theta[4])))
      df0 <- df0 + log(l1+l2)
    }
  return(df0)
}

theGWorker3 <- function(ind, theta, pat)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      w3 <- (1-theta[7]-theta[8])
      zs <- D[ind,subInd]
      szs <- unns[subInd]
      l1 <- theta[7]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[1],theta[4])))
      l2 <- theta[8]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[2],theta[5])))
      l3 <- w3*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[3],theta[6])))
      df0 <- df0 + log(l1+l2+l3)
    }
  return(df0)
}

# Other setup for the run:

mtol  <- convtol
done <- FALSE
iter <- 1
omix <- rep(-1,length(mix))

########################################
# Here is the Plan: E M2 [E M1] RET(Z) #
########################################

time0 <- proc.time()

#@# Step 1: E-Step

if(verbose > 0)
  cat("Begin Phase I (Initial E-Step) ...\n")

for(j in 1:L)
  Z[,j] <- exp(log(mix[j]) +
               sapply(1:p, FUN=theGManager, theta=theta, pat=j))
Z <- Z/(apply(Z,1,FUN=sum))
if(subFlag)
  Zra <- Z[randa,]

#@# Step 2: M2-Step

if(verbose > 0)
  cat("Begin Phase II (M2-Step) ...\n")

lowlim <- NULL
uplim <- NULL
if(G==1)
  {
    lowlim <- c(-Inf, 0)
    uplim <- c(Inf, Inf)
  }
if(G==2)
  {
    lowlim <- c(-Inf, -Inf, 0, 0, 0)
    uplim <- c(Inf, Inf, Inf, Inf, 1)
  }
if(G==3)
  {
    lowlim <- c(-Inf, -Inf, -Inf, 0, 0, 0, 0, 0)
    uplim <- c(Inf, Inf, Inf, Inf, Inf, Inf, 1, 1)
  }

con <- list()
if(!is.null(m2MaxIter))
  con <- list(maxfun=m2MaxIter)

if(subFlag)
{
  resu <- suppressWarnings(bobyqa(theta, theMiniLik, lower=lowlim, upper=uplim, control=con))
  theta <- resu$par
}
if(!subFlag)
{
  resu <- suppressWarnings(bobyqa(theta, theGLik, lower=lowlim, upper=uplim, control=con))
  theta <- resu$par
}

#@# Step 3: [E M1] Cycle

if(verbose > 0)
  cat("Begin Phase III ([E M1] Cycle) ...\n")

ZV <- matrix(0, p, L)
for(j in 1:L)
  ZV[,j] <- sapply(1:p, FUN=theGManager, theta=theta, pat=j)

while(!done)
{
  if(verbose > 0)
    cat(paste("Iteration:",iter,"\n"))

  for(j in 1:L)
    Z[,j] <- exp(log(mix[j]) + ZV[,j])
  Z <- Z/(apply(Z,1,FUN=sum))
  p.temp <- apply(Z,2,sum)
  mix <- p.temp/sum(p.temp)

  # Checking for convergence

  mis <- array(NaN, L)
  for(j in 1:L)
  {
    mis[j] <- abs(mix[j] - omix[j])
    if(!is.finite(mis[j])) mis[j] <- 0
  }

  done <- TRUE
  for(j in 1:L)
    if(mis[j] >= convtol)
      done <- FALSE
  iter <- iter +1
  if(enforceFloor)
  {
    if(mix[1] < 0.8)
    {
      if(verbose>1)
        cat("Floor enforced ...\n")
      if(L==2) mix <- c(0.8, 0.2)
      if(L>2)
      {
        ix <- mix[-1]
        ix <- ix*(0.2/sum(ix))
        mix <- c(0.8, ix)
      }
      for(j in 1:L)
        Z[,j] <- exp(log(mix[j]) + ZV[,j])
      Z <- Z/(apply(Z,1,FUN=sum))
      done <- TRUE
    }
  }
  omix <- mix
}

dimnames(Z)[[1]] <- dimnames(D)[[1]]
dimnames(Z)[[2]] <- c("EC",paste("DC",(1:(L-1)),sep=""))

hpObj <- NULL
if(G==1)
  hpObj <- list(G=G, MUS=theta[1], TAUS=theta[2], WEIGHTS=1)
if(G==2)
  hpObj <- list(G=G, MUS=theta[1:2], TAUS=theta[3:4],
                WEIGHTS=c(theta[5], 1-theta[5]))
if(G==3)
  hpObj <- list(G=G, MUS=theta[1:3], TAUS=theta[4:6],
                WEIGHTS=c(theta[7:8], 1-sum(theta[7:8])))

time1 <- proc.time()
cat(paste("One-Stepper Time:",(time1-time0)[3],"\n"))

out <- list(MODEL=list(MIX=mix, HPS=hpObj), POSTPROBS=Z)
return(out)
}

ebCoexpressFullTCAECM <- function(D, conditions, pattern, hpEsts,
                        controlOptions=list())
{

bool <- (sort(unique(conditions)) == sort(unique(unlist(pattern@patterns))))
if(sum(!bool) != 0)
  stop("Conditions and pattern are inconsistent!")
controlList <- list(convtol=0.0005, enforceFloor=TRUE, applyTransform=TRUE,
                    verbose=1, subsize=NULL, m2MaxIter=100)
cOnames <- names(controlOptions)
if(!is.null(cOnames))
{
  if("verbose" %in% cOnames)
    controlList$verbose <- controlOptions$verbose
  if("convtol" %in% cOnames)
    controlList$convtol <- controlOptions$convtol
  if("enforceFloor" %in% cOnames)
    controlList$enforceFloor <- controlOptions$enforceFloor
  if("applyTransform" %in% cOnames)
    controlList$applyTransform <- controlOptions$applyTransform
  if("subsize" %in% cOnames)
    controlList$subsize <- controlOptions$subsize
  if("m2MaxIter" %in% cOnames)
    controlList$m2MaxIter <- controlOptions$m2MaxIter
}

n <- length(conditions)         # No. of samples
k <- length(unique(conditions)) # No. of conditions
unns <- array(NaN, k)
for(j in 1:k)
  unns[j] <- sum(conditions==j) # No. of samples per condition
L <- length(pattern@patterns)   # No. of patterns

fZ <- function(rho) { return(0.5*log((1+rho)/(1-rho))) }
p <- dim(D)[1]
if(dim(D)[2] != k)
  stop("2nd dimension of D does not match number of unique conditions!")
if(controlList$applyTransform)
  D <- fZ(D)

convtol <- controlList$convtol
enforceFloor <- controlList$enforceFloor
verbose <- controlList$verbose
m2MaxIter <- controlList$m2MaxIter
subp <- controlList$subsize
subFlag <- FALSE
if(!is.null(subp))
  subFlag <- TRUE
if(subFlag)
  if(subp > p)
    stop("The provided subsize is greater than the number of pairs!")

### End of Setup ###

# Hyperparam setup

G <- hpEsts[[1]]     # It is critically vital that this be visible
eMus <- hpEsts[[2]]
eTaus <- hpEsts[[3]]
eWs <- hpEsts[[4]]
eWs <- eWs/sum(eWs)

theta <- NULL
if(G==1)
  { theta <- c(eMus, eTaus) ; tl <- 2 } # Length 2
if(G==2)
  { theta <- c(eMus, eTaus, eWs[1]) ; tl <- 5 } # Length 5
if(G==3)
  { theta <- c(eMus, eTaus, eWs[1:2]) ; tl <- 8 } # Length 8
if(!(G %in% c(1,2,3)))
  stop("Invalid hpEsts object!")

mix <- c(0.9,rep((1-0.9)/(L-1), L-1)) # Realistic mix
Z <- matrix(0,p,L) # Z is #pairs by #patterns
Zra <- NULL
randa <- NULL
if(subFlag)
  randa <- sort(sample(1:p, subp, replace=FALSE))

# EM functions

f0toC <- function(zvals, sizes, theta)
{
  return(as.numeric(.Call("f0worker",zvalsx=zvals,sizes3x=(sizes-3),
                          thetax=theta, dx=length(zvals), PACKAGE="EBcoexpress")))
}

theMiniLik <- function(theta)
{
  if(G==3)
    if(theta[7]+theta[8] > 1)   # This is the constraint for G=3
      return(Inf)
  tot <- 0
  # Note: randa (length=subp), a subset of 1:p, must be visible,
  #  as must Zra <- Z[randa,]
  pred.dens2 <- Zra # Something subp-by-L
  for(ind in 1:length(randa))
    for(pat in 1:L)
      pred.dens2[ind,pat] <- theGManager(randa[ind], theta, pat)
  tot <- sum(Zra * pred.dens2)
  return(-tot)
}

theGLik <- function(theta)
{
  if(G==3)
    if(theta[7]+theta[8] > 1)   # This is the constraint for G=3
      return(Inf)
  tot <- 0
  pred.dens <- Z # Something p-by-L
  for(ind in 1:p)
    for(pat in 1:L)
      pred.dens[ind,pat] <- theGManager(ind, theta, pat)
  tot <- sum(Z * pred.dens)
  return(-tot)
}

theGManager <- function(ind, theta, pat)
{
  if(G==1) return(theGWorker1(ind, theta, pat))
  if(G==2) return(theGWorker2(ind, theta, pat))
  if(G==3) return(theGWorker3(ind, theta, pat))
  stop("Unreachable line reached in G Manager!")
}

theGWorker1 <- function(ind, theta, pat)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      zs <- D[ind,subInd]
      szs <- unns[subInd]
      df0 <- df0 + f0toC(zvals=zs,sizes=szs,theta)
    }
  return(df0)
}

theGWorker2 <- function(ind, theta, pat)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      w2 <- (1-theta[5])
      zs <- D[ind,subInd]
      szs <- unns[subInd]
      l1 <- theta[5]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[1],theta[3])))
      l2 <- w2*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[2],theta[4])))
      df0 <- df0 + log(l1+l2)
    }
  return(df0)
}

theGWorker3 <- function(ind, theta, pat)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      w3 <- (1-theta[7]-theta[8])
      zs <- D[ind,subInd]
      szs <- unns[subInd]
      l1 <- theta[7]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[1],theta[4])))
      l2 <- theta[8]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[2],theta[5])))
      l3 <- w3*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[3],theta[6])))
      df0 <- df0 + log(l1+l2+l3)
    }
  return(df0)
}

# Other setup for the run:

mtol  <- convtol
done <- FALSE
iter <- 1
citer <- 1
omix <- rep(-1,length(mix))
cmix <- omix
ctheta <- rep(-10, length(theta))
cdone <- FALSE

##########################################
# Here is the Plan: [E M2 [E M1]] RET(Z) #
##########################################

time0 <- proc.time()
while(!cdone)
{
done <- FALSE
if(verbose > 0)
  { cat(paste("TCA-ECM Iteration:",citer,"\n")) ; cat("{\n") }

#@# Step 1: E-Step

if(verbose > 0)
  cat(" Begin Phase I (Initial E-Step) ...\n")

for(j in 1:L)
  Z[,j] <- exp(log(mix[j]) +
               sapply(1:p, FUN=theGManager, theta=theta, pat=j))
Z <- Z/(apply(Z,1,FUN=sum))
if(subFlag)
  Zra <- Z[randa,]

#@# Step 2: M2-Step

if(verbose > 0)
  cat(" Begin Phase II (M2-Step) ...\n")

lowlim <- NULL
uplim <- NULL
if(G==1)
  {
    lowlim <- c(-Inf, 0)
    uplim <- c(Inf, Inf)
  }
if(G==2)
  {
    lowlim <- c(-Inf, -Inf, 0, 0, 0)
    uplim <- c(Inf, Inf, Inf, Inf, 1)
  }
if(G==3)
  {
    lowlim <- c(-Inf, -Inf, -Inf, 0, 0, 0, 0, 0)
    uplim <- c(Inf, Inf, Inf, Inf, Inf, Inf, 1, 1)
  }

con <- list()
if(!is.null(m2MaxIter))
  con <- list(maxfun=m2MaxIter)

if(subFlag)
{
  resu <- suppressWarnings(bobyqa(theta, theMiniLik, lower=lowlim, upper=uplim, control=con))
  theta <- resu$par
}
if(!subFlag)
{
  resu <- suppressWarnings(bobyqa(theta, theGLik, lower=lowlim, upper=uplim, control=con))
  theta <- resu$par
}

#@# Step 3: [E M1] Cycle

if(verbose > 0)
  cat(" Begin Phase III ([E M1] Cycle) ...\n")

ZV <- matrix(0, p, L)
for(j in 1:L)
  ZV[,j] <- sapply(1:p, FUN=theGManager, theta=theta, pat=j)

iter <- 1

while(!done)
{
  if(verbose > 0)
    cat(paste("  Iteration:",iter,"\n"))

  for(j in 1:L)
    Z[,j] <- exp(log(mix[j]) + ZV[,j])
  Z <- Z/(apply(Z,1,FUN=sum))
  p.temp <- apply(Z,2,sum)
  mix <- p.temp/sum(p.temp)

  # Checking for convergence

  mis <- array(NaN, L)
  for(j in 1:L)
  {
    mis[j] <- abs(mix[j] - omix[j])
    if(!is.finite(mis[j])) mis[j] <- 0
  }

  done <- TRUE
  for(j in 1:L)
    if(mis[j] >= convtol)
      done <- FALSE
  iter <- iter +1
  if(enforceFloor)
  {
    if(mix[1] < 0.8)
    {
      if(verbose>1)
        cat("  Floor enforced ...\n")
      if(L==2) mix <- c(0.8, 0.2)
      if(L>2)
      {
        ix <- mix[-1]
        ix <- ix*(0.2/sum(ix))
        mix <- c(0.8, ix)
      }
      for(j in 1:L)
        Z[,j] <- exp(log(mix[j]) + ZV[,j])
      Z <- Z/(apply(Z,1,FUN=sum))
      done <- TRUE
    }
  }
  omix <- mix
}

# If we have gotten to here, then one TCA-ECM iteration has finished

if(verbose > 0)
  cat("}\n")

mis <- array(NaN, L)
for(j in 1:L)
{
  mis[j] <- abs(mix[j] - cmix[j])
  if(!is.finite(mis[j])) mis[j] <- 0
}
tis <- array(NaN, L)
for(j in 1:(length(theta)))
{
  tis[j] <- abs(theta[j] - ctheta[j])
  if(!is.finite(tis[j])) tis[j] <- 0
}

  cdone <- TRUE
  for(j in 1:L)
    if(mis[j] >= convtol)
      cdone <- FALSE
  for(j in 1:(length(theta)))
    if(tis[j] >= convtol)
      cdone <- FALSE

  citer <- citer +1
  cmix <- mix
  ctheta <- theta

} # End of the cdone while loop

dimnames(Z)[[1]] <- dimnames(D)[[1]]
dimnames(Z)[[2]] <- c("EC",paste("DC",(1:(L-1)),sep=""))

hpObj <- NULL
if(G==1)
  hpObj <- list(G=G, MUS=theta[1], TAUS=theta[2], WEIGHTS=1)
if(G==2)
  hpObj <- list(G=G, MUS=theta[1:2], TAUS=theta[3:4],
                  WEIGHTS=c(theta[5], 1-theta[5]))
if(G==3)
  hpObj <- list(G=G, MUS=theta[1:3], TAUS=theta[4:6],
                  WEIGHTS=c(theta[7:8], 1-sum(theta[7:8])))

time1 <- proc.time()
cat(paste("TCA-ECM Time:",(time1-time0)[3],"\n"))

out <- list(MODEL=list(MIX=mix, HPS=hpObj), POSTPROBS=Z)
return(out)
}

ebCoexpressMeta <- function(DList, conditionsList, pattern, hpEstsList,
                        controlOptions=list())
{

# Set-up

S <- length(DList)
S2 <- length(conditionsList)
if(S != S2)
  stop("DList and conditionsList must be the same length!")

controlList <- list(convtol=0.0005, enforceFloor=TRUE, applyTransform=TRUE,
                    verbose=1)
cOnames <- names(controlOptions)
if(!is.null(cOnames))
{
  if("verbose" %in% cOnames)
    controlList$verbose <- controlOptions$verbose
  if("convtol" %in% cOnames)
    controlList$convtol <- controlOptions$convtol
  if("enforceFloor" %in% cOnames)
    controlList$enforceFloor <- controlOptions$enforceFloor
  if("applyTransform" %in% cOnames)
    controlList$applyTransform <- controlOptions$applyTransform
}

ns <- array(NaN, S)
for(s in 1:S)
  ns[s] <- length(conditionsList[[s]])      # No. of samples over studies
k <- dim(DList[[1]])[2]                     # No. of conditions
unns <- matrix(NaN, k, S)
for(s in 1:S)
  for(j in 1:k)
    unns[j,s] <- sum(conditionsList[[s]]==j) # No. of samp per cond per study
L <- length(pattern@patterns)               # No. of patterns
p <- dim(DList[[1]])[1]                     # No. of pairs

fZ <- function(rho) { return(0.5*log((1+rho)/(1-rho))) }

if(controlList$applyTransform)
{
  for(s in 1:S)
    DList[[s]] <- fZ(DList[[s]])
}

convtol <- controlList$convtol
enforceFloor <- controlList$enforceFloor
verbose <- controlList$verbose

### End of Setup ###

# Hyperparam setup

G <- vector("list",S)     # Must be visible
eMus <- vector("list",S)
eTaus <- vector("list",S)
eWs <- vector("list",S)

for(s in 1:S)
{
G[[s]] <- hpEstsList[[s]][[1]]
eMus[[s]] <- hpEstsList[[s]][[2]]
eTaus[[s]] <- hpEstsList[[s]][[3]]
temp <- hpEstsList[[s]][[4]]
eWs[[s]] <- temp/sum(temp)
}

theta <- vector("list", S)

for(s in 1:S)
{
  if(G[[s]]==1)
    theta[[s]] <- c(eMus[[s]], eTaus[[s]])
  if(G[[s]]==2)
    theta[[s]] <- c(eMus[[s]], eTaus[[s]], eWs[[s]][1])
  if(G[[s]]==3)
    theta[[s]] <- c(eMus[[s]], eTaus[[s]], eWs[[s]][1:2])
}

mix <- c(0.9,rep((1-0.9)/(L-1), L-1)) # Realistic mix
Z <- matrix(0,p,L) # Z is #pairs by #patterns

# EM functions

f0toC <- function(zvals, sizes, theta)
{
  return(as.numeric(.Call("f0worker",zvalsx=zvals,sizes3x=(sizes-3),
                          thetax=theta, dx=length(zvals), PACKAGE="EBcoexpress")))
}

theMetaManager <- function(ind, theta, pat) # Just taking prods over s
{
  tot <- 0
  for(s in 1:S)
    tot <- tot + theGManager2(ind, theta, pat, s)
  return(tot)
}

theGManager2 <- function(ind, theta, pat, s)
{
  if(G[[s]]==1) return(theGWorker1m(ind, theta[[s]], pat, s))
  if(G[[s]]==2) return(theGWorker2m(ind, theta[[s]], pat, s))
  if(G[[s]]==3) return(theGWorker3m(ind, theta[[s]], pat, s))
  stop("Unreachable line reached in G Manager!")
}

theGWorker1m <- function(ind, theta, pat, s)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      zs <- DList[[s]][ind,subInd]
      szs <- unns[subInd,s]
      df0 <- df0 + f0toC(zvals=zs,sizes=szs,theta)
    }
  return(df0)
}

theGWorker2m <- function(ind, theta, pat, s)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      w2 <- (1-theta[5])
      zs <- DList[[s]][ind,subInd]
      szs <- unns[subInd,s]
      l1 <- theta[5]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[1],theta[3])))
      l2 <- w2*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[2],theta[4])))
      df0 <- df0 + log(l1+l2)
    }
  return(df0)
}

theGWorker3m <- function(ind, theta, pat, s)
{
  subK <- length(pattern@patterns[[pat]])
  df0 <- 0
  for(i in 1:subK)
    {
      subInd <- pattern@patterns[[pat]][[i]]
      w3 <- (1-theta[7]-theta[8])
      zs <- DList[[s]][ind,subInd]
      szs <- unns[subInd,s]
      l1 <- theta[7]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[1],theta[4])))
      l2 <- theta[8]*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[2],theta[5])))
      l3 <- w3*exp(f0toC(zvals=zs,sizes=szs,theta=c(theta[3],theta[6])))
      df0 <- df0 + log(l1+l2+l3)
    }
  return(df0)
}

# Other setup for the run:

done <- FALSE
iter <- 1
omix <- rep(-1,length(mix))

###################################
# Here is the Plan: [E M1] RET(Z) #
###################################

time0 <- proc.time()

if(verbose > 0)
  cat("Running the [E M1] Cycle ...\n")

ZV <- matrix(0, p, L)
for(j in 1:L)
  ZV[,j] <- sapply(1:p, FUN=theMetaManager, theta=theta, pat=j)

while(!done)
{
  if(verbose > 0)
    cat(paste("Iteration:",iter,"\n"))

  for(j in 1:L)
    Z[,j] <- exp(log(mix[j]) + ZV[,j])
  Z <- Z/(apply(Z,1,FUN=sum))
  p.temp <- apply(Z,2,sum)
  mix <- p.temp/sum(p.temp)

  # Checking for convergence

  mis <- array(NaN, L)
  for(j in 1:L)
  {
    mis[j] <- abs(mix[j] - omix[j])
    if(!is.finite(mis[j])) mis[j] <- 0
  }

  done <- TRUE
  for(j in 1:L)
    if(mis[j] >= convtol)
      done <- FALSE
  iter <- iter +1
  if(enforceFloor)
  {
    if(mix[1] < 0.8)
    {
      if(verbose>1)
        cat("Floor enforced ...\n")
      if(L==2) mix <- c(0.8, 0.2)
      if(L>2)
      {
        ix <- mix[-1]
        ix <- ix*(0.2/sum(ix))
        mix <- c(0.8, ix)
      }
      for(j in 1:L)
        Z[,j] <- exp(log(mix[j]) + ZV[,j])
      Z <- Z/(apply(Z,1,FUN=sum))
      done <- TRUE
    }
  }
  omix <- mix
}

dimnames(Z)[[1]] <- dimnames(DList[[1]])[[1]]
dimnames(Z)[[2]] <- c("EC",paste("DC",(1:(L-1)),sep=""))

time1 <- proc.time()
cat(paste("Meta Analysis Time:",(time1-time0)[3],"\n"))

out <- list(MODEL=list(MIX=mix, HPS=hpEstsList), POSTPROBS=Z)
return(out)
}

showPair <- function(pair, X, conditions, gsep="~", regLine=TRUE,
                      useBWMC=TRUE, colors=NULL, ...)
{
  dots <- list(...)
  args <- dots[names(dots) %in% c("lwd","lty")]
  args2 <- dots[!(names(dots) %in% c("lwd","lty"))]
  if(is.null(args$lwd))
    args$lwd <- 2
  if(is.null(args$lty))
    args$lty <- 1

  if(is.null(colors))
    colors <- palette()
  pa1 <- strsplit(pair, gsep)[[1]][1]
  pa2 <- strsplit(pair, gsep)[[1]][2]
  if(is.null(args2$xlab))
    args2$xlab <- pa1
  if(is.null(args2$ylab))
    args2$ylab <- pa2

  condM <- conditions
  K <- max(condM)
  monA <- X[pa1,]
  monB <- X[pa2,]

  getUsed <- function(aX, aY)
  {
    n <- length(aX)
    UXi <- (aX-median(aX))/(9*mad(aX, const=1))
    A <- UXi >= -1 & UXi <= 1
    VYi <- (aY-median(aY))/(9*mad(aY, const=1))
    B <- VYi >= -1 & VYi <= 1
    return(A & B)
  }

  args3 <- append(list(x=monA, y=monB, type="p", col=colors[condM]), args2)
  do.call(plot, args3)
  box()
  if(regLine)
  {
    if(useBWMC)
    {
      for(k in 1:K)
      {
        monUsedN <- getUsed(monA[condM==k], monB[condM==k])
        lm1 <- lm((monB[condM==k])[monUsedN] ~ (monA[condM==k])[monUsedN])
        co <- lm1$coef
        yint <- co[1] ; slo <- co[2]
        x0 <- min((monA[condM==k])[monUsedN])
        x1 <- max((monA[condM==k])[monUsedN])
        y0 <- x0*slo + yint
        y1 <- x1*slo + yint
        segargs <- append(list(x0=x0, y0=y0, x1=x1, y1=y1, col=colors[k]),args)
        do.call(segments, segargs)
      }
    }
    if(!useBWMC)
    {
      for(k in 1:K)
      {
        lm1 <- lm(monB[condM==k] ~ monA[condM==k])
        co <- lm1$coef
        yint <- co[1] ; slo <- co[2]
        x0 <- min((monA[condM==k]))
        x1 <- max((monA[condM==k]))
        y0 <- x0*slo + yint
        y1 <- x1*slo + yint
        segargs <- append(list(x0=x0, y0=y0, x1=x1, y1=y1, col=colors[k]),args)
        do.call(segments, segargs)
      }
    }
  }
  return(invisible(NULL))
}

showNetwork <- function(geneSet, D, condFocus,
  gsep="~", layout="kamada.kawai", seed=NULL, hidingThreshold=NULL, ...)
{
  library(graph)
  library(igraph)
  library(colorspace)

  k <- condFocus
  temp <- geneSet
  len <- length(temp)
  allPairNames <- array("", len*len)
  ca <- 0
  for(i in 1:len)
    for(j in 1:len)
    {
      ca <- ca + 1
      targ <- paste(temp[i],temp[j],sep=gsep)
      allPairNames[ca] <- targ
    }
  length(allPairNames)
  valid <- allPairNames %in% dimnames(D)[[1]]
  validPairNames <- allPairNames[valid]
  
  mat1 <- matrix(0, len, len)
  for(i in 1:len)
  for(j in 1:len)
  {
    targ <- paste(temp[i],temp[j],sep=gsep)
    if(targ %in% validPairNames)
      {
        val1 <- D[targ,k]
        mat1[i,j] <- val1
        mat1[j,i] <- val1
      }
  }
  dimnames(mat1) <- list(temp,temp)

  dummy <- graph.full(len, directed=FALSE)

  ecol7a <- array("",choose(len,2))

  ca <- 0
  for(i in 1:(len-1))
    for(j in (i+1):len)
    {
      ca <- ca + 1
      f1 <- temp[i]
      f2 <- temp[j]
      valA <- mat1[f1,f2]*2
      if(valA > 1) valA <- 1
      if(valA < -1) valA <- -1 # I have reversed the Red/Blue below

      ecol7a[ca] <- hex(RGB(1,1-abs(valA),1-abs(valA)))
      if(valA < 0) ecol7a[ca] <- hex(RGB(1-abs(valA),1-abs(valA),1))
      if(!is.null(hidingThreshold))
        if(valA < hidingThreshold && valA > -1*hidingThreshold)
          ecol7a[ca] <- "transparent"
    }

  if(!is.null(seed)) set.seed(seed)
  layo <- layout.circle(dummy) # Decoy, should be evaluated in the next line
  laytext <- paste("layo <- layout.",layout,"(dummy)",sep="")
  eval(parse(text=laytext))

  plot(dummy, layout=layo, vertex.label=temp, edge.arrow.mode=0,
    edge.color=ecol7a, ...)
  
  return(invisible(NULL))
}

rankMyGenes <- function(emOut, thresh=0.95, sep="~")
{
  pps <- emOut$POSTPROBS
  ethresh <- 1 - thresh
  dcpairNames <- (dimnames(pps)[[1]])[pps[,1] < ethresh]
  allNames <- unlist(strsplit(dcpairNames, split=sep))
  dcTable <- sort(table(allNames),decr=TRUE)
  return(dcTable)
}

# EOF