library(EBcoexpress)

data(fiftyGenes)
tinyCond <- c(rep(1,100),rep(2,25))
tinyPat <- ebPatterns(c("1,1","1,2"))
D <- makeMyD(fiftyGenes, tinyCond, useBWMC=TRUE)
set.seed(3)
initHP <- initializeHP(D, tinyCond)

zout <- ebCoexpressZeroStep(D, tinyCond, tinyPat, initHP)
par(mfrow=c(2,1))
priorDiagnostic(D, tinyCond, zout, 1)
priorDiagnostic(D, tinyCond, zout, 2)
par(mfrow=c(1,1))

oout <- ebCoexpressOneStep(D, tinyCond, tinyPat, initHP)
par(mfrow=c(2,1))
priorDiagnostic(D, tinyCond, oout, 1)
priorDiagnostic(D, tinyCond, oout, 2)
par(mfrow=c(1,1))

fout <- ebCoexpressFullTCAECM(D, tinyCond, tinyPat, initHP)
par(mfrow=c(2,1))
priorDiagnostic(D, tinyCond, fout, 1)
priorDiagnostic(D, tinyCond, fout, 2)
par(mfrow=c(1,1))

ii <- c(rep(1:24, times=24:1), rep(26:49, times=24:1))
base1 <- 2:25
for(j in 3:25) base1 <- c(base1, j:25)
base2 <- 27:50
for(j in 28:50) base2 <- c(base2, j:50)
jj <- c(base1, base2)
TP <- mapply("X",ii,"~","X",jj,FUN=paste,sep="")
names(TP) <- NULL
numDC <- length(TP)

result0 <- zout$POSTPROBS
result1 <- oout$POSTPROBS
resultF <- fout$POSTPROBS

###

thresh0 <- crit.fun(result0[,1], 0.05)
ppbDC0 <- 1-result0[,1]
crit_s <- crit.fun(result0[,1], 0.05)
kept_s <- ppbDC0[ppbDC0 >= crit_s]
kept_h <- ppbDC0[ppbDC0 >= 0.95]
nk_s <- length(kept_s)
klabs_s <- names(kept_s)
nk_h <- length(kept_h)
klabs_h <- names(kept_h)
nY_s <- sum(klabs_s %in% TP)
nY_h <- sum(klabs_h %in% TP)
(nk_s - nY_s)/nk_s           # Soft threshold Obs. FDR
(nk_h - nY_h)/nk_h           # Hard threshold Obs. FDR
nY_s/numDC                   # Soft threshold Obs. Power
nY_h/numDC                   # Hard threshold Obs. Power

thresh1 <- crit.fun(result1[,1], 0.05)
ppbDC1 <- 1-result1[,1]
crit_s <- crit.fun(result1[,1], 0.05)
kept_s <- ppbDC1[ppbDC1 >= crit_s]
kept_h <- ppbDC1[ppbDC1 >= 0.95]
nk_s <- length(kept_s)
klabs_s <- names(kept_s)
nk_h <- length(kept_h)
klabs_h <- names(kept_h)
nY_s <- sum(klabs_s %in% TP)
nY_h <- sum(klabs_h %in% TP)
(nk_s - nY_s)/nk_s           # Soft threshold Obs. FDR
(nk_h - nY_h)/nk_h           # Hard threshold Obs. FDR
nY_s/numDC                   # Soft threshold Obs. Power
nY_h/numDC                   # Hard threshold Obs. Power

threshF <- crit.fun(resultF[,1], 0.05)
ppbDCF <- 1-resultF[,1]
crit_s <- crit.fun(resultF[,1], 0.05)
kept_s <- ppbDCF[ppbDCF >= crit_s]
kept_h <- ppbDCF[ppbDCF >= 0.95]
nk_s <- length(kept_s)
klabs_s <- names(kept_s)
nk_h <- length(kept_h)
klabs_h <- names(kept_h)
nY_s <- sum(klabs_s %in% TP)
nY_h <- sum(klabs_h %in% TP)
(nk_s - nY_s)/nk_s           # Soft threshold Obs. FDR
(nk_h - nY_h)/nk_h           # Hard threshold Obs. FDR
nY_s/numDC                   # Soft threshold Obs. Power
nY_h/numDC                   # Hard threshold Obs. Power

# Visualization

twentyGeneNames <- dimnames(fiftyGenes)[[1]][c(1:10,26:35)]

showNetwork(twentyGeneNames, D, condFocus = 1, gsep = "~",
  layout = "kamada.kawai", seed = 5, vertex.shape="circle",
  vertex.label.cex=1, vertex.color="white", edge.width=2,
  vertex.frame.color="black", vertex.size=20,
  vertex.label.color="black", vertex.label.family="sans")
#
showNetwork(twentyGeneNames, D, condFocus = 2, gsep = "~",
  layout = "kamada.kawai", seed = 5, vertex.shape="circle",
  vertex.label.cex=1, vertex.color="white", edge.width=2,
  vertex.frame.color="black", vertex.size=20,
  vertex.label.color="black", vertex.label.family="sans")
#

showNetwork(twentyGeneNames, D, condFocus = 1, gsep = "~",
  layout = "kamada.kawai", seed = 5, vertex.shape="circle",
  vertex.label.cex=1, vertex.color="white", edge.width=2,
  vertex.frame.color="black", vertex.size=20,
  vertex.label.color="black", vertex.label.family="sans",
  hidingThreshold=0.3)
#
showNetwork(twentyGeneNames, D, condFocus = 2, gsep = "~",
  layout = "kamada.kawai", seed = 5, vertex.shape="circle",
  vertex.label.cex=1, vertex.color="white", edge.width=2,
  vertex.frame.color="black", vertex.size=20,
  vertex.label.color="black", vertex.label.family="sans",
  hidingThreshold=0.3)
#

showNetwork(dimnames(fiftyGenes)[[1]], D, condFocus = 1, gsep = "~",
  layout = "circle", seed = 5, vertex.shape="circle",
  vertex.label.cex=1, vertex.color="white", edge.width=2,
  vertex.frame.color="black", vertex.size=20,
  vertex.label.color="black", vertex.label.family="sans",
  hidingThreshold=0.5)

showNetwork(dimnames(fiftyGenes)[[1]], D, condFocus = 2, gsep = "~",
  layout = "circle", seed = 5, vertex.shape="circle",
  vertex.label.cex=1, vertex.color="white", edge.width=2,
  vertex.frame.color="black", vertex.size=20,
  vertex.label.color="black", vertex.label.family="sans",
  hidingThreshold=0.5)

showPair("X1~X2",fiftyGenes, conditions=tinyCond, pch=20,
  xlim=c(-4,4), ylim=c(-4,4))
#

showPair("X26~X35",fiftyGenes, conditions=tinyCond, pch=20,
  xlim=c(-4,4), ylim=c(-4,4))
#

showPair("X1~X35",fiftyGenes, conditions=tinyCond, pch=20,
  xlim=c(-4,4), ylim=c(-4,4))
#

hubs <- rankMyGenes(oout)
print(hubs)

# End of Demo
