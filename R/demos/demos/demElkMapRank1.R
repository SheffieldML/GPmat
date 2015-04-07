library(gpsim)

options(error = recover)

fileName <- paste("demElkMapLinear1.RData", sep="/")
if (file.exists(fileName)) {
  load(fileName)
} else {
  source("demElkMapLinear1.R")
}

expType <- "demElkMapRank"
expNo <- 1

loadopt <- 4
rankingData <- gpsimLoadElkData(loadopt)
ry <- rankingData$y
ryvar <- rankingData$yvar
rgenes <- rankingData$genes
rtimes <- rankingData$times
rscale <- rankingData$scale

rm(rankingData)

Nrep <- length(model$comp)

## For testing code
affInd <- rgenes[1:5]
Ntest <- length(affInd)

matchInd <- c()
# options(warn=1)
for ( i in 1:Ntest ) {
  ind <- grep(affInd[i], rgenes)
  index <- c()
  for ( j in seq(length=length(ind)) ) {
    if ( affInd[i] %in% rgenes[ind[j]] )
      index <- c(index, ind[j])      
  }
  if ( length(index) != 1 ) {
    warnmsg <- paste("Too many or too few matches for ", affInd[i], "!", sep="")
    warning(warnmsg)
  } else {
    matchInd <- c(matchInd, index)
  }
}
# options(warn=-1)

# testInd <- 4
# matchInd <- matchInd[testInd]

for ( i in seq(length=Nrep) ) {
  ry[[i]] <- as.matrix(ry[[i]][,matchInd])
  ryvar[[i]] <- as.matrix(ryvar[[i]][,matchInd])
}
rgenes <- rgenes[matchInd]
rscale <- rscale[matchInd]

ranks <- rgpsimMapRank(model, rtimes, ry, ryvar, rgenes)
rankOut <- list(genes=c(), llscore=c(), normalisedLl=c(), normLl2=c())
Ntest <- length(ranks$comp)
for ( i in seq(length=Ntest) ) {
  rankOut$genes <- c(rankOut$genes, ranks$comp[[i]]$geneName)
  rankOut$llscore <- c(rankOut$llscore, ranks$comp[[i]]$llscore)
  rankOut$normlisedLl <- c(rankOut$normlisedLl, ranks$comp[[i]]$normlisedLl)
  normLl2 <- -(ranks$comp[[i]]$llscore-ranks$comp[[i]]$trainingLl-ranks$comp[[i]]$flatLl)
  rankOut$normLl2 <- c(rankOut$normLl2, normLl2)
}

fileName <- paste(expType, expNo, ".Rdata", sep="")
save(model, ranks, expType, expNo, file=fileName)

rankout <- data.frame(rankOut$genes, rankOut$llscore, rankOut$normLl2)

outFileName <- paste(dataPath, "demBarencoRank1.csv", sep="/")
write.csv(rankout, file=outFileName)

# Show the profiles of interested individual genes
showGenes <- c("216396_s_at", "218627_at")

rgpsimRankSingleProfile(ranks, showGenes, expType, expNo)


