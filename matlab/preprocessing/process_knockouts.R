DATAPATH <- "~/data/eileen_data/mutant_data/"

FILES <- c("Mef2-loss-cDNA.filtered.txt", "bin-loss-INDAC.filtered.txt", "tin-loss-cDNA.filtered.txt", "twi-loss-INDAC.filtered.txt")

TFS <- c("mef2", "bin", "tin", "twi")

datas <- list()

for (i in 1:4) {
  datas[[i]] <- read.table(paste(DATAPATH, FILES[i], sep=""), sep="\t", quote="", comment.char="", header=TRUE)
  inds <- grep("fdr", names(datas[[i]]))
  inds2 <- grep("meanRatio", names(datas[[i]]))
  qvals <- datas[[i]][,inds]
  qvals[is.na(qvals)] <- 1
  ratios <- datas[[i]][,inds2]
  ratios[is.na(ratios)] <- 0
  qvals[abs(ratios) < .5] <- 1
  dims <- dim(qvals)
  minave <- apply((qvals[,-dims[2]] + qvals[,-1])/2, 1, min)
  names(minave) <- datas[[i]][,"accession"]
  q <- sort(minave)
  write.table(q, paste(DATAPATH, TFS[i], "_knockout_qvals.txt", sep=""), row.names=names(q), quote=FALSE, sep="\t")
#write.table()
}
