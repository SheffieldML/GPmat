library(puma)
library(affy)
expfiles <- paste("GSM38985", 1:6, ".CEL", sep="")
expdata <- ReadAffy(filenames=expfiles, celfile.path="/share/bayes/data/marta")
pData(expdata) <- data.frame("time.h" = c(0, 1, 2, 4, 7, 14), row.names=rownames(pData(expdata)))
# expdata@cdfName <- "dmgenome1dmrefseq"
mmgmos_exprs_refseq <- mmgmos(expdata)
# write.reslts(mmgmos_exprs_ug, file="mmgmos_exprs_ug")
