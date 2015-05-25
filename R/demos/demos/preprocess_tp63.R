library(puma)
library(affy)
expfiles <- paste("GSM2667", 80:92, ".CEL", sep="")
expdata <- ReadAffy(filenames=expfiles, celfile.path="/local/data/dibernardo")
pData(expdata) <- data.frame("time.h" = c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240), row.names=rownames(pData(expdata)))
# expdata@cdfName <- "dmgenome1dmrefseq"
mmgmos_exprs_refseq <- mmgmos(expdata)
# write.reslts(mmgmos_exprs_ug, file="mmgmos_exprs_ug")
