CELPATH <- paste(Sys.getenv("HOME"), "/data/dros_data/embryo_tc_array_data", sep="")

expfiles <- paste(rep(paste("embryo_tc", 2*2:4, sep="_"), each=12), "_", 1:12, ".CEL", sep="")
expdata <- ReadAffy(filenames=expfiles, celfile.path=CELPATH)
pData(expdata) <- data.frame("time.h" = rep(1:12, 3), row.names=rownames(pData(expdata)))
mmgmos_exprs <- mmgmos(expdata)
# write.reslts(mmgmos_exprs, file="mmgmos_exprs")

nameMapping <- "FLYBASE"
library(paste(annotation(mmgmos_exprs), '.db', sep=""), character.only=TRUE)
myMapping <- get(paste(annotation(mmgmos_exprs), nameMapping, sep=""))

genenames <- as.list(myMapping)

write.table(as.matrix(genenames), file='dros_fbgn_annotations.txt', quote=FALSE, col.names=FALSE, sep='\t')

#row.names(exprs(mmgmos_exprs)) <- genenames
#row.names(se.exprs(mmgmos_exprs)) <- genenames
#row.names(prcfive(mmgmos_exprs)) <- genenames
#row.names(prctwfive(mmgmos_exprs)) <- genenames
#row.names(prcfifty(mmgmos_exprs)) <- genenames
#row.names(prcsevfive(mmgmos_exprs)) <- genenames
#row.names(prcninfive(mmgmos_exprs)) <- genenames

write.reslts(mmgmos_exprs, file="mmgmos_exprs")

genenames <- as.list(drosgenome1SYMBOL)
write.table(as.matrix(genenames), file='dros_symbol_annotations.txt', quote=FALSE, col.names=FALSE, sep='\t')
