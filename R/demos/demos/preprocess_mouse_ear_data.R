library(puma)
library(affy)
expfiles <- c(paste("RA_exp_C", c(0, 1, 2, 4, 7, 14), ".CEL", sep=""),
              paste("RA_exp_RA", c(0, 1, 2, 4, 7, 14), ".CEL", sep=""))
expdata <- ReadAffy(filenames=expfiles, celfile.path="/Users/ahonkela/data/marta_data/orig")
pData(expdata) <- data.frame("time.days" = c(0, 1, 2, 4, 7, 14, 0, 1, 2, 4, 7, 14), "RA"=c(0,0,0,0,0,0,1,1,1,1,1,1), row.names=rownames(pData(expdata)))
mmgmos_exprs <- mmgmos(expdata)
