# MAKETP63DATASET Prepare the TP63 dataset.
# FORMAT
# DESC Loads the CEL files and pre-processes them with the PUMA package (RMA or
# MMGMOS normalization).
# Saves the final form of the dataset 
#
# SEEALSO : rma, mmgmos
#
# COPYRIGHT: Alfredo A. Kalaitzis, 2010, 2011
#
# DATASETS

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("puma")
## Load packages.
require(puma)

## Generate .CEL filenames
expfiles_tp63 = paste("GSM266", 780:792, ".CEL", sep="")
## Read CEL files and load into AffyBatch objects (your storage folder of .CEL files)
expdata_tp63 = ReadAffy(filenames=expfiles_tp63, celfile.path="/localhome/alkalait/mlprojects/datasets/R/GSE10562_RAW")
## Change annotation of samples according to their time of sampling
pData(expdata_tp63) = data.frame("time.h" =
  c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240),
  row.names=rownames(pData(expdata_tp63)))
## Convert AffyBatch objects to ExpressionSet (normalize)
## RMA
eset_tp63_RMA = rma(expdata_tp63)
## MMGMOS
# eset_tp63_mmgmos <- mmgmos(expdata_tp63) ## Takes a while.

## RMA normalization boxplots.
par(mfrow=c(1,2)) ## Sub-figures
## RMA
boxplot(data.frame(exprs(eset_tp63_RMA)),main="TP63 treatment")
## MMGMOS
# boxplot(data.frame(exprs(eset_tp63_mmgmos)),main="TP63 treatment")

## Column names (times).
## RMA
exprs_tp63_RMA = exprs(eset_tp63_RMA)
dimnames(exprs_tp63_RMA)[[2]] = seq(0,240,by=20)
## MMGMOS
# exprs_tp63_MMGMOS = exprs(eset_tp63_mmgmos)
# dimnames(exprs_tp63_MMGMOS)[[2]] = seq(0,240,by=20)

## Import Della Gatta data.
DGdata = read.csv("~/mlprojects/gprege/R/DellaGatta_SupTable1.csv", sep='')

## Create labels for the real dataset in Della Gatta et.al (2008).
## Assign '1' to each probeID in the dataset that is present in the TSNI ranking list.
DGatta_labels_byTSNI = dimnames(exprs_tp63_RMA)[[1]] %in% DGdata[,2]
DGatta_labels_byTSNItop100 = dimnames(exprs_tp63_RMA)[[1]] %in% DGdata[1:100,2]

## Extract gene symbols from gene-probe IDs. Takes a few minutes!
require(mouse430a2.db) ## Must have RSQLite package 0.7 or higher.
genesymbols = lapply(dimnames(exprs_tp63_RMA)[[1]], get, env=mouse430a2SYMBOL)
## Export data as a MATLAB .mat file.
# require(R.matlab)
# writeMat("TP63_aux.mat", exprs_tp63_RMA = exprs_tp63_RMA,
# #   exprs_tp63_MMGMOS = exprs_tp63_MMGMOS,
#   genesymbols = genesymbols,
#   DGlabels = DGlabels)

# install.packages('Rcompression', repos = "http://www.omegahat.org/R") # This is for reading a compressed MAT file.

BATSranking = matrix(0, dim(exprs_tp63_RMA)[1], 3)
for (i in 1:3) {
  tmp = read.table(url(paste('http://arxiv.org/src/1106.4333v1/anc/DGdat_p63_case',i,'_GL.txt',sep='')), skip=1) ## Read the gene numbers
  genenumbers = as.numeric(lapply( as.character(tmp[,2]), function(x) x=substr(x,2,nchar(x))))
  BATSranking[,i] <- tmp[sort(genenumbers, index.return=TRUE)$ix, 4] ## Sort rankings by gene numbers.
}

as.matrix(read.table(url('http://arxiv.org/src/1106.4333v1/anc/DGdat_p63_case1_GL.txt'), skip=1))

save(list = ls(all=TRUE), file = "~/mlprojects/datasets/R/DellaGattaData.RData")

## Gaussian process regression on the gene expression time-series.
# source('~/mlprojects/gprege/R/GPREGE.R')
# GPREGE(data=exprs_tp63_RMA, inputs=seq(0,240,by=20), explore=TRUE, labels=DGlabels)




