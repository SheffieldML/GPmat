\name{filterGenes}
\Rdversion{1.0}
\alias{filterGenes}
\alias{averageToSDRatio}
\alias{writeRatioData}
\title{Determining the ratio of the average and the standard deviation}
\description{
  'filterGenes' determines which genes have a ratio of at least the specified level.

  'averageToSDRatio' determines the ratio of the average and the standard deviation for each time point of each gene in the preprocessed data. The mean and median of the ratios are calculated for each gene.

  'writeRatioData' writes the means and medians of the ratios of genes to a text file.
}
\usage{
approvedGenes <- filterGenes(ratioData, filterLimit = 1.5, useMedians = TRUE)
ratioData <- averageToSDRatio(preprocData, searchedGenes = "100001_at", search = FALSE)
writeRatioData(ratioData, genes, fileName)
}
\arguments{
  \item{ratioData}{The data containing the ratios, means, medians and the relevant genes. It's used in 'filterGenes' for filtering the genes and in 'writeRatioData' for writing the data into a text file.}
  \item{filterLimit}{Genes with a filtering value above this figure are accepted after filtering. If this value is 0, all genes will be accepted.}
  \item{useMedians}{A logical value determining whether medians are used for filtering the genes; if false, means are used.}
  \item{preprocData}{The preprocessed data to be used.}
  \item{searchedGenes}{An array of the genes to be searched if the variable search is true.}
  \item{search}{A logical value determining whether only the ratios related to specific genes are calculated.}
  \item{genes}{The genes of the data given to the function writeRatioData.}
  \item{fileName}{The name of the text file in which the data is written.}
}
\value{
  'filterGenes' returns an array of the genes that passed the filtering.

  'averageToSDRatio' returns data containing an array of ratios, means and medians for each gene.
}
\seealso{
  \code{\link{GPrank}, \link{processData}}.
}
\examples{
## missing
}
\keyword{model}
