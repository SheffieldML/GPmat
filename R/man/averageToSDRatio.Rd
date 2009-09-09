\name{averageToSDRatio}
\Rdversion{1.0}
\alias{averageToSDRatio}
\alias{writeRatioData}
\title{Determining the ratio of the average and the standard deviation}
\description{
  'averageToSDRatio' determines the ratio of the average and the standard deviation for each time point of each gene in the preprocessed data. The mean and median of the ratios are calculated for each gene.

  'writeRatioData' writes the means and medians of the ratios of genes to a text file.
}
\usage{
ratioData <- averageToSDRatio(preprocData, searchedGenes = "100001_at", search = FALSE)
writeRatioData(ratioData, genes, fileName)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{searchedGenes}{An array of the genes to be searched if the variable search is true.}
  \item{search}{A logical value determining whether only the ratios related to specific genes are calculated.}
  \item{ratioData}{The data to be written into a text file, given as output by the function averageToSDRatio.}
  \item{genes}{The genes of the data given to the function writeRatioData.}
  \item{fileName}{The name of the text file in which the data is written.}
}
\value{
  'averageToSDRatio' returns data containing a list of ratios, means and medians for each gene.
}
\seealso{
  \code{\link{processData}}.
}
\examples{
## missing
}
\keyword{model}
