\name{processData}
\Rdversion{1.0}
\alias{processData}
\alias{searchData}
\title{Processing further data preprocessed by mmgmos}
\description{
  'processData' processes further data preprocessed by mmgmos in 'preprocess_affy'. The resulting data can then be used as input for the function 'GPrank'.

  'searchData' searches for the data of specified genes from the data preprocessed by mmgmos.
}
\usage{
preprocData <- processData(data, searchedGenes = "100001_at", search = FALSE)
newData <- searchData(yFull, genes, pcts, searchedGenes, amountOfColumns)
}
\arguments{
  \item{data}{The preprocessed data to be used.}
  \item{searchedGenes}{An array of the genes to be searched if the variable search is true.}
  \item{search}{A logical value determining whether the preprocessed data is searched for the data of the specified genes.}
  \item{yFull}{The expression data searched.}
  \item{genes}{An array of all genes.}
  \item{pcts}{An array of the percentiles at the 5 used time points.}
  \item{amountOfColumns}{The amount of columns in the expression data, in other words the amount of time points.}
}
\details{
  In 'processData', the expression data and percentiles are normalized. A normal distribution is then fitted into the data with distfit. The expression data is then rescaled so that average standard deviation of curves is 1.
}
\value{
  'processData' returns the new processed data containing the expression data, variances, time points, genes and scale values.

  'searchData' returns the data of the specified genes.
}
\seealso{
  \code{\link{GPrank}, \link{preprocess_affy}}.
}
\examples{
## missing
}
\keyword{model}
