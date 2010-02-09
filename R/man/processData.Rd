\name{processData}
\Rdversion{1.0}
\alias{processData}
\title{Processing further data preprocessed by mmgmos}
\description{
  'processData' processes further data preprocessed by mmgmos in 'preprocess_affy'. The resulting data can then be used as input for the function 'GPrank'.
}
\usage{
preprocData <- processData(data, replicates = NULL, searchedGenes = NULL)
}
\arguments{
  \item{data}{The preprocessed data from mmgMOS to be used.}
  \item{replicates}{The replicate structure of the data indicating
    which expression data points belong to which replicates.
    By default all the data points
    are assumed to be from same replicate.}
  \item{searchedGenes}{If set, an array of the genes to be include in
    the result.}
}
\details{
  The expression data and percentiles are normalized. A normal
  distribution is then fitted into the data with distfit. The expression
  data is then rescaled so that average standard deviation of curves is
  1. The ratio of the average to the standard deviation of each value in
  the expression data is also included, as well as a mean and median for
  each gene, to be used for filtering genes.
}
\value{
  'processData' returns the new processed data containing the expression data, variances, time points, genes and scale values.
}
\author{
  Jonatan Ropponen, Antti Honkela
}
\seealso{
  \code{\link{GPLearn}, \link{GPrankTargets}}.
}
\examples{
## missing
}
\keyword{model}
