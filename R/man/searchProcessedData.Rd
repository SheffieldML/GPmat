\name{searchProcessedData}
\Rdversion{1.0}
\alias{searchProcessedData}
\title{Searching preprocessed data for the data of specific genes}
\description{
  Searches for the data of specific genes from given preprocessed data and returns a new set of data containing only the desired data.
}
\usage{
newData <- searchProcessedData(preprocData, searchedGenes)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{searchedGenes}{An array of the genes to be searched.}
}
\value{
  Returns a new set of preprocessed data containing only the data of the specified genes.
}
\seealso{
  \code{\link{GPrank}}.
}
\examples{
## missing
}
\keyword{model}
