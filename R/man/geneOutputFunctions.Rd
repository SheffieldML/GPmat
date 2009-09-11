\name{geneOutputFunctions}
\Rdversion{1.0}
\alias{printScoreListGenes}
\alias{writeScoreListGenes}
\alias{writeGenes}
\title{Printing and writing genes and LLs}
\description{
  'printScoreListGenes' prints each gene and its log-likelihood in order of log-likelihood starting from the highest. 

  'writeScoreListGenes' writes each gene and its log-likelihood in order of log-likelihood starting from the highest into a text file.

  'writeGenes' writes the genes of a given array into a text file.
}
\usage{
printScoreListGenes(scoreList)
writeScoreListGenes(scoreList, fileName)
writeGenes(genes, fileName)
}
\arguments{
  \item{scoreList}{A score list containing the genes and related data. Usually received as output from the function GPrankTargets or GPrankTFs.}
  \item{fileName}{The name of the text file in which the genes and LLs are written.}
  \item{genes}{An array of genes to be written in a text file.}
}
\seealso{
  \code{\link{GPrankTargets}, \link{GPrankTFs}, \link{readRegulators}}.
}
\examples{
## missing
}
\keyword{model}
