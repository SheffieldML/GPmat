\name{gpsimMapCreate}
\Rdversion{1.0}
\alias{gpsimMapCreate}
\title{Create a GPSIMMAP model.}
\description{
  creates a model for single input motifs with Gaussian processes.
}
\usage{
model <- gpsimMapCreate(numGenes, numProteins, times, geneVals,
                        geneVars, options)
}
\arguments{
  \item{numGenes}{number of genes to be modelled in the system.}
  \item{numProteins}{number of proteins to be modelled in the system.}
  \item{times}{the time points where the data is to be modelled.}
  \item{geneVals}{the values of each gene at the different time points.}
  \item{geneVars}{the varuabces of each gene at the different time points.}
  \item{options}{options structure, the default options can be generated using gpsimMapOptions.}
}
\value{
  \item{model}{model structure containing default parameterisation.}
}
\seealso{
\code{\link{modelCreate}, \link{gpsimMapOptions}}.
}
\examples{
## missing
}
\keyword{model}
