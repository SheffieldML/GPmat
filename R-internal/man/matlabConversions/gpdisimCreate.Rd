\name{gpdisimCreate}
\Rdversion{1.0}
\alias{gpdisimCreate}
\title{Create a GPDISIM model.}
\description{
  creates a model for single input motifs with Gaussian processes.
}
\usage{
model <- gpdisimCreate(numGenes, numProteins, times, geneVals,
                       geneVars, options, annotation)
}
\arguments{
  \item{numGenes}{number of genes to be modelled in the system.}
  \item{numProteins}{number of proteins to be modelled in the system.}
  \item{times}{the time points where the data is to be modelled.}
  \item{geneVals}{the values of each gene at the different time points.}
  \item{geneVars}{the varuabces of each gene at the different time points.}
  \item{options}{options structure, the default options can be generated using gpsimOptions.}
  \item{annotation}{annotation for the data (gene names, etc.) that is stored with the model. (Optional)}
}
\value{
  \item{model}{model structure containing default parameterisation.}
}
\seealso{
\code{\link{modelCreate}, \link{gpsimOptions}}.
}
\examples{
## missing
}
\keyword{model}
