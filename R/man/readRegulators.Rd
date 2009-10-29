\name{readRegulators}
\Rdversion{1.0}
\alias{readRegulators}
\alias{readPMIDRegulators}
\alias{readProbeRegulators}
\alias{removeDuplicates}
\title{Converting alias regulators in a text file into probes}
\description{
  'readRegulators' reads regulators in alias (common gene symbol identifiers) form in a text file and converts them into probe (manufacturer identifier) form. 

  'readPMIDRegulators' reads regulators in PMID (PubMed identifier) form in a text file and converts them into probe form.

  'readProbeRegulators' reads regulators in probe form from a text file.

  'removeDuplicates' takes an array of regulators and removes all duplicates and NAs from it.
}
\usage{
regulators <- readRegulators(fileName)
regulators <- readPMIDRegulators(fileName)
regulators <- readProbeRegulators(fileName)
regulators <- removeDuplicates(fileName)
}
\arguments{
  \item{fileName}{The name of the text file to be read.}
  \item{regulators}{A list of regulators for removeDuplicates.}
}
\value{
  Each function returns an array of the regulators in probe form.
}
\author{
  Jonatan Ropponen
}
\seealso{
\code{\link{geneOutputFunctions}, \link{mgu74av2ALIAS2PROBE}, \link{mgu74av2PMID2PROBE}}.
}
\examples{
## missing
}
\keyword{model}
