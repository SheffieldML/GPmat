\name{gpDataIndices}
\Rdversion{1.0}
\alias{gpDataIndices}
\title{Return indices of present data.}
\description{
  returns the indices of data which is not missing for a given dimension in the GP-LVM and a block number in the PITC approximation.
}
\usage{
gpDataIndices(model, dimNo, blockNo)
}
\arguments{
  \item{model}{the model for which the indices are being returned.}
  \item{dimNo}{the dimension for which the presence of missing data is being looked at.}
  \item{blockNo}{the block number in the PITC approximation for which the indices are required.}
}
\value{
  \item{ind}{indices of training data along that dimension which isn't missing.}
}
\seealso{
\code{\link{gpCreate}}.
}
\examples{
## missing
}
\keyword{model}
