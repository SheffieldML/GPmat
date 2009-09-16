\name{multiKernCacheBlock}
\Rdversion{1.0}
\alias{multiKernCacheBlock}
\title{}
\description{
  
}
\usage{
K <- multiKernCacheBlock(kern, fhandle, arg, i, j, X1, X2)
multiKernCacheBlock()
}
\arguments{
  \item{kern}{the structure containing the kernel.}
  \item{fhandle}{function handle for computing the kernel}
  \item{arg}{cell array of arguments to fhandle}
  \item{i}{the row of the block of the kernel to be computed.}
  \item{j}{the column of the block of the kernel to be computed.}
  \item{X1}{first set of kernel inputs.}
  \item{X2}{second set of kernel inputs.}
}
\details{
  \code{K <- multiKernCacheBlock(kern, fhandle, arg, i, j, X1, X2)}
  Maintains a cache of values of multiKern blocks with constant parameters and returns values

  \code{multiKernCacheBlock()}
  Clears the cache.

}
\value{
  \item{K}{the kernel matrix for the given inputs.}
}
\seealso{
\code{\link{multiKernCreate}, \link{multiKernCompute}, \link{multiKernComputeBlock}, \link{multiKernFixBlocks}}.
}
\examples{
## missing
}
\keyword{model}
