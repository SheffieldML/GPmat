\name{disimXdisimKernCompute}
\Rdversion{1.0}
\alias{disimXdisimKernCompute}
\title{Compute a cross kernel between two DISIM kernels.}
\description{
  Compute a cross kernel between two DISIM kernels.
}
\usage{
k <- disimXdisimKernCompute(disimKern1, disimKern2, t)
k <- disimXdisimKernCompute(disimKern1, disimKern2, t1, t2)
}
\arguments{
  \item{disimKern1}{the kernel structure associated with the first DISIM kernel.}
  \item{disimKern2}{the kernel structure associated with the second DISIM kernel.}
  \item{t}{inputs for which kernel is to be computed.}
  \item{t1}{row inputs for which kernel is to be computed.}
  \item{t2}{column inputs for which kernel is to be computed.}
}
\details{
  \code{k <- disimXdisimKernCompute(disimKern1, disimKern2, t)}
  computes cross kernel terms between two DISIM kernels for the multiple output kernel.

  \code{k <- disimXdisimKernCompute(disimKern1, disimKern2, t1, t2)}
  computes cross kernel terms between two DISIM kernels for the multiple output kernel.

}
\value{
  \item{k}{block of values from kernel matrix.}
}
\seealso{
\code{\link{multiKernParamInit}, \link{multiKernCompute}, \link{disimKernParamInit}}.
}
\examples{
## missing
}
\keyword{model}
