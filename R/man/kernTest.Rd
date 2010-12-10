\name{kernTest}
\Rdversion{1.0}
\alias{kernTest}
\title{Run some tests on the specified kernel.}
\description{
  runs some tests on a given kernel structure to ensure it is correctly implemented.
}
\usage{
  values <- kernTest(kernType, numIn=4, tieParamNames=list(), testindex=NULL)
}
\arguments{
  \item{kernType}{type of kernel to test. For example, 'rbf' or {'cmpnd', 'rbf', 'lin', 'white'}.}
  \item{numIn}{the number of input dimensions (default is 4).}
  \item{tieParamNames}{list of regular expressions for parameter names that should be tied (default is none).}
%  \item{flags}{vector with a series of binary flags indicating the type of the kernel used. Currently there are only two flags indicating: (1) the normalisation of the kernel; (2) its stationarity. Note that for some kernels these flags may be pointless (e.g. there is not a normalised version of the white kernel and it is always stationary). In these cases the kernel computation should simply ignore the flags. By default they are set to false.}
}
\value{
  \item{kern}{the kernel as it was used in the tests. }
}
\seealso{
\code{\link{kernCreate}}.
}
\examples{
## missing
}
\keyword{model}
