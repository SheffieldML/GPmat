\name{lnDiffErfs}
\Rdversion{1.0}
\alias{lnDiffErfs}
\title{Helper function for computing the log of difference}
\description{
  Helper function for computing the log of difference
}
\usage{
values <- lnDiffErfs(x1, x2)
v <- lnDiffErfs(x1, x2)
}
\arguments{
  \item{x1}{argument of the positive erf}
  \item{x2}{argument of the negative erf}
}
\details{
  \code{v, s <- lnDiffErfs(x1, x2)}
  computes the log of the difference of two erfs in a numerically stable manner.

  \code{v <- lnDiffErfs(x1, x2)}
  computes the log of the difference of two erfs in a numerically stable manner.

}
\value{
  \item{v}{log(abs(erf(x1) - erf(x2)))}
  \item{s}{sign(erf(x1) - erf(x2))}
}
\seealso{
\code{\link{gradLnDiffErfs}}.
}
\examples{
## missing
}
\keyword{model}
