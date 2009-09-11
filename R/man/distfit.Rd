\name{distfit}
\Rdversion{1.0}
\alias{distfit}
\alias{distfit_obj}
\title{Transforming the distributions over log-expression to means and variances of absolute expression}
\description{
  The distributions over log-expression are transformed to means and variances of absolute expression by minimizing the squared deviation of 5 %, 25 %, 50 %, 75 % and 95 % quantiles of a Gaussian or Gamma distribution to corresponding quantiles from mmgMOS.

  'distfit' calls the optimization and determines the used distribution.

  'distfit_obj' is the optimized function, returning the squared deviation of the quantiles of the distribution and the data.
}
\usage{
t <- distfit(data, dist = "normal")
r <- distfit_obj(theta, y, cdf)
}
\arguments{
  \item{data}{The preprocessed data given to 'distfit'.}
  \item{dist}{The distribution used in optimizing the quantiles. Possible values for this variable are "normal", referring to the Gaussian distribution, and "gamma", referring to the Gamma distribution. Other values result in an error.}
  \item{theta}{An array consisting of two values: the first is the mean of the distribution and the second the standard deviation.}
  \item{y}{The preprocessed data given to 'distfit_obj'.}
  \item{cdf}{The quantile function of the desired distribution, determined by the variable dist.}
}
\value{
  'distfit' returns the means and variances of absolute expression.

  'distfit_obj' returns the squared deviation of the quantiles of the distribution and the data.
}
\seealso{
  \code{\link{processData}}.
}
\examples{
## missing
}
\keyword{model}
