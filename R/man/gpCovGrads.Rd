\name{gpCovGrads}
\Rdversion{1.0}
\alias{gpCovGrads}
\title{Sparse objective function gradients wrt Covariance functions for inducing variables.}
\description{
  gives the gradients of the log likelihood with respect to the components of the sparse covariance (or the full covariance for the ftc case).
}
\usage{
gpCovGrads(model, M)
}
\arguments{
  \item{model}{the model for which the gradients are to be computed.}
  \item{M}{The training data for which the computation is to be made}
}
\value{
  \item{gK_uu}{the gradient of the likelihood with respect to the elements of K_uu (or in the case of the 'ftc' criterion the gradients with respect to the kernel).}
  \item{gK_uf}{the gradient of the likelihood with respect to the elements of K_uf.}
  \item{gLambda}{the gradient of the likelihood with respect to the diagonal term in the fitc approximation and the blocks of the pitc approximation.}
  \item{gBeta}{the gradient with respect to the beta term in the covariance structure. }
}
\seealso{
\code{\link{gpCreate}, \link{gpLogLikeGradients}}.
}
\examples{
## missing
}
\keyword{model}
