\name{linmod}
\alias{linmod}
\title{Create a simple linear model}
\usage{
  linmod(y, X, V = diag(length(y)))
}
\description{
  This creates a simple linear model.
}
\examples{
X <- matrix(c(rep(1, 20), 1:20), ncol = 2)
y <- X \%*\% c(1, 0.7) + rnorm(20)
o <- linmod(y, X)
}

