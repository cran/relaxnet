\name{relaxnet-package}
\alias{relaxnet-package}

\docType{package}

\title{
Relaxation (as in Relaxed Lasso, Meinshausen 2007) Applied to glmnet Models
}

\description{
Extends the glmnet package with "relaxation", done by running glmnet once on the entire predictor matrix, then again on each different subset of variables from along the regularization path. Penalty may be lasso (alpha = 1) or elastic net (0 < alpha < 1). For this version, family may be "gaussian" or "binomial" only. Takes advantage of fast fortran code from the glmnet package.
}

\author{
Stephan Ritter, with design contributions from Alan Hubbard.

Much of the code (and some help file content) is adapted from the \pkg{glmnet} package, whose authors are Jerome Friedman, Trevor Hastie and Rob Tibshirani.

Maintainer: Stephan Ritter <sritter@berkeley.edu>
}

\references{
Stephan Ritter and Alan Hubbard, Tech report (forthcoming).

Jerome Friedman, Trevor Hastie, Rob Tibshirani (2010) \dQuote{Regularization Paths for Generalized Linear Models via Coordinate Descent.} \emph{Journal of Statistical Software} \bold{33}(1)

Nicolai Meinshausen (2007) \dQuote{Relaxed Lasso} \emph{Computational Statistics and Data Analysis} \bold{52}(1), 374-393
}

% \keyword{ package }

\seealso{
\code{\link{relaxnet}}, \code{\link{cv.relaxnet}}, \code{\link{glmnet}}
}

\note{
This is a preliminary release and several additional features are planned for later versions.
}