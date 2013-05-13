\name{print.relaxnet}
\alias{print.relaxnet}

\title{
Print Method for relaxnet Objects
}

\description{
This function just calls \code{print(summary(x))}. See \code{\link{summary.relaxnet}}.
}

\usage{
\method{print}{relaxnet}(x, digits, ...)
}

\arguments{
  \item{x}{
The \code{"relaxnet"} object to be printed.
}
  \item{digits}{
Passed to \code{print.summary.relaxnet}.
}
  \item{\dots}{
Passed to \code{print.summary.relaxnet}.
}
}
\value{
Returns \code{x} invisibly.
}

\author{
Stephan Ritter, with design contributions from Alan Hubbard.

Much of the code (and some help file content) is adapted from the \pkg{glmnet} package, whose authors are Jerome Friedman, Trevor Hastie and Rob Tibshirani.
}

\seealso{
\code{\link{relaxnet}}, \code{\link{summary.relaxnet}}
}
