\name{summary.relaxnet}
\alias{summary.relaxnet}
\alias{print.summary.relaxnet}

\title{
Generate and print summaries of class \code{"relaxnet"} objects.
}

\description{
Print summaries of relaxnet objects.
}

\usage{
\method{summary}{relaxnet}(object, ...)
\method{print}{summary.relaxnet}(x, digits = 3, ...)
}

\arguments{
  \item{object}{
The \code{"relaxnet"} object to summarize.
}
  \item{x}{
An object of class \code{"summary.relaxnet"}.
}
  \item{digits}{
Passed to \code{print.default}.
}
  \item{\dots}{
For the summary method: ignored.
For the print method: passed to \code{print.default}.
}
}

\value{
For the summary method: an object of type \code{"summary.relaxnet"} containing a subset of the elements of \code{object}.

For the print method: returns \code{x} invisibly.
}

\author{
Stephan Ritter, with design contributions from Alan Hubbard.

Much of the code (and some help file content) is adapted from the \pkg{glmnet} package, whose authors are Jerome Friedman, Trevor Hastie and Rob Tibshirani.
}

\seealso{
\code{\link{relaxnet}}, \code{\link{print.relaxnet}}
}
