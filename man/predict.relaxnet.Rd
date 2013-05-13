\name{predict.relaxnet}
\alias{predict.relaxnet}

\title{
Predict Method for \code{"relaxnet"} Objects
}

\description{
Similar to other predict methods, this functions predicts fitted values, logits,
coefficients and more from a fitted \code{"relaxnet"} object.
}

\usage{
\method{predict}{relaxnet}(object,
        newx,
        which.model,
        s = NULL,
        type = c("link", "response", "coefficients", "nonzero", "class"),
        exact = FALSE,
	...)
}

\arguments{
  \item{object}{
The \code{"relaxnet"} object from which to get predictions.
}
  \item{newx}{
Matrix of new values for \code{x} at which predictions are to be made. Must be a matrix; can be sparse as in \code{Matrix} package. This argument is not used for \code{type=c("coefficients","nonzero")}
}
  \item{which.model}{
Specifies the submodel from which predictions are required. \code{"main"} indicates the main glmnet model, while an integer indicates one of the relaxed models.
}
  \item{s}{
Value(s) of the penalty parameter \code{lambda} at which predictions are required. Default is the entire sequence of lambda values for the model specified by \code{which.model}.
}
  \item{type}{
See \code{link[glmnet]{predict.glmnet}}.
}
  \item{exact}{
Only the default, \code{FALSE}, is supported. See \code{link[glmnet]{predict.glmnet}}.
}
  \item{\dots}{
Further arguments passed to \code{predict.glmnet}. In the current version, these are not guaranteed to work correctly (for example, \code{offset} has not yet been implemented for relaxnet).
}
}

\value{
The object returned depends on type.
}

\author{
Stephan Ritter, with design contributions from Alan Hubbard.

Much of the code (and some help file content) is adapted from the \pkg{glmnet} package, whose authors are Jerome Friedman, Trevor Hastie and Rob Tibshirani.
}

\seealso{
\code{\link{relaxnet}}, \code{\link{predict.glmnet}}
}
