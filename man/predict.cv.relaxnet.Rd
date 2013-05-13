\name{predict.cv.relaxnet}
\alias{predict.cv.relaxnet}
\alias{predict.cv.alpha.relaxnet}

\title{
Predict Methods for cv.relaxnet and cv.alpha.relaxnet Objects
}

\description{
Similar to other predict methods, this functions predicts fitted values, logits, coefficients and more from a fitted \code{"cv.relaxnet"} or \code{"cv.alpha.relaxnet"} object. By default, predictions are made at those values of the tuning parameters which "won" the cross-validation.
}

\usage{
\method{predict}{cv.relaxnet}(object,
        newx,
        which.model = object$which.model.min,
        s = object$overall.lambda.min,
        type = c("link", "response", "coefficients", "nonzero", "class"),
        exact = FALSE,
        ...)

\method{predict}{cv.alpha.relaxnet}(object,
        newx,
        alpha.val = object$which.alpha.min,
        type = c("link", "response", "coefficients", "nonzero", "class"),
        ...)
}

\arguments{
  \item{object}{
The object from which predictions are to be made.
}
  \item{newx}{
Matrix of new values for \code{x} at which predictions are to be made. Must be a matrix; can be sparse as in \code{Matrix} package. This argument is not used for \code{type=c("coefficients","nonzero")}
}
  \item{alpha.val}{
Value of alpha at which predictions are to be made. Default is to use that value which "won" the cross-validation.
}
  \item{which.model}{
Specifies the submodel from which predictions are required. \code{"main"} indicates the main glmnet model, while an integer indicates one of the relaxed models. Default for both functions is to use the submodel which "won" the cross-validation.
}
  \item{s}{
Value(s) of the penalty parameter \code{lambda} at which predictions are required. Default for both functions is to use that value which "won" the cross-validation.
}
  \item{type}{
See \code{link[glmnet]{predict.glmnet}}.
}
  \item{exact}{
Only the default, \code{FALSE}, is supported. See \code{link[glmnet]{predict.glmnet}}.
}
  \item{\dots}{
Further arguments passed to predict.relaxnet or to predict.cv.relaxnet (for the alpha version).
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
\code{\link{relaxnet}}, \code{\link{cv.relaxnet}}, \code{\link{cv.alpha.relaxnet}}, \code{\link{predict.relaxnet}}, \code{\link{predict.glmnet}}
}

