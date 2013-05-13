\name{relaxnet}
\alias{relaxnet}
\title{Relaxation (as in Relaxed Lasso, Meinshausen 2007) applied to glmnet Models}
\description{
Runs glmnet once on the full x matrix, then again on each distinct subset of columns from along the solution path. The penalty may be lasso (\code{alpha} = 1) or elastic net (0 < \code{alpha} < 1). The outcome (\code{y}) may be continuous or binary.
}
\usage{
relaxnet(x, y, family = c("gaussian", "binomial"),
         nlambda = 100,
         alpha = 1,
         relax = TRUE,
         relax.nlambda = 100,
         relax.max.vars = min(nrow(x), ncol(x)) * 0.8,
         lambda = NULL,
         relax.lambda.index = NULL,
         relax.lambda.list = NULL,
         ...)
}
\arguments{
  \item{x}{
Input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format (inherit from class \code{"sparseMatrix"} as in package \code{Matrix}). Must have unique colnames.
}
  \item{y}{
response variable. Quantitative for \code{family="gaussian"}. For \code{family="binomial"} should be either a factor with two levels, or a two-column matrix of counts or proportions.
}
  \item{family}{
Response type (see above).
}
  \item{nlambda}{
The number of \code{lambda} values - default is 100. Determines how fine the grid of lambda values should be.
}
  \item{alpha}{
Elastic net mixing parameter (see \code{\link{glmnet}}).
}
  \item{relax}{
Should the model be relaxed. If FALSE, only the main glmnet model is run and no relaxed models are.
}
  \item{relax.nlambda}{
Like nlambda but for secondary (relaxed) models.
}
  \item{relax.max.vars}{
Maximum number of variables for relaxed models. No relaxation will be done for subsets along the regularization path with number of variables greater than relax.max.vars. If \code{ncol(x)} > \code{nrow(x)} and \code{alpha} < 1, it may make sense to use a value > \code{nrow(x)}, but this may lead to increased computation time.
}
  \item{lambda}{
See (see \code{\link{glmnet}}). Optional and meant primarily for use by \code{\link{cv.relaxnet}}.
}
  \item{relax.lambda.index}{
Vector which indexes the lambda argument and specifyies the values at which a relaxed model should be fit. Optional and meant primarily for use by \code{\link{cv.relaxnet}}. Ignored if lambda argument is \code{NULL}.
}
  \item{relax.lambda.list}{
List of lambda values to use for the relaxed models. Optional and meant primarily for use by \code{\link{cv.relaxnet}}. Ignored if lambda argument is \code{NULL}.
}
  \item{\dots}{
Further aruments passed to glmnet. Use with caution as this has not yet been tested. For example, setting \code{standardize = FALSE} will probably work correctly, but setting an offset probably won't.
}
}

\value{
Object of class code{"relaxnet"} with the following components:

  \item{call}{
A copy of the call which produced this object
}
  \item{main.glmnet.fit}{
The object resulting from running \code{\link{glmnet}} on the entire \code{x} matrix.
}
  \item{relax}{
The value of the relax argument. If this is \code{FALSE}, then several of the other elements of this result will be set to \code{NA}.
}
  \item{relax.glmnet.fits}{
A list containing the secondary \code{\link{glmnet}} fits gotten by running \code{\link{glmnet}} on the distinct subsets of the columns of \code{x} resulting along the solution path of lambda values.
}
  \item{relax.num.vars}{
Vector giving the number of variables in each "relaxed" model.
}
  \item{relax.lambda.index}{
This vector indexes result$main.glmnet.fit$lambda and gives the lambda values at which the relax.glmnet.fits were obtained.
}
  \item{total.time}{
Total time in seconds to produce this result.
}
  \item{main.fit.time}{
Time in seconds to produce the main glmnet fit.
}
  \item{relax.keep}{
In certain cases some of the relaxed models are removed after fitting. \code{relax.fit.times} records times for these removed models as well. This logical vector shows which of the models whose timings are given in \code{relax.fit.times} were actually kept and have results given in \code{relax.glmnet.fits} above. Hopefully this will not be necessary in later versions.
}
  \item{relax.fit.times}{
Vector of times in seconds to produce secondary "relaxed" models.
}
}

\references{
Stephan Ritter and Alan Hubbard, Tech report (forthcoming).

Jerome Friedman, Trevor Hastie, Rob Tibshirani (2010) \dQuote{Regularization Paths for Generalized Linear Models via Coordinate Descent.} \emph{Journal of Statistical Software} \bold{33}(1)

Nicolai Meinshausen (2007) \dQuote{Relaxed Lasso} \emph{Computational Statistics and Data Analysis} \bold{52}(1), 374-393
}

\author{
Stephan Ritter, with design contributions from Alan Hubbard.

Much of the code (and some help file content) is adapted from the \pkg{glmnet} package, whose authors are Jerome Friedman, Trevor Hastie and Rob Tibshirani.
}

\seealso{
\code{\link[glmnet]{glmnet}}, \code{\link{cv.relaxnet}}, \code{\link{predict.relaxnet}}
}

\note{
This is a preliminary release and several additional features are planned for later versions.
}

\examples{
## generate predictor matrix

nobs <- 100
nvars <- 200

set.seed(23)

x <- matrix(rnorm(nobs * nvars), nobs, nvars)

## make sure it has unique colnames

colnames(x) <- paste("x", 1:ncol(x), sep = "")

## let y depend on first 5 columns plus noise

y <- rowSums(x[, 1:5]) + rnorm(nrow(x))

## default is family = "gaussian"

result1 <- relaxnet(x, y)

summary(result1)

## now fit family = "binomial" model

y.bin <- rbinom(nrow(x), 1, prob = plogis(0.2 * rowSums(x[, 1:5])))

result2 <- relaxnet(x, y.bin, family = "binomial")

summary(result2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
