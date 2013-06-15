\name{cv.relaxnet}
\alias{cv.relaxnet}
\alias{cv.alpha.relaxnet}
\title{
Cross-Validation for relaxnet Models
}
\description{
Both of these functions will perform v-fold cross-validation to select tuning parameters for relaxnet models. \code{cv.relaxnet} will cross-validate on the value of lambda for both the main model and for the relaxed models, a two dimensional cross-validation. \code{cv.alpha.relaxnet} will in addition cross-validate on the value of alpha. For each value of alpha, relaxnet is run once on whole data set, then it is run again v times on subsets of the rows of the data.
}

\usage{
cv.relaxnet(x, y, family = c("gaussian", "binomial"),
            nlambda = 100,
            alpha = 1,
            relax = TRUE,
            relax.nlambda = 100,
            relax.max.vars = min(nrow(x), ncol(x)) * 0.8,
            lambda = NULL,
            relax.lambda.index = NULL,
            relax.lambda.list = NULL,
            nfolds = 10,
            foldid,
            multicore = FALSE,
            mc.cores,
            mc.seed = 123,
            ...)

cv.alpha.relaxnet(x, y, family = c("gaussian", "binomial"),
                  nlambda = 100,
                  alpha = c(.1, .3, .5, .7, .9),
                  relax = TRUE,
                  relax.nlambda = 100,
                  relax.max.vars = min(nrow(x), ncol(x)) * 0.8,
                  lambda = NULL,
                  relax.lambda.index = NULL,
                  relax.lambda.list = NULL,
                  nfolds = 10,
                  foldid,
                  multicore = FALSE,
                  mc.cores,
                  mc.seed = 123,
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
Elastic net mixing parameter (see \code{\link{glmnet}}). For \code{cv.relaxnet}, this should be a single value. For \code{cv.alpha.relaxnet} it should be a vector of values.
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
See (see \code{\link{glmnet}}). Optional: default is to let \code{glmnet} choose its own sequence.
}
  \item{relax.lambda.index}{
Vector which indexes the lambda argument and specifyies the values at which a relaxed model should be fit. Optional: default is to let \code{relaxnet} determine these values based on the \code{beta} matrix from the main \code{glmnet} fit. Ignored if \code{lambda} argument is \code{NULL}.
}
  \item{relax.lambda.list}{
List of lambda values to use for the relaxed models. Optional: default is to let \code{relaxnet} determine these values. Ignored if lambda argument is \code{NULL}.
}
  \item{nfolds}{
Number of folds - default is 10. Although \code{nfolds} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is \code{nfolds=3}.
}
  \item{foldid}{
An optional vector of values between 1 and \code{nfold} identifying what fold each observation is in. If supplied, \code{nfolds} can be missing.
}
  \item{multicore}{
Should execution be parallelized over cv folds (for \code{cv.relaxnet}) or over alpha values (for \code{cv.alpha.relaxnet}) using multicore functionality from R's parallel package?
}
  \item{mc.cores}{
Number of cores/cpus to be used for multicore processing. Processing will be most efficient if nfolds (for \code{cv.relaxnet}) or the length of alpha (for \code{cv.alpha.relaxnet}) is a multiple of mc.cores. Ignored if \code{multicore} is \code{FALSE}
}
  \item{mc.seed}{
Integer value with which to seed the RNG when using parallel processing (internally, \code{\link[base]{RNGkind}} will be called to set the RNG to \code{"L'Ecuyer-CMRG"}). Will be ignored if \code{multicore} is \code{FALSE}. If \code{mulicore} is \code{FALSE}, one should be able to get reprodicible results by setting the seed normally (with \code{\link[base]{set.seed}}) prior to running.
}
  \item{\dots}{
Further aruments passed to glmnet. Use with caution as this has not yet been tested. For example, setting \code{standardize = FALSE} will probably work correctly, but setting an offset probably won't.
}
}

\details{


cv.glmnet's \code{type.measure} argument has not yet been implemented. For type = gaussian models, mean squared error is used, and for type = binomial, binomial deviance is used.}

\value{
For \code{cv.relaxnet} -- an object of class \code{"cv.relaxnet"} containing the following slots:

  \item{call}{
A copy of the call which produced this object
}
  \item{relax}{
The value of the relax argument. If this is \code{FALSE}, then several of the other elements of this result will be set to \code{NA}.
}
  \item{lambda}{
lambda sequence used for this fit
}
  \item{cvm}{
The mean cross-validated error - a vector of length \code{length(lambda)}. For main model.
}
  \item{cvsd}{
estimate of standard error of \code{cvm} for main model.
}
  \item{cvup}{
upper curve = \code{cvm+cvsd} for main model.
}
  \item{cvlo}{
lower curve = \code{cvm-cvsd} for main model.
}
  \item{nzero}{
number of non-zero coefficients at each \code{lambda} for main model.
}
  \item{name}{
a text string indicating type of measure.
}
  \item{relaxnet.fit}{
Fitted relaxnet object for the full data.
}
  \item{relax.cvstuff.list}{
List containing cvm and cvsd for each of the relaxed models.
}
  \item{relax.lambda.list.trunc}{
List containing the values of lambda used for cross-validation for each relaxed model
}
  \item{which.model.min}{
This will have value \code{"main"} if the main model "won" the cross-validation, and if not, it will be an integer specifying which relaxed model won (i.e. which element of relaxnet.fit$relax.glmnet.fits).
}
  \item{overall.lambda.min}{
The value of lambda with overall min cvm (i.e. from the submodel specified by \code{which.model.min}).
}
  \item{min.cvm}{
The overall minimum value of cvm
}
  \item{main.lambda.min}{
lambda.min, restricted to main model only.
}
  \item{main.lambda.1se}{
lambda.1se, restricted to main model only (see \code{\link{cv.glmnet}}).
}
  \item{main.min.cvm}{
Minimum of cvm, restricted to main model only.
}
  \item{total.time}{
Time in seconds to fit this object (full data fit plus all cross-validation)
}
  \item{full.data.fit.time}{
Time in seconds to fit the relaxnet.fit.
}
  \item{cv.fit.times}{
Time in seconds to fit the models for each set of cv folds.
}


For \code{cv.alpha.relaxnet} -- an object of class \code{"cv.alpha.relaxnet"} containing the following slots:

  \item{call}{
A copy of the call which produced this object
}
  \item{relax}{
The value of the relax argument. If this is \code{FALSE}, then several of the other elements of this result will be set to \code{NA}.
}
  \item{alpha}{
The alpha values used.
}
  \item{cv.relaxnet.results}{
List of cv.relaxnet objects, one for each alpha value.
}
  \item{which.alpha.min}{
The alpha value which "won" the cross-validation.
}
  \item{total.time}{
Time in seconds to fit this object.
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
\code{\link{relaxnet}}, \code{\link{predict.cv.relaxnet}} 
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

## run cv.relaxnet

cv.result <- cv.relaxnet(x, y)

predict(cv.result, type = "nonzero")

## very few false positives compared to glmnet alone

## glmnet min rule

predict(cv.result$relaxnet.fit$main.glmnet.fit,
        type = "nonzero",
        s = cv.result$main.lambda.min)

## glmnet 1se rule

predict(cv.result$relaxnet.fit$main.glmnet.fit,
        type = "nonzero",
        s = cv.result$main.lambda.1se)

## get values of the coefs for cv.relaxnet's chosen fit

coefs <- drop(predict(cv.result, type = "coef"))

coefs[coefs != 0]
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
