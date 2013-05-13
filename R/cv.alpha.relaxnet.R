################################################################################

## cross-validation (on both lambda and alpha) for relaxnet models
## adapted from the cv.glmnet function from package glmnet

cv.alpha.relaxnet <- function(x, y, family = c("gaussian", "binomial"),
                        nlambda = 100,
                        alpha = c(.1, .3, .5, .7, .9),

                        relax = TRUE,
                        relax.nlambda = 100,
                        relax.max.vars = min(nrow(x), ncol(x)) * 0.8,

                        lambda = NULL,
                        relax.lambda.index = NULL,
                        relax.lambda.list = NULL,

                        ##type.measure=c("mse","deviance","class","auc","mae"),
                        ## just set it in code for now
                        
                        nfolds = 10,
                        foldid,

                        ...) {

  start.time <- Sys.time()

  family = match.arg(family)

  ## do something with this if including type.measure as an arg
  
  ## if(family == "gaussian") type.measure <- "mse"
  ## if(family == "binomial") type.measure <- "deviance"

  N=nrow(x)

  if(missing(foldid)) {

    ## check nfolds here

    foldid <- sample(rep(seq(nfolds), length=N))

  }## else {

    ## check foldid

##  }
    
  if(relax && ncol(x) == 1) {

    warning("x has only one column, setting relax to FALSE")
    relax <- FALSE
  }

  cv.relaxnet.results <-

    lapply(alpha,
           function(alpha.val) {

             cv.relaxnet(x, y, family,
                         nlambda,
                         alpha.val,
                         relax,
                         relax.nlambda,
                         relax.max.vars,
                         lambda,
                         relax.lambda.index,
                         relax.lambda.list,
                         foldid = foldid,
                         ...)
           })

  alpha.min.index <- which.min(sapply(cv.relaxnet.results, function(cv.result) cv.result$min.cvm))

  end.time <- Sys.time()
  
  obj <- list(call = match.call(),
              relax = relax,
              alpha = alpha,
              cv.relaxnet.results = cv.relaxnet.results,
              which.alpha.min = alpha[alpha.min.index],
              total.time = as.double(difftime(end.time, start.time,
                                              units = "secs")))

  class(obj) <- "cv.alpha.relaxnet"

  obj
}
