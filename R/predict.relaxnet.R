
## predict function for relaxnet objects -- currently uses colnames
## of beta matrix in order to subset newx -- may need to change this

predict.relaxnet <- function(object,
                             newx,
                             which.model,
                             s = NULL,
                             type = c("link", "response", "coefficients",
                                     "nonzero", "class"),
                             exact = FALSE,
                             ##offset,
                             ...) {

  type = match.arg(type)

  if(!identical(exact, FALSE))
    stop("For the current version, the only supported value of\n",
         "the exact argument is FALSE")
  
  if(which.model == "main") {
    return(predict(object$main.glmnet.fit,
                   newx = newx,
                   s = s,
                   type = type,
                   exact = exact,
                   ## offset = offset,
                   ...))
  }

  ## otherwise we want predictions from one of the relaxed models
  

  ## intercept only model
  
  if(inherits(object$relax.glmnet.fits[[which.model]],
              "relaxnet.intercept.only")) {

    int.model <- object$relax.glmnet.fits[[which.model]]

    if(type == "coefficients") {

      if(length(s) != 1) stop("for type = coefficients, s must have length 1")

      return(predict(object$main.glmnet.fit,
                     s = object$main.glmnet.fit$lambda[1],
                     type = "coefficients",
                     ...))
    }

    return(switch(class(object$main.glmnet.fit)[1],
                  elnet =
                    if(type %in% c("link", "response")) {
                      return(rep(int.model,
                                 nrow(newx)))
                    } else {
                           stop("type = ", type,
                                "has not been implemented yet\n",
                                "for intercept only elnet models")
                    },

                  lognet =
                    switch(type,
                           link = rep(int.model,
                             nrow(newx)),
                           
                           response = 1 / (1 +
                             exp(-rep(int.model,
                                      nrow(newx)))),
                           ## check to make sure I'm returning the
                           ## class the right way (i.e. as 0's/1's)
                           class = rep(ifelse(int.model>=0, 1, 0),
                                       nrow(newx)),
                           stop("type = ", type,
                                "has not been implemented yet\n",
                                "for intercept only lognet models"))))

  }

  ## if type is coefs, we need to make sure that length of result is number of
  ## vars in main fit, not relaxed fit, and fill in rest of vars with 0's
  ## just allow one value of s for now
  
  if(type == "coefficients") {

    if(length(s) != 1) stop("for type = coefficients, s must have length 1")

    ## just do this to get a single column matrix of the right size

    result <- predict(object$main.glmnet.fit, type = "coefficients",
                      s = object$main.glmnet.fit$lambda[1])

    ## now get the actual coefs

    ## subset the columns of newx to conform with the relaxed model

    actual.coefs <- predict(object$relax.glmnet.fits[[which.model]],
                    newx = newx[,
                      rownames(object$relax.glmnet.fits[[which.model]]$beta)],
                    s = s,
                    type = type,
                    exact = exact,
                    ##offset = offset,
                    ...)

    ## set all to zero, keeping dimnames
    
    result[] <- 0

    ## replace relevant ones with actual
    
    result[rownames(actual.coefs), ] <- actual.coefs

    return(result)
  }

  ## otherwise type is not coefs and length will be OK
  ## (i.e. not different from main)
  
  predict(object$relax.glmnet.fits[[which.model]],
          newx = newx[,
            rownames(object$relax.glmnet.fits[[which.model]]$beta)],
          s = s,
          type = type,
          exact = exact,
          ...)  
}
