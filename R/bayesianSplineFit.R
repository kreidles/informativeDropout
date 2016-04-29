#
# Supporting functions for the Spline model / RJMCMC loop
#
#
#

#' Data stored with each iteration of the Bayesian spline model
#' during the RJMCMC run
#'
#' @param knots list of knot positions by group
#' @param Theta list of spline coefficients by group
#' @param betas.covariates regression coefficients related to covariates
#' @param sigma.spline covariance of the spline coefficients
#' @param sigma.randomIntercept variance of the random intercepts
#' @param sigma.randomSlope variance of the random slopes
#' @param sigma.randomInterceptSlope covariance of the random intercepts and slopes
#' @param sigma.error for Gaussian outcomes, the residual error
#' @param sigma.error.shape for Gaussian outcomes, the shape hyperparamter of the inverse gamma
#' distribution of the residual error
#' @param sigma.error.rate for Gaussian outcomes, the rate hyperparamter of the inverse gamma
#' distribution of the residual error

#' @param shape.tau
#'
#'
#' @exportClass bayes.splines.iteration
#'
bayes.splines.iteration <- function(knots=NULL, Theta=NULL, betas.covariates=NULL,
                                    sigma.spline = 1,
                                    sigma.randomIntercept = 1,
                                    sigma.randomSlope = 1, sigma.randomInterceptSlope = 0,
                                    sigma.error = 1, 
                                    sigma.error.shape = 0.001, sigma.error.rate = 0.001) {
  mi = list(
    
    #
    # The following items are separated into a list 
    # with each item representing the value for a treatment group
    #
    
    # knots per group
    knots = knots,
    # per group intercept and spline coefficients
    Theta = Theta,
    
    #
    # The following items are shared across the groups
    #
    
    # regression coefficients associated with the shared covariates 
    # (includes group offsets)
    betas.covariates = betas.covariates,
    
    # variance of spline coefficients
    sigma.spline = sigma.spline,
    # variance / covariance of the random effects
    sigma.randomIntercept = sigma.randomIntercept,
    sigma.randomSlope = sigma.randomSlope,
    sigma.randomInterceptSlope = sigma.randomInterceptSlope,
    
    # residual error variance
    sigma.error = sigma.error,
    # hyperparameters describing the inverse gamma distribution of the variance components
    sigma.error.shape = sigma.error.shape,
    sigma.error.rate = sigma.error.rate,
    
    # indicates which changes to the model were proposed at this iteration
    proposed=list(knot.add=FALSE, knot.remove=FALSE, knot.move=FALSE,
                  fixedEffects=FALSE, fixedEffectsCovariates=FALSE),
    # incidates which changes were accepted at this iteration
    accepted=list(knot.add=FALSE, knot.remove=FALSE, knot.move=FALSE,
                  fixedEffects=FALSE, fixedEffectsCovariates=FALSE)
  )
  
  class(mi) <- append(class(mi), "bayes.splines.iteration")
  return(mi)
  
}


#'
#' Simulation and model options for the natural b-spline model
#' 
#' @param iterations number of iterations for the MCMC simulation
#' @param burnin burn in period for the simulation, i.e. the number of 
#' iterations to throw away at the beginning of the simulation
#' @param thin thinning interval, i.e. if thin=n, only keep every nth iteration
#' @param knots.prob.birth probability of adding a new knot to the model on a given iteration
#' @param knots.min minimum number of knots in the model. Must be greater than or equal to 1.
#' @param knots.max maximum number of knots in the model.
#' 
#' @exportClass bayes.splines.model.options
#' 
bayes.splines.model.options = function(iterations=10000, burnin=500, thin=NA,
                                       knots.prob.birth=0.5, knots.min=1, knots.max=NA,
                                       knots.positions.start=NULL, knots.positions.candidate=NULL,
                                       prob.min=0.00001, prob.max=0.99999,
                                       dropout.estimationTimes=NULL,
                                       shape.tau=NULL, rate.tau=NULL,
                                       sigma.beta=NULL, lambda.numKnots=NULL,
                                       sigma.error=0,
                                       sigma.error.df = NULL,
                                       sigma.error.scale = NULL,
                                       eta.null=NULL) {
  # TODO: add na.rm to ignore subjects with missing outcomes or covariates
  # validate the mcmc options
  if (is.na(iterations) || is.null(iterations) || iterations <= 1) {
    stop("model options error :: invalid number of iterations")
  }
  if (!is.na(burnin) && !is.null(iterations) && burnin >= iterations) {
    stop("model options error :: the burn in period must be less than the total number of iterations")
  }
  if (!is.na(thin) && !is.null(thin)  && thin >= iterations) {
    stop("model options error :: thinning value must be less that the iterations")
  }
  
  # validate the knot options.
  if (is.na(knots.prob.birth) || is.null(knots.prob.birth) || knots.prob.birth <= 0 || knots.prob.birth >= 1) {
    stop("model options error :: knot birth probability must be a value between 0 and 1.")
  }
  if (is.na(knots.min) || is.null(knots.min) || knots.min <= 0) {
    stop("model options error :: The minimum number of knots must be 1 or greater.")
  }
  if (is.na(knots.max) || is.null(knots.max) || knots.max <= knots.min) {
    stop("model options error :: The maximum number of knots must be greater than the minimum number of knots.")
  }
  
  
  # validate the prior options
  if (is.na(shape.tau) || is.null(shape.tau) || shape.tau <= 0) {
    stop("Prior options error :: shape.tau must be greater than 0")
  }
  if (is.na(rate.tau) || is.null(rate.tau) || rate.tau <= 0) {
    stop("Prior options error :: rate.tau must be greater than 0")
  }
  if (is.na(sigma.error.df) || is.null(sigma.error.df) || sigma.error.df <= 0) {
    stop("Prior options error :: sigmaError.df must be greater than 0")
  }
  if (is.null(sigma.error.scale) || 
      !is.positive.definite(sigma.error.scale)) {
    stop("Prior options error :: sigma.error.scale must be a positive definite matrix")
  }
  
  opts = list(
    # mcmc iterations
    iterations = iterations,
    # burn in period
    burnin = burnin,
    # thinning interval
    thin = thin,
    # probability of adding a knot on a given iteration
    knots.prob.birth = knots.prob.birth,
    # minimum number of knots
    knots.min = knots.min,
    # maximum number of knots,
    knots.max = knots.max,
    # starting positions for the knots
    knots.positions.start = knots.positions.start,
    # candidate positions for the knots
    knots.positions.candidate = knots.positions.candidate,
    # minimum/maximum probability for Metropolis Hastings for binary outcomes
    # these control numeric instability with taking logs
    prob.min <- prob.min,
    prob.max <- prob.max,
    
    # times at which the dropout time dependent slopes will be estimated
    dropout.estimationTimes = dropout.estimationTimes,
    
    # hyperparameters for the shape parameter
    shape.tau = shape.tau, 
    #
    rate.tau = rate.tau,
    sigma.error.df = sigma.error.df,
    sigma.error.scale = sigma.error.scale
  )
  
  class(opts) <- append(class(opts), "bayes.splines.model.options")
  return(opts)
  
}

#'
#' For binary outcomes, add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model.
#'
#'
#'
#'
#'
#'
addKnot.binary <- function(knots.previous, outcomes, 
                           times.dropout, times.observation, 
                           covariates, X.previous, Theta.previous,
                           Z, alpha, betaCovariates, 
                           sigma.beta, sigma.residual, lambda.numKnots, eta.null,
                           model.options) {
  
  # add a knot by randomly selecting a candidate knot
  candidatesPositions = model.options$candidatePositions[! model.options$candidatePositions %in% knots.previous]
  newKnot.value = sample(candidatesPositions, 1)
  knots.star <- sort(c(knots.previous, newKnot.value))
  # get the interior and boundary knots, and grab the position of the knot that
  # was just added
  knots.boundary = range(knots.star)
  knots.interior = knots.star[-c(1,length(knots.star))] 
  newKnot.position = which(knots.star == newKnot.value)
  
  # Calculate spline transformation of dropout time and create the proposed X matrix
  X.star <- cbind(
    rep(1,length(times.observation)),
    ns(times.dropout, knots=knots.interior, Boundary.knots=knots.boundary, intercept=T) * times.observation
  )
  
  # Calculate y-random effects for least squares calculations
  y = outcomes
  if (!is.null(covariates)) {
    cBeta = as.vector(as.matrix(covariates) %*% betaCovariates)
  }
  zAlpha = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
  
  # calculate residuals
  eta.wls <- eta0 + zAlpha + ifelse(!is.null(covariates), cBeta, 0)
  Theta.LSXprev <- wls(y, X, eta.wls, model.options)
  Theta.LSXstar <- wls(y, X, eta.wls, model.options)
  Theta.LSresid <- Theta.previous - Theta.LSXprev
  
  #Draw a residual for the added coefficient and calculate coefficient transformation
  residual <- rnorm(1, 0, sqrt(sigma.residual))
  # adjust position by 1 since first coefficient is the group intercept
  beta.newKnot <- Theta.LSXstar[(newKnot.position+1)] + residual
  beta.other <- Theta.LSXstar[-(newKnot.position+1)] + Theta.LSresid
  # insert the new beta value in the correct position
  if (newKnot.position == 1) {
    Theta.star = c(beta.other[1], beta.newKnot, beta.other[2:length(beta.other)])
  } else if (newKnot.position == length(knots.star)) {
    Theta.star = c(beta.other, beta.newKnot)
  } else {
    Theta.star = c(beta.other[1:(newKnot.position)], beta.newKnot, 
                   beta.other[(newKnot.position+1):length(beta.other)])
  }
  
  # Calculate birth and death probabilities                                                        
  probBirth <- ifelse(length(knots.previous) == model.options$min, 1, model.options$birthProbability)
  probDeath <- ifelse(length(knots.previous) == model.options$max - 1, 1, 1 - model.options$birthProbability)
  
  # calculate the previous eta and associated probability
  eta.previous = as.vector(X.previous %*% Theta.previous + cBeta + zAlpha)
  prob.previous = inv.logit(eta.previous)
  # adjust probabilities within tolerance levels
  prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
  prob.previous[prob.previous < model.options$prob.max] <-  model.options$prob.max
  loglikelihood.previous <- sum(log((1 - prob.previous[outcomes==0]))) + sum(log(prob.previous[outcomes==1]))    
  
  
  # calculate the proposal eta and associated probability
  eta.star = as.vector(X.star %*% Theta.star + cBeta + zAlpha)
  prob.star = inv.logit(eta.star)
  # adjust probabilities within tolerance levels
  prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
  prob.star[prob.star < model.options$prob.max] <-  model.options$prob.max
  loglikelihood.star <- sum(log((1 - prob.star[outcomes==0]))) + sum(log(prob.star[outcomes==1]))    
  
  #Calculate Acceptance Probability  
  rho <- (log(lambda.numKnots) - 
            log(length(knots.previous)) + 
            log(probDeath) - log(probBirth) + 
            log(sqrt(sigma.residual)) - 
            0.5 * log(sigma.beta) +
            (residual^2 / (2 * sigma.residual)) + 
            ((crossprod(Theta.previous) - crossprod(Theta.star)) / (2 * sigma.beta)) +
            loglikelihood.star - loglikelihood.previous
  )
  
  if (rho > log(runif(1))) {
    return (list(X=X.star, knots=knots.star, Theta=Theta.star, accepted=TRUE))
  } else {
    return (list(X=X.previous, knots=knots.previous, Theta=Theta.previous, accepted=FALSE))          
  }
}




#' Add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param 
#' @param 
#' @return 
#' @examples
#' 
addKnot.gaussian <- function(knots.previous, outcomes, 
                             times.dropout, times.observation, 
                             covariates, X.previous, Theta.previous,
                             Z, alpha, betaCovariates, 
                             sigma.beta, sigma.residual, sigma.error, lambda.numKnots) {
  
  # add a knot by randomly selecting a candidate knot
  candidatesPositions = model.options$candidatePositions[! model.options$candidatePositions %in% knots.previous]
  newKnot.value = sample(candidatesPositions, 1)
  knots.star <- sort(c(knots.previous, newKnot.value))
  # get the interior and boundary knots, and grab the position of the knot that
  # was just added
  knots.boundary = range(knots.star)
  knots.interior = knots.star[-c(1,length(knots.star))] 
  newKnot.position = which(knots.star == newKnot.value)
  
  # Calculate spline transformation of dropout time and create the proposed X matrix
  X.star <- cbind(
    rep(1,length(times.observation)),
    ns(times.dropout, knots=knots.interior, Boundary.knots=knots.boundary, intercept=T) * times.observation
  )
  
  # Calculate y-random effects for least squares calculations
  y = outcomes
  if (!is.null(covariates)) {
    cBeta = as.vector(as.matrix(covariates) %*% betaCovariates)
  }
  zAlpha = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
  
  # get the residuals
  yls <- as.vector(y - zAlpha)
  if (!is.null(covariates)) {
    yls <- (yls - cBeta)
  } 
  # Calculate least squares estimates for coefficients and differences between LS and current coefficients
  Theta.LSXprev <- ginv(crossprod(X.previous))%*%(crossprod(X.previous,yls))
  Theta.LSXstar <- ginv(crossprod(X.star))%*%(crossprod(X.star,yls))
  
  Theta.LSresid <- Theta.previous - Theta.LSXprev
  
  #Draw a residual for the added coefficient and calculate coefficient transformation
  residual <- rnorm(1, 0, sqrt(sigma.residual))
  # adjust position by 1 since first coefficient is the group intercept
  beta.newKnot <- Theta.LSXstar[(newKnot.position+1)] + residual
  beta.other <- Theta.LSXstar[-(newKnot.position+1)] + Theta.LSresid
  # insert the new beta value in the correct position
  if (newKnot.position == 1) {
    Theta.star = c(beta.other[1], beta.newKnot, beta.other[2:length(beta.other)])
  } else if (newKnot.position == length(knots.star)) {
    Theta.star = c(beta.other, beta.newKnot)
  } else {
    Theta.star = c(beta.other[1:(newKnot.position)], beta.newKnot, 
                   beta.other[(newKnot.position+1):length(beta.other)])
  }
  
  # Calculate birth and death probabilities                                                        
  probBirth <- ifelse(length(knots.previous) == model.options$min, 1, model.options$birthProbability)
  probDeath <- ifelse(length(knots.previous) == model.options$max - 1, 1, 1 - model.options$birthProbability)
  
  # Calculate residuals for likelihood ratio
  if (!is.null(covariates)) {
    LRresid.star <- as.vector(y - X.star %*% Theta.star - cBeta - zAlpha)
    LRresid.prev <- as.vector(y - X.previous %*% Theta.previous - cBeta - zAlpha)
  } else {
    LRresid.star <-as.vector(y - X.star %*% Theta.star - zAlpha)
    LRresid.prev <- as.vector(y - X.previous %*% Theta.previous - zAlpha)
  }
  
  #Calculate Acceptance Probability                                                        
  rho <- (log(lambda.numKnots) - 
            log(length(knots.previous)) + 
            log(probDeath) - log(probBirth) + 
            log(sqrt(sigma.residual)) - 
            0.5 * log(sigma.beta) +
            (residual^2 / (2 * sigma.residual)) + 
            ((crossprod(Theta.previous) - crossprod(Theta.star)) / (2 * sigma.beta)) +
            ((crossprod(LRresid.prev) - crossprod(LRresid.star)) / (2 * sigma.error))
  )
  
  if (rho > log(runif(1))) {
    return (list(X=X.star, knots=knots.star, Theta=Theta.star, accepted=TRUE))
  } else {
    return (list(X=X.previous, knots=knots.previous, Theta=Theta.previous, accepted=FALSE))          
  }
}


#' Remove a knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param 
#' @param 
#' @return 
#' @examples
#' 
removeKnot.binary <- function(dist, knots.previous, outcomes, times.dropout, 
                              times.observation, covariates,
                              X.previous, Theta.previous, Z, alpha, betaCovariates, 
                              sigma.beta, sigma.residual, lambda.numKnots, eta.null,
                              model.options) {
  
  # randomly remove an existing knot
  index = sample(1:length(knots.previous), 1)
  knots.star <- knots.previous[-index]
  
  knots.boundary = range(knots.star)
  knots.interior = knots.star[-c(1,length(knots.star))] 
  
  if (length(knots.star) > 1 ) {  
    X.star<-cbind(
      rep(1, length(times.dropout)),
      ns(times.dropout, knots=knots.interior, Boundary.knots=knots.boundary, 
         intercept=T) * times.observation
    )
  } else {
    X.star<-cbind(rep(1, length(times.observation)), times.observation)
  }
  
  # Calculate residuals
  y = as.matrix(outcomes)
  if (!is.null(covariates)) {
    cBeta = as.vector(as.matrix(covariates) %*% betaCovariates)
  }
  zAlpha = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
  
  yls <- as.vector(y - zAlpha)
  if (!is.null(covariates)) {
    yls <- (yls - cBeta)
  } 
  
  # calculate residuals
  eta.wls <- eta0 + zAlpha + ifelse(!is.null(covariates), cBeta, 0)
  Theta.LSXprev <- wls(y, X, eta.wls, model.options)
  Theta.LSXstar <- wls(y, X, eta.wls, model.options)
  Theta.LSresid <- Theta.previous - Theta.LSXprev
  
  # update the coefficients
  residual.deletedKnot <- Theta.LSresid[index]
  Theta.star <- Theta.LSXstar + Theta.LSresid[-index]
  
  
  # calculate the previous eta and associated probability
  eta.previous = as.vector(X.previous %*% Theta.previous + cBeta + zAlpha)
  prob.previous = inv.logit(eta.previous)
  # adjust probabilities within tolerance levels
  prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
  prob.previous[prob.previous < model.options$prob.max] <-  model.options$prob.max
  loglikelihood.previous <- sum(log((1 - prob.previous[outcomes==0]))) + sum(log(prob.previous[outcomes==1]))    
  
  # calculate the proposal eta and associated probability
  eta.star = as.vector(X.star %*% Theta.star + cBeta + zAlpha)
  prob.star = inv.logit(eta.star)
  # adjust probabilities within tolerance levels
  prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
  prob.star[prob.star < model.options$prob.max] <-  model.options$prob.max
  loglikelihood.star <- sum(log((1 - prob.star[outcomes==0]))) + sum(log(prob.star[outcomes==1]))    
  
  # Calculate birth and death probabilities                                                        
  probBirth <- ifelse(length(knots.star) == model.options$min, 
                      1, model.options$birthProbability)
  probDeath <- ifelse(length(knots.previous) == model.options$max, 
                      1, 1 - model.options$birthProbability)
  
  #Calculate Acceptance Probability  
  rho <- (-log(lambda.numKnots) + 
            log(length(knots.star)) -  
            log(probDeath) + log(probBirth) - 
            log(sqrt(sigma.residual)) + 
            0.5 * log(sigma.beta) -
            (residual.deletedKnot^2 / (2 * sigma.residual)) + 
            ((crossprod(Theta.previous) - crossprod(Theta.star)) / (2 * sigma.beta)) +
            loglikelihood.star - loglikelihood.previous
  )
  
  if (rho > log(runif(1))) {
    return (list(X=X.star, knots=knots.star, Theta=Theta.star, accepted=TRUE))
  } else {
    return (list(X=X.previous, knots=knots.previous, Theta=Theta.previous, accepted=FALSE))         
  }
  
}




#' Remove a knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param 
#' @param 
#' @return 
#' @examples
#' 
removeKnot.gaussian <- function(dist, knots.previous, model.options, 
                                outcomes, times.dropout, times.observation, covariates,
                                X.previous, Theta.previous, 
                                Z, alpha, betaCovariates, sigma.residual,  
                                sigma.error, sigma.beta, lambda.numKnots) {
  
  # randomly remove an existing knot
  index = sample(1:length(knots.previous), 1)
  knots.star <- knots.previous[-index]
  
  knots.boundary = range(knots.star)
  knots.interior = knots.star[-c(1,length(knots.star))] 
  
  if (length(knots.star) > 1 ) {  
    X.star<-cbind(
      rep(1, length(times.dropout)),
      ns(times.dropout, knots=knots.interior, Boundary.knots=knots.boundary, 
         intercept=T) * times.observation
    )
  } else {
    X.star<-cbind(rep(1, length(times.observation)), times.observation)
  }
  
  # Calculate residuals
  y = as.matrix(outcomes)
  if (!is.null(covariates)) {
    cBeta = as.vector(as.matrix(covariates) %*% betaCovariates)
  }
  zAlpha = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
  
  yls <- as.vector(y - zAlpha)
  if (!is.null(covariates)) {
    yls <- (yls - cBeta)
  } 
  
  # Calculate least squares estimates for coefficients and differences 
  # between LS and current coefficients
  Theta.LSXprev <- ginv(crossprod(X.previous))%*%(crossprod(X.previous, yls))
  Theta.LSXstar <- ginv(crossprod(X.star))%*%(crossprod(X.star, yls))
  Theta.LSresid <- Theta.previous - Theta.LSXprev
  # update the coefficients
  residual.deletedKnot <- Theta.LSresid[index]
  Theta.star <- Theta.LSXstar + Theta.LSresid[-index]
  
  # Calculate residuals for likelihood ratio
  if (!is.null(covariates)) {
    LRresid.star <- as.vector(y - X.star %*% Theta.star - cBeta - zAlpha)
    LRresid.prev <- as.vector(y - X.previous %*% Theta.previous - cBeta - zAlpha)
  } else {
    LRresid.star <-as.vector(y - X.star %*% Theta.star - zAlpha)
    LRresid.prev <- as.vector(y - X.previous %*% Theta.previous - zAlpha)
  }
  
  # Calculate birth and death probabilities                                                        
  probBirth <- ifelse(length(knots.star) == model.options$min, 
                      1, model.options$birthProbability)
  probDeath <- ifelse(length(knots.previous) == model.options$max, 
                      1, 1 - model.options$birthProbability)
  
  #Calculate Acceptance Probability                                                        
  rho <- (-log(lambda.numKnots) + 
            log(length(knots.star)) -  
            log(probDeath) + log(probBirth) - 
            log(sqrt(sigma.residual)) + 
            0.5 * log(sigma.beta) -
            (residual.deletedKnot^2 / (2 * sigma.residual)) + 
            ((crossprod(Theta.previous) - crossprod(Theta.star)) / (2 * sigma.beta)) +
            ((crossprod(LRresid.prev) - crossprod(LRresid.star)) / (2 * sigma.error))
  )
  
  if (rho > log(runif(1))) {
    return (list(X=X.star, knots=knots.star, Theta=Theta.star, accepted=TRUE))
  } else {
    return (list(X=X.previous, knots=knots.previous, Theta=Theta.previous, accepted=FALSE))         
  }
  
}


#' Add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
moveKnot.binary <- function(dist, knots.previous, knots.stepSize, knots.candidatePositions,
                            outcomes, times.dropout, times.observation, covariates,
                            X.previous, Theta.previous, 
                            Z, alpha, betaCovariates, eta.null) {
  
  #Pick a knot to move 
  knotToMove <- sample(knots.previous, 1) 
  # get index of knot to be moved
  index <- which(knots.previous == knotToMove) 
  # get the knots that are staying in the same place
  knotsToKeep <- knots.previous[-index] 
  
  # find a new location from the potential knot locations
  # TODO: floating point issue not recognizing existing knot positions
  potentialLocations <- 
    knots.candidatePositions[!(knots.candidatePositions %in% knots.previous)]
  # here we only allow movement within some small window
  potentialLocations <- 
    potentialLocations[potentialLocations > (knotToMove - knots.stepSize) & 
                         potentialLocations < (knotToMove + knots.stepSize)] 
  
  if (length(potentialLocations) > 0) {
    # pick a new location
    p <- sample(potentialLocations,1)
    # get the new knots
    knots.star <- sort(c(knotsToKeep, p))
    
    if (length(knots.star) > 1) {
      knotstar.b <- range(knots.star)
      knotstar.i <- knots.star[!(knots.star %in% knotstar.b)]
      # Calculate spline transformation of dropout time and proposed X matrix and residuals
      X.star <- cbind(
        rep(1, length(times.observation)),
        ns(times.dropout, knots=knotstar.i, Boundary.knots=knotstar.b, intercept=T) * times.observation
      )
    } else {
      X.star <- cbind(rep(1, length(times.observation)), times.observation)
    }
    
    # Calculate residuals for likelihood ratio
    y = as.matrix(outcomes)
    if (!is.null(covariates)) {
      cBeta = as.vector(as.matrix(covariates) %*% betaCovariates)
    }
    zAlpha = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
    
    # calculate the previous eta and associated probability
    eta.previous = as.vector(X.previous %*% Theta.previous + cBeta + zAlpha)
    prob.previous = inv.logit(eta.previous)
    # adjust probabilities within tolerance levels
    prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
    prob.previous[prob.previous < model.options$prob.max] <-  model.options$prob.max
    loglikelihood.previous <- sum(log((1 - prob.previous[outcomes==0]))) + sum(log(prob.previous[outcomes==1]))    
    
    # calculate the proposal eta and associated probability
    eta.star = as.vector(X.star %*% Theta.star + cBeta + zAlpha)
    prob.star = inv.logit(eta.star)
    # adjust probabilities within tolerance levels
    prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
    prob.star[prob.star < model.options$prob.max] <-  model.options$prob.max
    loglikelihood.star <- sum(log((1 - prob.star[outcomes==0]))) + sum(log(prob.star[outcomes==1]))    
    
    # Calculate the acceptance probability
    rho <- (log(length(potentialLocations)) - log(length(knots.star)) + loglikelihood.star - loglikelihood.previous)
    if (rho > log(runif(1))) {
      return (list(knots=knots.star, X=X.star, proposed=TRUE, accepted=TRUE))
    } else {
      return (list(knots=knots.previous, X=X.previous, proposed=TRUE, accepted=FALSE))
    }
  } else {
    # nowhere to move to
    return (list(knots=knots.previous, X=X.previous, proposed=FALSE, accepted=FALSE))
  }
  
}

#' Add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
moveKnot.gaussian <- function(dist, knots.previous, knots.stepSize, knots.candidatePositions,
                              outcomes, times.dropout, times.observation, covariates,
                              X.previous, Theta.previous, 
                              Z, alpha, betaCovariates, sigma.error) {
  #Pick a knot to move 
  knotToMove <- sample(knots.previous, 1) 
  # get index of knot to be moved
  index <- which(knots.previous == knotToMove) 
  # get the knots that are staying in the same place
  knotsToKeep <- knots.previous[-index] 
  
  # find a new location from the potential knot locations
  # TODO: floating point issue not recognizing existing knot positions
  potentialLocations <- 
    knots.candidatePositions[!(knots.candidatePositions %in% knots.previous)]
  # here we only allow movement within some small window
  potentialLocations <- 
    potentialLocations[potentialLocations > (knotToMove - knots.stepSize) & 
                         potentialLocations < (knotToMove + knots.stepSize)] 
  
  if (length(potentialLocations) > 0) {
    # pick a new location
    p <- sample(potentialLocations,1)
    # get the new knots
    knots.star <- sort(c(knotsToKeep, p))
    
    if (length(knots.star) > 1) {
      knotstar.b <- range(knots.star)
      knotstar.i <- knots.star[!(knots.star %in% knotstar.b)]
      # Calculate spline transformation of dropout time and proposed X matrix and residuals
      X.star <- cbind(
        rep(1, length(times.observation)),
        ns(times.dropout, knots=knotstar.i, Boundary.knots=knotstar.b, intercept=T) * times.observation
      )
    } else {
      X.star <- cbind(rep(1, length(times.observation)), times.observation)
    }
    
    # Calculate residuals for likelihood ratio
    y = as.matrix(outcomes)
    if (!is.null(covariates)) {
      cBeta = as.vector(as.matrix(covariates) %*% betaCovariates)
    }
    zAlpha = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
    
    
    if (!is.null(covariates)) {
      LRresid.star <- as.vector(y - X.star %*% Theta.previous - cBeta - zAlpha)
      LRresid.prev <- as.vector(y - X.previous %*% Theta.previous - cBeta - zAlpha)
    } else {
      LRresid.star <-as.vector(y - X.star %*% Theta.previous - zAlpha)
      LRresid.prev <- as.vector(y - X.previous %*% Theta.previous - zAlpha)
    }
    
    # Calculate the acceptance probability
    rho <- (log(length(potentialLocations)) - log(length(knots.star)) + 
              (crossprod(LRresid.prev)-crossprod(LRresid.star))/2/sigma.error)
    if (rho > log(runif(1))) {
      return (list(knots=knots.star, X=X.star, proposed=TRUE, accepted=TRUE))
    } else {
      return (list(knots=knots.previous, X=X.previous, proposed=TRUE, accepted=FALSE))
    }
  } else {
    # nowhere to move to
    return (list(knots=knots.previous, X=X.previous, proposed=FALSE, accepted=FALSE))
  }
  
}

#' Add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
updateFixedEffects.binary <- function(knots.previous, outcomes, times.dropout, times.observation, covariates,
                                      X.previous, Theta.previous, 
                                      Z, alpha, betaCovariates, 
                                      sigma.beta, lambda.numKnots, eta.null, model.options) {
  
  # build components of eta
  y <- as.matrix(outcomes)
  if (!is.null(covariates)) {
    cBeta = as.vector(as.matrix(covariates) %*% betaCovariates)
  }
  zAlpha = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
  
  XTheta.previous = X.previous %*% Theta.previous 
  # calculate the previous eta and associated probability
  eta.previous = as.vector(XTheta.previous + ifelse(!is.null(covariates), cBeta, 0) + zAlpha)
  prob.previous = inv.logit(eta.previous)
  # adjust probabilities within tolerance levels
  prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
  prob.previous[prob.previous < model.options$prob.max] <-  model.options$prob.max
  loglikelihood.previous <- sum(log((1 - prob.previous[y==0]))) + sum(log(prob.previous[y==1]))    
  
  # create the covariance matrix for the intercept and theta coefficients for the splines - R0
  covarIntTheta <- diag(rep(sigma.beta, (length(knots.previous)+1))) 
  covarIntThetaInverse <- diag(rep(1/sigma.beta, (length(knots.previous)+1))) 
  # build y-tilde
  yt = XTheta.previous + (y - prob.previous) * (1 / (prob.previous * (1 - prob.previous)))
  # build the weight matrix
  weight <- Diagonal(x = prob.previous * (1-prob.previous))
  # build the covariance of the beta coefficients and make sure it is positive definite
  covar <- as.matrix(nearPD(solve(covarIntThetaInverse + crossprod(X.previous, weight %*% X.previous)))$mat)
  # build the mean
  mu <- covar  %*% (crossprod(X.previous, weight %*% yt))
  # draw the proposed coefficients for the fixed effects
  Theta.star <- rmvnorm(1, mu, covar)
  
  # get proposal probabilities
  XTheta.star = X.previous %*% t(Theta.star)
  eta.star <- XTheta.star + ifelse(!is.null(covariates), cBeta, 0) + zAlpha
  prob.star = inv.logit(eta.previous)
  # adjust probabilities within tolerance levels
  prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
  prob.star[prob.star < model.options$prob.max] <-  model.options$prob.max
  loglikelihood.star <- sum(log((1 - prob.star[y == 0]))) + sum(log(prob.star[y == 1]))    
  
  y.star <- XTheta.star + (y - prob.star) * (1 / (prob.star * (1 - prob.star)))
  weight.star <- Diagonal(x=prob.star * (1 - prob.star))
  covar.star <- as.matrix(nearPD(solve(covarIntThetaInverse + crossprod(X.previous, weight.star %*% X.previous)))$mat)
  mu.star <- covar.star %*% (crossprod(Xt,weightstar%*%ystar))
  
  # Metropolis hastings step
  rho <- (loglikelihood.star - loglikelihood.previous + 
            log(dmvnorm(as.vector(ct[[i]]), as.vector(mntstar), covtstar)) - 
            log(dmvnorm(pbeta, as.vector(mnt),covt)) +
            (crossprod(as.vector(ct[[i]])) - tcrossprod(pbeta))/ (2 * sigma.beta))
  if (rho > log(runif(1))) {
    return (list(Theta=Theta.star, accepted=TRUE))
  } else {
    return (list(Theta=Theta.previous, accepted=FALSE))
  }
  
}


#' Add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
updateFixedEffects.gaussian <- function(dist, knots.previous, model.options, 
                                        outcomes, times.dropout, times.observation, covariates,
                                        X.previous, Theta.previous, 
                                        Z, alpha, betaCovariates,  
                                        sigma.error, sigma.beta, lambda.numKnots) {
  
  # create the covariance matrix for the intercept and theta coefficients for the splines - R0
  covarIntTheta <- diag(rep(sigma.beta, (length(knots.previous)+1))) 
  covarIntThetaInverse <- diag(rep(1/sigma.beta, (length(knots.previous)+1))) 
  # Calculate residuals 
  y = as.matrix(outcomes)
  if (!is.null(covariates)) {
    cBeta = as.vector(as.matrix(covariates) %*% betaCovariates)
  }
  zAlpha = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
  
  if (!is.null(covariates)) {
    LRresid <- as.vector(y - cBeta - zAlpha)
  } else {
    LRresid <- as.vector(y - zAlpha)
  }
  
  # calculate the covariance and mean of the proposal distribution
  proposedCovariance <- ginv(covarIntThetaInverse + (1/sigma.error) * crossprod(X.previous))
  proposedMean <- proposedCovariance %*% ((1/sigma.error) * crossprod(X.previous, LRresid))
  
  # Scale Cov to adjust acceptance rate
  adjust = ifelse(!is.null(mcmc.options$fixedEffectAcceptRateAdjust), 
                  mcmc.options$fixedEffectAcceptRateAdjust, 1)
  proposedCovariance <- adjust * proposedCovariance 
  
  # ensure the covariance is positive definite
  proposedCovariance <- as.matrix(nearPD(proposedCovariance)$mat) 
  # draw a proposed set of coefficients
  Theta.star <- t(as.matrix(rmvnorm(1, proposedMean, proposedCovariance)))
  
  # Calculate residuals for likelihood ratio
  resid.star <- LRresid - X.previous %*% Theta.star
  resid.prev <- LRresid - X.previous %*% Theta.previous
  
  # Calculate acceptance probability
  rho <- (log(dmvnorm(as.vector(Theta.previous), as.vector(proposedMean), proposedCovariance)) - 
            log(dmvnorm(as.vector(Theta.star), as.vector(proposedMean), proposedCovariance)) + 
            (crossprod(Theta.previous)-crossprod(Theta.star))/(2 * sigma.beta) + 
            (crossprod(resid.prev)-crossprod(resid.star))/(2 * sigma.error))
  
  if (rho > log(runif(1))) {
    return (list(Theta=Theta.star, accepted=TRUE))
  } else {
    return (list(Theta=Theta.previous, accepted=FALSE))
  }
}

#'
#' Update the regression coefficients related to common
#' covariates and group effects
#'
updateFixedEffectsCovariates.binary <- function(outcomes, covariates, 
                                                X, Theta, Z, alpha, betaCovariates.previous,  
                                                sigma.error, sigma.beta) {
  # make sure there are covariates to update 
  if (is.null(covariates)) {
    return (list(betaCovariates=NULL, accepted=FALSE))
  }
  
  # build components of eta
  y <- as.matrix(outcomes)
  C = as.matrix(covariates)
  cBeta = as.vector(C %*% betaCovariates.previous)
  zAlpha = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
  XTheta.previous = X.previous %*% Theta.previous 
  # calculate the previous eta and associated probability
  eta.previous = as.vector(XTheta.previous + cBeta + zAlpha)
  prob.previous = inv.logit(eta.previous)
  # adjust probabilities within tolerance levels
  prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
  prob.previous[prob.previous < model.options$prob.max] <-  model.options$prob.max
  loglikelihood.previous <- sum(log((1 - prob.previous[y==0]))) + sum(log(prob.previous[y==1]))    
  
  # create the covariance matrix for the intercept and theta coefficients for the splines - R0
  covarIntTheta <- diag(rep(sigma.beta, (length(knots.previous)+1))) 
  covarIntThetaInverse <- diag(rep(1/sigma.beta, (length(knots.previous)+1))) 
  # build y-tilde
  yt = cBeta + (y - prob.previous) * (1 / (prob.previous * (1 - prob.previous)))
  # build the weight matrix
  weight <- Diagonal(x = prob.previous * (1-prob.previous))
  # build the covariance of the beta coefficients and make sure it is positive definite
  covar <- as.matrix(nearPD(solve(covarIntThetaInverse + crossprod(C, weight %*% C)))$mat)
  # build the mean
  mu <- covar  %*% (crossprod(C, weight %*% yt))
  # draw the proposed coefficients for the fixed effects
  betaCovariates.star <- rmvnorm(1, mu, covar)
  
  # get proposal probabilities
  eta.star <- XTheta.previous + C %*% betaCovariates.star + zAlpha
  prob.star = inv.logit(eta.star)
  # adjust probabilities within tolerance levels
  prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
  prob.star[prob.star < model.options$prob.max] <-  model.options$prob.max
  loglikelihood.star <- sum(log((1 - prob.star[y == 0]))) + sum(log(prob.star[y == 1]))    
  
  y.star <- cBeta.star + (y - prob.star) * (1 / (prob.star * (1 - prob.star)))
  weight.star <- Diagonal(x=prob.star * (1 - prob.star))
  covar.star <- as.matrix(nearPD(solve(covarIntThetaInverse + crossprod(C, weight.star %*% C)))$mat)
  mu.star <- covar.star %*% (crossprod(C,weight.star%*%y.star))
  
  # Metropolis hastings step
  rho <- (loglikelihood.star - loglikelihood.previous + 
            log(dmvnorm(as.vector(betaCovariates.star), as.vector(mu.star), covar.star)) - 
            log(dmvnorm(pbeta, as.vector(mnt),covt)) +
            (crossprod(as.vector(betaCovariates.star)) - tcrossprod(prob.star))/ (2 * sigma.beta))
  
  if (rho > log(runif(1))) {
    return (list(betaCovariates=as.vector(betaCovariates.star), accepted=TRUE))
  } else {
    return (list(betaCovariates=betaCovariates.previous, accepted=FALSE))
  }
}



#'
#' Update the regression coefficients related to common
#' covariates and group effects
#'
updateFixedEffectsCovariates.gaussian <- function(outcomes, covariates, 
                                                  X, Theta, Z, alpha, betaCovariates.previous,  
                                                  sigma.error, sigma.beta) {
  
  # build X * beta for each group and combine into a single vector
  cBeta = vector()
  for(i in 1:length(X)) {
    cBeta.group = X[[i]] %*% Theta[[i]]
    cBeta <- c(cBeta, cBeta.group) 
  }
  
  # calculate the residuals
  residuals <- outcomes - cBeta - Z$intercept * alpha$intercept - Z$slope * alpha$slope
  residuals = residuals[,outcomes.var]
  # get the proposed mean/variance of the fixed effects associated with covariates
  covariates.matrix = as.matrix(covariates)
  covarIntTheta <- diag(rep(sigma.beta, ncol(covariates))) 
  covarIntThetaInverse <- diag(rep(1/sigma.beta, ncol(covariates))) 
  proposedCovariance <- ginv(covarIntThetaInverse + 
                               (1/sigma.error) * t(covariates.matrix) %*% covariates.matrix)
  proposedMean <- proposedCovariance %*% 
    ((1/sigma.error) * t(covariates.matrix) %*% as.matrix(residuals))
  
  # Scale Cov to adjust acceptance rate
  proposedCovariance <- proposedCovariance * 
    ifelse(!is.null(mcmc.options$fixedEffectAcceptRateAdjust), 
           mcmc.options$fixedEffectAcceptRateAdjust, 1)
  
  # ensure the covariance is positive definite
  proposedCovariance <- as.matrix(nearPD(proposedCovariance)$mat) 
  
  # draw a proposed set of coefficients for the covariates
  betaCovariates.star <- rmvnorm(1, proposedMean, proposedCovariance)
  
  # calculate the residuals
  resid.star <- residuals - covariates.matrix %*% t(betaCovariates.star)
  resid.previous <- residuals - covariates.matrix %*% (as.matrix(betaCovariates.previous))
  
  # calculate the acceptance probability
  rho<-sum(log(dnorm(as.vector(resid.star), 0, sqrt(sigma.error)))) + 
    log(dmvnorm(betaCovariates.star, rep(0, ncol(covariates)), covarIntTheta)) + 
    log(dmvnorm(betaCovariates.previous, proposedMean, proposedCovariance)) - 
    sum(log(dnorm(resid.star, 0, sqrt(sigma.error)))) - 
    log(dmvnorm(betaCovariates.previous, rep(0, ncol(covariates)), covarIntTheta)) - 
    log(dmvnorm(betaCovariates.star, proposedMean, proposedCovariance))
  
  if (rho > log(runif(1))) {
    return (list(betaCovariates=as.vector(betaCovariates.star), accepted=TRUE))
  } else {
    return (list(betaCovariates=betaCovariates.previous, accepted=FALSE))
  }
}


#' Update the random effects
#' 
#' @param 
#' @return 
#' @examples
#'
updateRandomEffects.binary <- function(dist, numSubjects, numObservations, firstObsPerSubject,
                                       subjectsPerGroup, ids, outcomes, times.observation, 
                                       covariates, X, Theta, 
                                       Z, alpha, betaCovariates,  
                                       sigma.randomIntercept, sigma.randomSlope,
                                       sigma.randomInterceptSlope, sigma.error) {
  
  # Update B1 and B2 using MH step
  
  l<-ls<-rep(0,nsub)   # store subject-level log-likes
  
  # Current Observation-Level Log-Likelihood
  Lt<-Lstar<-rep(0, N)
  Lt[y==0]<-log(1-pt[y==0])  
  Lt[y==1]<-log(pt[y==1])
  # Draw Candidate from Symmetric Rand Walk Proposal 
  bs<-cbind(B1,B2)+rmvt(nsub,sigma=diag(2),3)      # or b+rmvnorm(n,sigma=Sigma)
  etastar<-Xt%*%c(ct[[i]])+rep(bs[,1],nobs)+rep(bs[,2],nobs)*t  
  pstar<-inv.logit(etastar)
  pstar[pstar<0.00001]<-0.00001
  pstar[pstar>0.99999]<-0.99999
  Lstar[y==0]<-log(1-pstar[y==0])
  Lstar[y==1]<-log(pstar[y==1])
  # Subject-Level Log-Likelihood
  l<-tapply(Lt,patid,sum)
  ls<-tapply(Lstar, patid, sum)
  
  ratio<-ls+dmvnorm(bs,c(0,0),matrix((Sigma.bs[i-1,]), nrow=2, ncol=2),log=TRUE)-(l+dmvnorm(cbind(B1,B2),c(0,0),matrix((Sigma.bs[i-1,]), nrow=2, ncol=2),log=TRUE))
  
  un<-1*(log(runif(nsub))<ratio)
  B1[un==1]<-bs[un==1,1]
  B2[un==1]<-bs[un==1,2]
  etat<-Xt%*%c(ct[[i]])+rep(B1,nobs)+rep(B2,nobs)*t 
  pt<-inv.logit(etat)
  loglt<-sum(log((1-pt[y==0])))+sum(log(pt[y==1]))
}



#' Update the random effects
#' 
#' @param 
#' @return 
#' @examples
#'
updateRandomEffects.gaussian <- function(dist, numSubjects, numObservations, firstObsPerSubject,
                                         subjectsPerGroup, ids, outcomes, times.observation, 
                                         covariates, X, Theta, 
                                         Z, alpha, betaCovariates,  
                                         sigma.randomIntercept, sigma.randomSlope,
                                         sigma.randomInterceptSlope, sigma.error) {
  # get the random intercepts, one per subject
  alpha.slopeOnePerSubject = alpha$slope[firstObsPerSubject]
  alpha.interceptOnePerSubject = alpha$intercept[firstObsPerSubject]
  
  # calculate rho
  rho = (sigma.randomInterceptSlope / sqrt(sigma.randomIntercept * sigma.randomSlope))
  
  # some convenience variables
  tau.error = 1 / sigma.error
  tau.randomIntercept = 1 / sigma.randomIntercept
  tau.randomSlope = 1 / sigma.randomSlope
  tau.randomInterceptSlope = 1 / sigma.randomInterceptSlope
  
  # build the residuals
  cBeta = vector()
  for(i in 1:length(X)) {
    cBeta.group = X[[i]] %*% Theta[[i]]
    cBeta <- c(cBeta, cBeta.group) 
  }
  # calculate the residuals
  residuals <- outcomes - cBeta - Z$slope * alpha$slope
  if (!is.null(covariates)) {
    residuals <- residuals - as.matrix(covariates) %*% as.matrix(betaCovariates)
  }
  residuals = residuals[,outcomes.var]
  
  # get the conditional distribution of the random intercept, given the random slope
  variance.randomIntercept <- 1 / (tau.error * numObservations + tau.randomIntercept * 
                                     1/(1 - rho*rho))
  mean.randomIntercept <- (variance.randomIntercept * 
                             (tau.error * as.vector(tapply(residuals,ids,sum)) + 
                                alpha.slopeOnePerSubject * (rho / (1 - rho*rho)) * 
                                (1 / sqrt(sigma.randomIntercept * sigma.randomSlope))))
  # draw a new random effect sample
  randomIntercepts <-rnorm(numSubjects, mean.randomIntercept, sqrt(variance.randomIntercept))
  
  # get the conditional distribution of the random slope -- FIX RESIDUALS!
  # update the residuals
  residuals = outcomes - cBeta - Z$intercept * alpha$intercept 
  if (!is.null(covariates)) {
    residuals <- residuals - as.matrix(covariates) %*% as.matrix(betaCovariates)
  }
  residuals = residuals[,outcomes.var]
  
  variance.randomSlope <- 1 / (tau.error * as.vector(tapply(times.observation^2, ids, sum)) + 
                                 tau.randomSlope * 1/(1 - rho*rho))
  mean.randomSlope <- (variance.randomSlope * 
                         (tau.error * as.vector(tapply(times.observation*residuals,ids,sum)) + 
                            alpha.interceptOnePerSubject * (rho / (1 - rho*rho)) * 
                            (1 / sqrt(sigma.randomIntercept * sigma.randomSlope))))
  
  
  randomSlopes <-rnorm(numSubjects, mean.randomSlope, sqrt(variance.randomSlope))
  
  # build the new random effects matrix
  alpha.total = data.frame(intercept=randomIntercepts, slope=randomSlopes)
  
  alpha = alpha.total[rep(seq_len(nrow(alpha.total)), numObservations),]
  # expand back out to complete alpha matrix
  return (alpha)
  
}

#' Add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
#' 
updateCovarianceParameters.binary <- function(numSubjects, firstObsPerSubject,
                                              alpha, model.options) {
  
  
  # sample from an inverse wishart to update the covariance of the random effects
  perSubjectAlpha = alpha[firstObsPerSubject,]
  sigma.alpha <- riwish(prior.options$sigmaError.df + numSubjects, 
                        prior.options$sigmaError.scaleMatrix + crossprod(as.matrix(perSubjectAlpha)))
  
  sigma.randomIntercept = sigma.alpha[1,1]
  sigma.randomSlope = sigma.alpha[2,2]
  sigma.randomInterceptSlope = sigma.alpha[1,2]
  
  return (list(sigma.randomIntercept = sigma.alpha[1,1],
               sigma.randomSlope = sigma.alpha[2,2],
               sigma.randomInterceptSlope = sigma.alpha[1,2]))
}

#' Add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
#' 
updateCovarianceParameters.gaussian <- function(totalObservations, numSubjects, 
                                                firstObsPerSubject,
                                                outcomes, covariates, X, Theta, 
                                                Z, alpha, betaCovariates,
                                                sigma.error,
                                                prior.options) {
  
  # sample from an inverse wishart to update the covariance of the random effects
  perSubjectAlpha = alpha[firstObsPerSubject,]
  sigma.alpha <- riwish(prior.options$sigmaError.df + numSubjects, 
                        prior.options$sigmaError.scaleMatrix + crossprod(as.matrix(perSubjectAlpha)))
  
  sigma.randomIntercept = sigma.alpha[1,1]
  sigma.randomSlope = sigma.alpha[2,2]
  sigma.randomInterceptSlope = sigma.alpha[1,2]
  
  # update sigma.error
  # build the residuals
  cBeta = vector()
  for(i in 1:length(X)) {
    cBeta.group = X[[i]] %*% Theta[[i]]
    cBeta <- c(cBeta, cBeta.group) 
  }
  # calculate the residuals
  residuals <- outcomes - cBeta - Z$intercept * alpha$intercept - Z$slope * alpha$slope
  if (!is.null(covariates)) {
    residuals <- residuals - as.matrix(covariates) %*% as.matrix(betaCovariates)
  }
  residuals = residuals[,outcomes.var]
  
  # sample from an inverse Gamma to update sigma.error
  shape <- prior.options$shape.tau + (totalObservations / 2) 
  rate <- prior.options$rate.tau + (crossprod(residuals) / 2)
  sigma.error <- 1 / rgamma(1, shape, rate)
  
  return (list(sigma.error = sigma.error,
               sigma.randomIntercept = sigma.alpha[1,1],
               sigma.randomSlope = sigma.alpha[2,2],
               sigma.randomInterceptSlope = sigma.alpha[1,2]))
  
}

#'
#' Calculate the marginal slope at the given iteration
#'
#'
calculateMarginalSlope <- function(knotsByGroup, ThetaByGroup, subjectsPerGroup,
                                   times.dropout) {
  
  marginal <- vector(0, mode="list")
  startRow <- 1
  for(i in 1:length(knotsByGroup)) {
    knots <- knotsByGroup[[i]]
    times.dropout.group = times.dropout[startRow:(startRow + subjectsPerGroup[[i]] - 1)]
    # get the current spline coefficients minus the intercept
    Theta <- ThetaByGroup[[i]]
    ThetaNoInt <- Theta[2:length(Theta)] 
    
    if (length(knots) > 1) {
      knots.boundary = range(knots)
      knots.interior = knots[-c(1,length(knots))] 
      # Calculate spline transformation of dropout time and create the proposed X matrix
      spline <- ns(times.dropout.group, knots=knots.interior, Boundary.knots=knots.boundary,
                   intercept=T) 
      
      # randomly select the proportion of subjects dropping out at each time
      propDroppedOut <- rdirichlet(1, rep(1,subjectsPerGroup[[i]]))
      # calculate the marginal slope for the current group
      marginal[[i]] <- sum(propDroppedOut * t((spline) %*% ThetaNoInt))
      
    } else {
      marginal[[i]] <- ThetaNoInt
    }
    
    startRow = startRow + subjectsPerGroup[[i]]
  }
  
  return(marginal)
}


calculateDropoutTimeSpecificSlope <- function(dropoutEstimationTimes, knotsByGroup, ThetaByGroup, 
                                              subjectsPerGroup,
                                              times.observation, times.dropout) {
  
  dropoutSpecificSlopes <- vector(0, mode="list")
  for(i in 1:length(knotsByGroup)) {
    knots <- knotsByGroup[[i]]
    
    # get the current spline coefficients minus the intercept
    Theta <- ThetaByGroup[[i]]
    ThetaNoInt <- Theta[2:length(Theta)] 
    
    if (length(knots) > 1) {
      knots.boundary = range(knots)
      knots.interior = knots[-c(1,length(knots))] 
      # Calculate spline transformation at specified dropout times
      spline <- ns(dropoutEstimationTimes, knots=knots.interior, Boundary.knots=knots.boundary, intercept=T)
      # calculate the marginal slope for the current group
      dropoutSpecificSlopes[[i]] <- t((spline) %*% ThetaNoInt)
      
    } else {
      dropoutSpecificSlopes[[i]] <- rep(ThetaNoInt, length(dropoutEstimationTimes))
    } 
  }
  
  return(dropoutSpecificSlopes)
}

#'
#'
#'
#'
#'
getInitialEstimatesTheta <- function(dist, groupList, X, outcomes) {
  data.theta = cbind(outcomes, X)
  formula = as.formula(paste(c(paste(names(outcomes), "~"), 
                               paste(names(X)[2:length(names(X))], 
                                     collapse=" + ")), 
                             collapse=" "))
  if (dist == 'gaussian') {
    fit.Theta <- lm(formula, data=data.theta)
    return (lapply(1:length(groupList), function(i) {
      return (as.vector(coef(fit.Theta)))
    }))
  } else {
    family = ifelse(dist == 'poisson', poisson(link = "log"), binomial(link = "logit"))
    fit.Theta <- glm(formula, family=family, data=data.theta)
    return (lapply(1:length(groupList), function(i) {
      return (as.vector(coef(fit.Theta)))
    }))
  }
}

#'
#'
#'
#'
#'
getInitialEstimatesCovariates <- function(dist, covariates, outcomes) {
  if (is.null(covariates) || ncol(covariates) == 0) {
    return (NULL)
  }
  
  data.covar = cbind(outcomes, covariates)
  formula = as.formula(paste(c(paste(names(outcomes), "~"), 
                               paste(names(covariates), collapse=" + ")), collapse=" "))
  if (dist == 'gaussian') {
    fit.beta <- lm(formula, data=data.covar)
    return (as.vector(coef(fit.beta))[-1])
  } else if (dist == 'binary') {
    fit.beta <- glm(formula, family=binomial, data=data.covar)
    return (as.vector(coef(fit.beta))[-1])
  } else {
    stop("unsupported distribution")
  }
}

#'
#' Fit a varying coefficient model for longitudinal studies with
#' informative dropout.  Uses a Bayesian approach with a spline fit
#' to model the relationship between dropout times and slope
#'
#' 
#'
#'
#'
informativeDropout.bayes.splines <- function(data, ids.var, outcomes.var, groups.var, 
                                             covariates.var, times.dropout.var, times.observation.var, 
                                             dist, model.options) {
  
  
  # create some reasonable defaults for the start positions and candidate positions
  # if not specified
  if (is.null(model.options$startPositions) || is.null(model.options$candidatePositions)) {
    dropout.min = min(data[,times.dropout])
    dropout.max = max(data[,times.dropout])
    if (is.null(model.options$candidatePositions)) {
      model.options$candidatePositions = seq(dropout.min, dropout.max, 1)
    }
    if (is.null(model.options$startPositions)) {
      model.options$startPositions = sample(model.options$candidatePositions, model.options$min)
    }
  }
  
  # sort the data by group, id, observation time
  data <- data[order(data[,groups.var], data[,ids.var], data[,times.observation.var]),]
  rownames(data) <- NULL
  # get the list of treatment groups
  groupList <- unique(data[,groups.var])
  # number of subjects
  numSubjects = length(unique(data[,ids.var]))
  # number of subjects per group
  subjectsPerGroup = sapply(groupList, function(group) { 
    # get the covariates
    return(length(unique(data[data[,groups.var] == group, ids.var])))
  })
  # start/end rows for groups
  numObsPerGroup = sapply(groupList, function(group) {
    return(nrow(data[data[,groups.var] == group,]))
  })
  startRowByGroup = c(1, 1 + cumsum(numObsPerGroup)[-length(numObsPerGroup)])
  endRowByGroup = cumsum(numObsPerGroup)
  
  # number of observations per subject
  numObservations = sapply(unique(data[,ids.var]), function(id) {
    return(nrow(data[data[,ids.var] == id,]))
  })
  # index of the first observation per subject
  firstObsPerSubject = c(1,1+cumsum(numObservations)[-length(numObservations)])
  
  # initialize the random effects design matrix
  # note, this is not a true full 
  Z = data.frame(groups=data[,groups.var], intercept=rep(1,nrow(data)), times=data[,times.observation.var])
  names(Z) = c(groups.var, "intercept", "slope")
  # initialize the random effects
  alpha = data.frame(intercept=rep(0, nrow(data)), slope=rep(0, nrow(data)))
  names(alpha) = c("intercept", "slope")
  
  # initialize the X matrix - this is split by group since each group may
  X = lapply(groupList, function(group) { 
    groupTimes = data[data[,groups.var] == group,times.observation.var]
    groupDropout = data[data[,groups.var] == group,times.dropout.var]
    knots.boundary = range(model.options$startPositions)
    knots.interior = model.options$startPositions[-c(1,length(model.options$startPositions))] 
    return(as.matrix(cbind(
      rep(1,length(groupTimes)),
      ns(groupDropout, knots=knots.interior, Boundary.knots=knots.boundary, intercept=T) * groupTimes
    )))
  })
  # combine into an initial complete X matrix
  X.full = as.data.frame(abind(X, along=1))
  names(X.full) <- sapply(0:(ncol(X.full)-1), function(i) { return(paste("theta", i, sep=''))})
  
  # extract the outcomes
  outcomes = as.data.frame(data[,outcomes.var])
  names(outcomes) = c(outcomes.var)
  # extract the covariates
  covariates = NULL
  if (!is.null(covariates.var)) {
    covariates = as.data.frame(data[,covariates.var])
    names(covariates) = c(covariates.var)
  }
  
  # use a simple linear model fit to obtain initial estimates for 
  # the spline coefficients
  Theta.init = getInitialEstimatesTheta(dist, groupList, X.full, outcomes)
  # use a simple linear model fit to obtain initial estimates for the
  # covariate coefficients
  betaCovariate.init = getInitialEstimatesCovariates(dist, covariates, outcomes)
  
  # for binary case, initialize eta
  if (dist == 'binary') {
    if (is.null(model.options$eta.null)) {
      eta.null = model.options$eta.null
    } else {
      eta.null = lapply(subjectsPerGroup, function(N) {
        return(logit(rep(sum(y)/N, N)))
      })
    }
  }
  
  
  # initialize the first model iteration
  # TODO: mcmc options specify starting values
  modelIterationList <- vector(mode = "list", length = mcmc.options$iterations)
  modelIterationList[[1]] = 
    rjmcmc.iteration(
      knots=lapply(1:length(groupList), function(i) { return (model.options$startPositions); }), 
      Theta=Theta.init, betaCovariate=betaCovariate.init,
      sigma.error = lapply(1:length(groupList), function(i) { return (model.options$sigma.error); }),
      sigma.spline = lapply(1:length(groupList), function(i) { return (model.options$sigma.spline); }),
      sigma.randomIntercept = lapply(1:length(groupList), function(i) { return (model.options$sigma.randomIntercept); }),
      sigma.randomSlope = lapply(1:length(groupList), function(i) { return (model.options$sigma.randomSlope); }),
      sigma.randomInterceptSlope = lapply(1:length(groupList), function(i) { return (model.options$sigma.randomInterceptSlope); }),
      shape.tau = lapply(1:length(groupList), function(i) { return (model.options$shape.tau); }),
      rate.tau = lapply(1:length(groupList), function(i) { return (model.options$rate.tau); })
    )
  #
  # Run the reversible jump MCMC
  #
  for (i in 2:mcmc.options$iterations) {
    print(paste("ITER = ", i, sep=""))
    
    model.previous = modelIterationList[[i-1]]
    # make a copy which will be modified as we move through the iteration
    model.current = model.previous
    
    for (group.index in 1:length(groupList)) {
      group = groupList[group.index]
      # get the subset of data for this group
      groupData = data[data[,groups.var] == group,]
      group.startRow = startRowByGroup[group.index]
      group.endRow = endRowByGroup[group.index]
      # group specific variables (these are unchanged throughout the iteration)
      group.outcomes = groupData[,outcomes.var]
      group.times.dropout = groupData[,times.dropout.var]
      group.times.observation = groupData[,times.observation.var]
      
      group.covariates = NULL
      if (!is.null(covariates)) {
        group.covariates = groupData[, covariates.var]
      }
      
      group.Z = Z[Z[,groups.var] == group,c("intercept", "slope")]
      group.alpha = alpha[(group.startRow:group.endRow),]
      
      # randomly decide to add/remove a knot
      u = runif(1)
      if ((u < model.options$birthProbability && 
           length(model.current$knots[[group.index]]) < model.options$max) || 
          (length(model.current$knots[[group.index]]) <= model.options$min)) {
        if (dist == "gaussian") {
          # add a knot
          result = addKnot.gaussian(knots.previous=model.current$knots[[group.index]], 
                                    model.options=model.options, 
                                    outcomes=group.outcomes, 
                                    times.dropout=group.times.dropout, 
                                    times.observation=group.times.observation, 
                                    covariates=group.covariates,
                                    X.previous=X[[group.index]], 
                                    Theta.previous=model.current$Theta[[group.index]],
                                    Z=group.Z, alpha=group.alpha, 
                                    betaCovariates=model.current$betaCovariates,  
                                    sigma.residual=mcmc.options$sigma.residual,
                                    sigma.error=model.current$sigma.error, 
                                    sigma.beta=prior.options$sigma.beta, 
                                    lambda.numKnots=prior.options$lambda.numKnots)
        } else {
          # binary case
          result = addKnot.binary(knots.previous=model.current$knots[[group.index]], 
                                  model.options=model.options, 
                                  outcomes=group.outcomes, 
                                  times.dropout=group.times.dropout, 
                                  times.observation=group.times.observation, 
                                  covariates=group.covariates,
                                  X.previous=X[[group.index]], 
                                  Theta.previous=model.current$Theta[[group.index]],
                                  Z=group.Z, alpha=group.alpha, 
                                  betaCovariates=model.current$betaCovariates,  
                                  sigma.residual=mcmc.options$sigma.residual,
                                  sigma.error=model.current$sigma.error, 
                                  sigma.beta=prior.options$sigma.beta, 
                                  lambda.numKnots=prior.options$lambda.numKnots)
        }
        
        # update the model iteration
        X[[group.index]] = result$X
        model.current$knots[[group.index]] = result$knots
        model.current$Theta[[group.index]] = result$Theta
        model.current$proposed$knot.add = TRUE
        model.current$accepted$knot.add = result$accepted
        
      } else {
        if (dist == "gaussian") {
          
          # remove a knot
          result = removeKnot.gaussian(dist=dist, knots.previous=model.current$knots[[group.index]], 
                                       model.options=model.options, 
                                       outcomes=group.outcomes, 
                                       times.dropout=group.times.dropout, 
                                       times.observation=group.times.observation, 
                                       covariates=group.covariates,
                                       X.previous=X[[group.index]], 
                                       Theta.previous=model.current$Theta[[group.index]],
                                       Z=group.Z, alpha=group.alpha, 
                                       betaCovariates=model.current$betaCovariates,  
                                       sigma.residual=mcmc.options$sigma.residual,
                                       sigma.error=model.current$sigma.error, 
                                       sigma.beta=prior.options$sigma.beta, 
                                       lambda.numKnots=prior.options$lambda.numKnots)  
        } else {
          # remove a knot
          result = removeKnot.binary(dist=dist, knots.previous=model.current$knots[[group.index]], 
                                     model.options=model.options, 
                                     outcomes=group.outcomes, 
                                     times.dropout=group.times.dropout, 
                                     times.observation=group.times.observation, 
                                     covariates=group.covariates,
                                     X.previous=X[[group.index]], 
                                     Theta.previous=model.current$Theta[[group.index]],
                                     Z=group.Z, alpha=group.alpha, 
                                     betaCovariates=model.current$betaCovariates,  
                                     sigma.residual=mcmc.options$sigma.residual,
                                     sigma.error=model.current$sigma.error, 
                                     sigma.beta=prior.options$sigma.beta, 
                                     lambda.numKnots=prior.options$lambda.numKnots)  
        }
        # update the model iteration
        X[[group.index]] = result$X
        model.current$knots[[group.index]] = result$knots
        model.current$Theta[[group.index]] = result$Theta
        model.current$proposed$knot.remove = TRUE
        model.current$accepted$knot.remove = result$accepted
      }
      print ("TOTAL KNOTS")
      print (length(model.current$knots[[group.index]]))
      # Move knots
      if (dist == "gaussian") {
        
        
        result = moveKnot.gaussian(dist=dist, knots.previous=model.current$knots[[group.index]], 
                                   knots.stepSize=model.options$stepSize, 
                                   knots.candidatePositions=model.options$candidatePositions,
                                   outcomes=group.outcomes, 
                                   times.dropout=group.times.dropout, 
                                   times.observation=group.times.observation, 
                                   covariates=group.covariates,
                                   X.previous=X[[group.index]], 
                                   Theta.previous=model.current$Theta[[group.index]],
                                   Z=group.Z, alpha=group.alpha, 
                                   betaCovariates=model.current$betaCovariates,  
                                   sigma.error=model.current$sigma.error)
      } else {
        result = moveKnot.binary(dist=dist, knots.previous=model.current$knots[[group.index]], 
                                 knots.stepSize=model.options$stepSize, 
                                 knots.candidatePositions=model.options$candidatePositions,
                                 outcomes=group.outcomes, 
                                 times.dropout=group.times.dropout, 
                                 times.observation=group.times.observation, 
                                 covariates=group.covariates,
                                 X.previous=X[[group.index]], 
                                 Theta.previous=model.current$Theta[[group.index]],
                                 Z=group.Z, alpha=group.alpha, 
                                 betaCovariates=model.current$betaCovariates,  
                                 sigma.error=model.current$sigma.error)
      }
      X[[group.index]] = result$X
      model.current$knots[[group.index]] = result$knots
      model.current$proposed$knot.move = result$proposed
      model.current$accepted$knot.move = result$accepted
      
      # update fixed effects (includes coefficients for covariates and time varying slopes)
      print("FIXED")
      if (dist == "gaussian") {
        result = updateFixedEffects.gaussian(dist=dist,
                                             knots.previous=model.current$knots[[group.index]], 
                                             model.options=model.options, 
                                             outcomes=group.outcomes, 
                                             times.dropout=group.times.dropout, 
                                             times.observation=group.times.observation, 
                                             covariates=group.covariates,
                                             X.previous=X[[group.index]], 
                                             Theta.previous=model.current$Theta[[group.index]],
                                             Z=group.Z, 
                                             alpha=group.alpha, 
                                             betaCovariates=model.current$betaCovariates,  
                                             sigma.error=model.current$sigma.error, 
                                             sigma.beta=prior.options$sigma.beta, 
                                             lambda.numKnots=prior.options$lambda.numKnots)
      } else {
        result = updateFixedEffects.binary(dist=dist,
                                           knots.previous=model.current$knots[[group.index]], 
                                           model.options=model.options, 
                                           outcomes=group.outcomes, 
                                           times.dropout=group.times.dropout, 
                                           times.observation=group.times.observation, 
                                           covariates=group.covariates,
                                           X.previous=X[[group.index]], 
                                           Theta.previous=model.current$Theta[[group.index]],
                                           Z=group.Z, 
                                           alpha=group.alpha, 
                                           betaCovariates=model.current$betaCovariates,  
                                           sigma.error=model.current$sigma.error, 
                                           sigma.beta=prior.options$sigma.beta, 
                                           lambda.numKnots=prior.options$lambda.numKnots)
      }
      model.current$Theta[[group.index]] = result$Theta
      model.current$proposed$fixedEffects = TRUE
      model.current$accepted$fixedEffects = result$accepted
      
    }  
    
    # update fixed effects associated with covariates
    if (!is.null(covariates)) {
      if (dist == "gaussian") {
        result = updateFixedEffectsCovariates.gaussian(dist=dist, 
                                                       outcomes=outcomes, 
                                                       covariates=covariates, 
                                                       X=X, 
                                                       Theta=model.current$Theta, 
                                                       Z=Z, alpha=alpha, 
                                                       betaCovariates.previous=model.current$betaCovariates,  
                                                       sigma.error=model.current$sigma.error, 
                                                       sigma.beta=prior.options$sigma.beta)
      } else {
        result = updateFixedEffectsCovariates.binary(dist=dist, 
                                                     outcomes=outcomes, 
                                                     covariates=covariates, 
                                                     X=X, 
                                                     Theta=model.current$Theta, 
                                                     Z=Z, alpha=alpha, 
                                                     betaCovariates.previous=model.current$betaCovariates,  
                                                     sigma.error=model.current$sigma.error, 
                                                     sigma.beta=prior.options$sigma.beta)
      }
      model.current$betaCovariates = result$betaCovariates
      model.current$proposed$fixedEffectsCovariates = TRUE
      model.current$accepted$fixedEffectsCovariates = result$accepted
    }
    
    # update random effects
    if (dist == "gaussian") {
      alpha = 
        updateRandomEffects.gaussian(dist=dist, 
                                     numSubjects=numSubjects, 
                                     numObservations=numObservations, 
                                     firstObsPerSubject=firstObsPerSubject,
                                     subjectsPerGroup=subjectsPerGroup,
                                     ids=data[,ids.var], outcomes=outcomes, 
                                     times.observation=data[,times.observation.var],
                                     covariates=covariates, 
                                     X=X, Theta=model.current$Theta, 
                                     Z=Z, alpha=alpha, 
                                     betaCovariates=model.current$betaCovariates,
                                     sigma.randomIntercept=model.current$sigma.randomIntercept, 
                                     sigma.randomSlope=model.current$sigma.randomSlope,
                                     sigma.randomInterceptSlope=model.current$sigma.randomInterceptSlope,
                                     sigma.error=model.current$sigma.error)
    } else {
      alpha = 
        updateRandomEffects.binary(dist=dist, 
                                   numSubjects=numSubjects, 
                                   numObservations=numObservations, 
                                   firstObsPerSubject=firstObsPerSubject,
                                   subjectsPerGroup=subjectsPerGroup,
                                   ids=data[,ids.var], outcomes=outcomes, 
                                   times.observation=data[,times.observation.var],
                                   covariates=covariates, 
                                   X=X, Theta=model.current$Theta, 
                                   Z=Z, alpha=alpha, 
                                   betaCovariates=model.current$betaCovariates,
                                   sigma.randomIntercept=model.current$sigma.randomIntercept, 
                                   sigma.randomSlope=model.current$sigma.randomSlope,
                                   sigma.randomInterceptSlope=model.current$sigma.randomInterceptSlope,
                                   sigma.error=model.current$sigma.error)
    }
    
    # update variance components
    if (dist == "gaussian") {
      result = updateCovarianceParameters.gaussian(dist=dist,
                                                   totalObservations=nrow(data), 
                                                   numSubjects=numSubjects,
                                                   firstObsPerSubject=firstObsPerSubject,
                                                   outcomes=outcomes, 
                                                   covariates=covariates,
                                                   X=X, Theta=model.current$Theta,
                                                   Z=Z, alpha=alpha,
                                                   betaCovariates=model.current$betaCovariates,
                                                   sigma.error=sigma.error,
                                                   prior.options=prior.options)
    } else {
      result = updateCovarianceParameters.binary(dist=dist,
                                                 totalObservations=nrow(data), 
                                                 numSubjects=numSubjects,
                                                 firstObsPerSubject=firstObsPerSubject,
                                                 outcomes=outcomes, 
                                                 covariates=covariates,
                                                 X=X, Theta=model.current$Theta,
                                                 Z=Z, alpha=alpha,
                                                 betaCovariates=model.current$betaCovariates,
                                                 sigma.error=sigma.error,
                                                 prior.options=prior.options)
    }
    model.current$sigma.error = result$sigma.error
    model.current$sigma.randomIntercept = result$sigma.randomIntercept
    model.current$sigma.randomSlope = result$sigma.randomSlope
    model.current$sigma.randomInterceptSlope = result$sigma.randomInterceptSlope
    
    # calculate marginal slope
    result = calculateMarginalSlope(knotsByGroup = model.current$knots, 
                                    ThetaByGroup = model.current$Theta, 
                                    subjectsPerGroup=subjectsPerGroup,
                                    times.dropout=data[firstObsPerSubject,times.dropout.var])
    model.current$slope.marginal = result
    
    # calculate dropout time specific slopes
    result = calculateDropoutTimeSpecificSlope(
      dropoutEstimationTimes=dropoutEstimationTimes, 
      knotsByGroup = model.current$knots, 
      ThetaByGroup = model.current$Theta, 
      subjectsPerGroup=subjectsPerGroup,
      times.observation=data[,times.observation.var],
      times.dropout=data[firstObsPerSubject,times.dropout.var]
    )
    model.current$slope.dropoutSpecific = result
    
    # save the current iteration
    modelIterationList[[i]] = model.current
  }
  
  # calculate the final estimates as the mean across the different iterations
  
  
  # return the estimates, with distributions, and the model results from each iteration
  return (modelIterationList)
  
}


