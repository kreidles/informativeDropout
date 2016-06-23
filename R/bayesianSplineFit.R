#
# Supporting functions for the Spline model / RJMCMC loop
#
#
#
#' @include generics.R

#' Data stored with each iteration of the Bayesian spline model
#' during the RJMCMC run
#'
#' @title bayes.splines.iteration
#' @param knots list of knot positions by group
#' @param Theta list of spline coefficients by group
#' @param betas.covariates regression coefficients related to covariates
#' @param sigma.residual covariance of the spline coefficients
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
#' @export bayes.splines.iteration
#'
bayes.splines.iteration <- function(knots=NULL, Theta=NULL, betas.covariates=NULL,
                                    sigma.residual = 1,
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
    
    # variance of residuals
    sigma.residual = sigma.residual,
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


#' bayes.splines.model.options
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
#' @export bayes.splines.model.options
#' 
bayes.splines.model.options = function(iterations=10000, burnin=500, thin=1, print=1,
                                       knots.prob.birth=0.5, knots.min=1, knots.max=NA, knots.stepSize=3,
                                       knots.positions.start=NULL, knots.positions.candidate=NULL,
                                       prob.min=0.00001, prob.max=0.99999, accept.rate.adjust=1,
                                       dropout.estimationTimes=NULL,
                                       sigma.beta=NULL, sigma.residual=NULL,
                                       sigma.error=NULL,
                                       sigma.error.shape.tau=NULL, sigma.error.rate.tau=NULL,
                                       lambda.numKnots=NULL,
                                       sigma.randomIntercept = NULL,
                                       sigma.randomSlope = NULL,
                                       sigma.randomInterceptSlope = NULL,
                                       sigma.randomEffects.df = NULL,
                                       sigma.randomEffects.scale = NULL,
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
  if (is.na(sigma.error.shape.tau) || is.null(sigma.error.shape.tau) || sigma.error.shape.tau <= 0) {
    stop("Prior options error :: shape.tau must be greater than 0")
  }
  if (is.na(sigma.error.rate.tau) || is.null(sigma.error.rate.tau) || sigma.error.rate.tau <= 0) {
    stop("Prior options error :: rate.tau must be greater than 0")
  }
  if (is.na(sigma.randomEffects.df) || is.null(sigma.randomEffects.df) || sigma.randomEffects.df <= 0) {
    stop("Prior options error :: sigmaError.df must be greater than 0")
  }
  if (is.null(sigma.randomEffects.scale) || 
      !is.positive.definite(sigma.randomEffects.scale)) {
    stop("Prior options error :: sigma.randomEffects.scale must be a positive definite matrix")
  }
  
  opts = list(
    # mcmc iterations
    iterations = iterations,
    # burn in period
    burnin = burnin,
    # thinning interval
    thin = thin,
    # printing interval
    print = print,
    # probability of adding a knot on a given iteration
    knots.prob.birth = knots.prob.birth,
    # minimum number of knots
    knots.min = knots.min,
    # maximum number of knots,
    knots.max = knots.max,
    # step size for moving a knot 
    knots.stepSize = knots.stepSize,
    # starting positions for the knots
    knots.positions.start = knots.positions.start,
    # candidate positions for the knots
    knots.positions.candidate = knots.positions.candidate,
    # minimum/maximum probability for Metropolis Hastings for binary outcomes
    # these control numeric instability with taking logs
    prob.min = prob.min,
    prob.max = prob.max,
    # acceptance rate adjustment
    accept.rate.adjust = accept.rate.adjust,
    # times at which the dropout time dependent slopes will be estimated
    dropout.estimationTimes = dropout.estimationTimes,
    
    # covariance of regression coefficients associated with fixed effects
    sigma.beta = sigma.beta,
    # covariance of residuals
    sigma.residual = sigma.residual,
    # for gaussian outcomes, the residual error starting value
    sigma.error = sigma.error,
    # hyperpriors for inverse gamma distribution of sigma.error
    sigma.error.shape.tau = sigma.error.shape.tau, 
    sigma.error.rate.tau = sigma.error.rate.tau,
    # prior for the poisson distribution of the number of knots
    lambda.numKnots = lambda.numKnots,
    # starting values for random effects covariance parameters
    sigma.randomIntercept = sigma.randomIntercept,
    sigma.randomSlope = sigma.randomSlope,
    sigma.randomInterceptSlope = sigma.randomInterceptSlope,
    # priors for the inverse wishart distribution of the covariance of random effects
    sigma.randomEffects.df = sigma.randomEffects.df,
    sigma.randomEffects.scale = sigma.randomEffects.scale,
    # for binary outcomes, eta null
    eta.null = eta.null
    
  )
  
  class(opts) <- append(class(opts), "bayes.splines.model.options")
  return(opts)
  
}

#' bayes.splines.fit
#' 
#' Model fit for a Bayesian spline model run
#' 
#' @param model.options the original model options
#' @param dist the distribution of the outcome ("gaussian" or "binary")
#' @param groups the list of groups
#' @param covariates.var list of covariate column names
#' @param iterations the model run iterations after removing burn-in iterations and thinning
#'
#' @export bayes.splines.fit
#' 
bayes.splines.fit <- function(model.options, dist, groups, covariates.var, iterations) {
  fit = list(
    
    # store the model options for reference
    model.options = model.options,
    # distribution of the outcome    
    dist = dist,
    # list of group names
    groups = groups,
    # list of covariate names
    covariates.var = covariates.var,
    # list of dirichlet.iteration objects
    iterations = iterations
  )
  
  class(fit) <- append(class(fit), "bayes.splines.fit")
  return(fit)
}

#' addKnot.binary
#' 
#' For binary outcomes, add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model.
#' 
#' @param model.options model options
#' @param knots.previous previous set of knots
#' @param outcomes vector of outcomes
#' @param times.dropout vector of dropout times
#' @param times.observation vector of observation times 
#' @param covariates data frame of covariates
#' @param X.previous previous X matrix
#' @param Theta.previous previous Theta values
#' @param Z random effects design matrix
#' @param alpha random effects
#' @param betaCovariates regression coefficients for covariates
#' 
#' @return list containing updated X, knots, Theta, and boolean indicating if change accepted
#' 
addKnot.binary <- function(model.options, knots.previous, outcomes, 
                           times.dropout, times.observation, 
                           covariates, X.previous, Theta.previous,
                           Z, alpha, betaCovariates) {
  
  # get the relevant priors from the model options
  sigma.residual <- model.options$sigma.residual
  sigma.beta <- model.options$sigma.beta
  lambda.numKnots <- model.options$lambda.numKnots
  eta.null <- model.options$eta.null
  
  # add a knot by randomly selecting a candidate knot
  candidatesPositions = model.options$knots.positions.candidate[! model.options$knots.positions.candidate %in% knots.previous]
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
  eta.wls <- eta.null + zAlpha + ifelse(!is.null(covariates), cBeta, 0)
  Theta.LSXprev <- wls.binary(y, X.previous, eta.wls, model.options)
  Theta.LSXstar <- wls.binary(y, X.star, eta.wls, model.options)
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
  probBirth <- ifelse(length(knots.previous) == model.options$knots.min, 1, model.options$knots.prob.birth)
  probDeath <- ifelse(length(knots.previous) == model.options$knots.max - 1, 1, 1 - model.options$knots.prob.birth)
  
  # calculate the previous eta and associated probability
  eta.previous = as.vector(X.previous %*% Theta.previous + zAlpha + ifelse(!is.null(covariates), cBeta, 0))
  prob.previous = inv.logit(eta.previous)
  # adjust probabilities within tolerance levels
  prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
  prob.previous[prob.previous > model.options$prob.max] <-  model.options$prob.max
  loglikelihood.previous <- sum(log((1 - prob.previous[outcomes==0]))) + sum(log(prob.previous[outcomes==1]))    
  
  
  # calculate the proposal eta and associated probability
  eta.star = as.vector(X.star %*% Theta.star + zAlpha + ifelse(!is.null(covariates), cBeta, 0))
  prob.star = inv.logit(eta.star)
  # adjust probabilities within tolerance levels
  prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
  prob.star[prob.star > model.options$prob.max] <-  model.options$prob.max
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


#' addKnot.gaussian
#' 
#' Add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param model.options model options
#' @param knots.previous previous set of knots
#' @param outcomes vector of outcomes
#' @param times.dropout vector of dropout times
#' @param times.observation vector of observation times 
#' @param covariates data frame of covariates
#' @param X.previous previous X matrix
#' @param Theta.previous previous Theta values
#' @param Z random effects design matrix
#' @param alpha random effects
#' @param betaCovariates regression coefficients for covariates
#' @param sigma.error residual variance
#' 
#' @importFrom MASS ginv
#' @return list containing updated X, knots, Theta, and boolean indicating if change accepted
#' 
addKnot.gaussian <- function(model.options, knots.previous, outcomes, 
                             times.dropout, times.observation, 
                             covariates, X.previous, Theta.previous,
                             Z, alpha, betaCovariates, sigma.error) {
  
  # get the relevant priors from the model options
  sigma.residual <- model.options$sigma.residual
  sigma.beta <- model.options$sigma.beta
  lambda.numKnots <- model.options$lambda.numKnots
  
  # add a knot by randomly selecting a candidate knot
  candidatePositions = model.options$knots.positions.candidate[!(model.options$knots.positions.candidate %in% knots.previous)]
  newKnot.value = sample(candidatePositions, 1)
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
  probBirth <- ifelse(length(knots.previous) == model.options$knots.min, 1, model.options$knots.prob.birth)
  probDeath <- ifelse(length(knots.previous) == model.options$knots.max - 1, 1, 1 - model.options$knots.prob.birth)
  
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
removeKnot.binary <- function(model.options, knots.previous, outcomes, times.dropout, 
                              times.observation, covariates,
                              X.previous, Theta.previous, Z, alpha, betaCovariates) {
  
  # get the relevant priors from the model options
  sigma.residual <- model.options$sigma.residual
  sigma.beta <- model.options$sigma.beta
  lambda.numKnots <- model.options$lambda.numKnots
  eta.null <- model.options$eta.null
  
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
  eta.wls <- eta.null + zAlpha + ifelse(!is.null(covariates), cBeta, 0)
  Theta.LSXprev <- wls.binary(y, X.previous, eta.wls, model.options)
  Theta.LSXstar <- wls.binary(y, X.star, eta.wls, model.options)
  Theta.LSresid <- Theta.previous - Theta.LSXprev
  
  # update the coefficients
  residual.deletedKnot <- Theta.LSresid[index]
  Theta.star <- Theta.LSXstar + Theta.LSresid[-index]
  
  
  # calculate the previous eta and associated probability
  eta.previous = as.vector(X.previous %*% Theta.previous + ifelse(!is.null(covariates), cBeta, 0) + zAlpha)
  prob.previous = inv.logit(eta.previous)
  # adjust probabilities within tolerance levels
  prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
  prob.previous[prob.previous > model.options$prob.max] <-  model.options$prob.max
  loglikelihood.previous <- sum(log((1 - prob.previous[outcomes==0]))) + sum(log(prob.previous[outcomes==1]))    
  
  # calculate the proposal eta and associated probability
  eta.star = as.vector(X.star %*% Theta.star + ifelse(!is.null(covariates), cBeta, 0) + zAlpha)
  prob.star = inv.logit(eta.star)
  # adjust probabilities within tolerance levels
  prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
  prob.star[prob.star > model.options$prob.max] <-  model.options$prob.max
  loglikelihood.star <- sum(log((1 - prob.star[outcomes==0]))) + sum(log(prob.star[outcomes==1]))    
  
  # Calculate birth and death probabilities                                                        
  probBirth <- ifelse(length(knots.star) == model.options$knots.min, 
                      1, model.options$knots.prob.birth)
  probDeath <- ifelse(length(knots.previous) == model.options$knots.max, 
                      1, 1 - model.options$knots.prob.birth)
  
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
#' @importFrom MASS ginv
#' @return 
#' @examples
#' 
removeKnot.gaussian <- function(model.options, knots.previous, 
                                outcomes, times.dropout, times.observation, covariates,
                                X.previous, Theta.previous, 
                                Z, alpha, betaCovariates, sigma.error) {
  
  # get the relevant priors from the model options
  sigma.residual <- model.options$sigma.residual
  sigma.beta <- model.options$sigma.beta
  lambda.numKnots <- model.options$lambda.numKnots

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
  probBirth <- ifelse(length(knots.star) == model.options$knots.min, 
                      1, model.options$knots.prob.birth)
  probDeath <- ifelse(length(knots.previous) == model.options$knots.max, 
                      1, 1 - model.options$knots.prob.birth)
  
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
moveKnot.binary <- function(model.options, knots.previous, 
                            outcomes, times.dropout, times.observation, covariates,
                            X.previous, Theta.previous, 
                            Z, alpha, betaCovariates) {
  
  # get the relevant priors from the model options
  # priors
  sigma.residual <- model.options$sigma.residual
  sigma.beta <- model.options$sigma.beta
  lambda.numKnots <- model.options$lambda.numKnots
  eta.null <- model.options$eta.null
  # candidate positions
  potentialLocations = model.options$knots.positions.candidate[! model.options$knots.positions.candidate %in% knots.previous]
  # step size for moving a knot
  knots.stepSize = model.options$knots.stepSize
  
  #Pick a knot to move 
  knotToMove <- sample(knots.previous, 1) 
  # get index of knot to be moved
  index <- which(knots.previous == knotToMove) 
  # get the knots that are staying in the same place
  knotsToKeep <- knots.previous[-index] 
  
  # find a new location from the potential knot locations
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
    eta.previous = as.vector(X.previous %*% Theta.previous + ifelse(!is.null(covariates), cBeta, 0) + zAlpha)
    prob.previous = inv.logit(eta.previous)
    # adjust probabilities within tolerance levels
    prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
    prob.previous[prob.previous > model.options$prob.max] <-  model.options$prob.max
    loglikelihood.previous <- sum(log((1 - prob.previous[outcomes==0]))) + sum(log(prob.previous[outcomes==1]))    
    
    # calculate the proposal eta and associated probability
    eta.star = as.vector(X.star %*% Theta.previous + ifelse(!is.null(covariates), cBeta, 0) + zAlpha)
    prob.star = inv.logit(eta.star)
    # adjust probabilities within tolerance levels
    prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
    prob.star[prob.star > model.options$prob.max] <-  model.options$prob.max
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
moveKnot.gaussian <- function(model.options, knots.previous, 
                              outcomes, times.dropout, times.observation, covariates,
                              X.previous, Theta.previous, 
                              Z, alpha, betaCovariates, sigma.error) {
  
  # get the relevant priors from the model options
  # priors
  sigma.residual <- model.options$sigma.residual
  sigma.beta <- model.options$sigma.beta
  lambda.numKnots <- model.options$lambda.numKnots
  # candidate positions
  potentialLocations = model.options$knots.positions.candidate[! model.options$knots.positions.candidate %in% knots.previous]
  # step size for moving a knot
  knots.stepSize = model.options$knots.stepSize
  
  #Pick a knot to move 
  knotToMove <- sample(knots.previous, 1) 
  # get index of knot to be moved
  index <- which(knots.previous == knotToMove) 
  # get the knots that are staying in the same place
  knotsToKeep <- knots.previous[-index] 
  
  # find a new location from the potential knot locations
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
updateFixedEffects.binary <- function(model.options, knots.previous, 
                                      outcomes, times.dropout, times.observation, covariates,
                                      X.previous, Theta.previous, 
                                      Z, alpha, betaCovariates) {
  
  # get the relevant priors from the model options
  # priors
  sigma.beta <- model.options$sigma.beta
  lambda.numKnots <- model.options$lambda.numKnots
  eta.null <- model.options$eta.null
  
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
  prob.previous[prob.previous > model.options$prob.max] <-  model.options$prob.max
  loglikelihood.previous <- sum(log((1 - prob.previous[y==0]))) + sum(log(prob.previous[y==1]))    
  
  # create the covariance matrix for the intercept and theta coefficients for the splines - R0
  covarIntTheta <- diag(rep(sigma.beta, (length(knots.previous)+1))) 
  covarIntThetaInverse <- diag(rep(1/sigma.beta, (length(knots.previous)+1))) 
  # build y-tilde
  yt = XTheta.previous + (y - prob.previous) * (1 / (prob.previous * (1 - prob.previous)))
  # build the weight matrix
  weight <- Diagonal(x = prob.previous * (1-prob.previous))
  # build the covariance of the beta coefficients and make sure it is positive definite
  covar <- as.matrix(nearPD(solve(covarIntThetaInverse + crossprod(X.previous, as.matrix(weight %*% X.previous))))$mat)
  # build the mean
  mu <- covar  %*% (crossprod(X.previous, as.matrix(weight %*% yt)))
  # draw the proposed coefficients for the fixed effects
  Theta.star <- rmvnorm(1, mu, covar)
  
  # get proposal probabilities
  XTheta.star = X.previous %*% t(Theta.star)
  eta.star <- XTheta.star + ifelse(!is.null(covariates), cBeta, 0) + zAlpha
  prob.star = inv.logit(eta.previous)
  # adjust probabilities within tolerance levels
  prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
  prob.star[prob.star > model.options$prob.max] <-  model.options$prob.max
  loglikelihood.star <- sum(log((1 - prob.star[y == 0]))) + sum(log(prob.star[y == 1]))    
  
  y.star <- XTheta.star + (y - prob.star) * (1 / (prob.star * (1 - prob.star)))
  weight.star <- Diagonal(x=prob.star * (1 - prob.star))
  covar.star <- as.matrix(nearPD(solve(covarIntThetaInverse + crossprod(X.previous, as.matrix(weight.star %*% X.previous))))$mat)
  mu.star <- covar.star %*% (crossprod(X.previous, as.matrix(weight.star %*% y.star)))
  
  # Metropolis hastings step
  rho <- (loglikelihood.star - loglikelihood.previous + 
            log(dmvnorm(Theta.previous, as.vector(mu.star), covar.star)) - 
            log(dmvnorm(Theta.star, as.vector(mu),covar)) +
            (crossprod(Theta.previous) - tcrossprod(Theta.star))/ (2 * sigma.beta))
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
#' @importFrom MASS ginv
#' 
#' @examples
#' add(1, 1)
#' add(10, 1)
updateFixedEffects.gaussian <- function(model.options, knots.previous, 
                                        outcomes, times.dropout, times.observation, covariates,
                                        X.previous, Theta.previous, 
                                        Z, alpha, betaCovariates, sigma.error) {
  
  # get the relevant priors from the model options
  # priors
  sigma.beta <- model.options$sigma.beta
  lambda.numKnots <- model.options$lambda.numKnots
  
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
  adjust = ifelse(!is.null(model.options$accept.rate.adjust), 
                  model.options$accept.rate.adjust, 1)
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
updateFixedEffectsCovariates.binary <- function(model.options, outcomes, covariates, 
                                                X, Theta, Z, alpha, betaCovariates.previous) {
  
  # grab the relevant model options
  sigma.beta <- model.options$sigma.beta
  
  # make sure there are covariates to update 
  if (is.null(covariates)) {
    return (list(betaCovariates=NULL, accepted=FALSE))
  }
  
  # build components of eta
  y <- as.matrix(outcomes)
  C = as.matrix(covariates)
  cBeta = C %*% betaCovariates.previous
  # build X * beta for each group and combine into a single vector
  XTheta = vector()
  for(i in 1:length(X)) {
    XTheta.group = X[[i]] %*% Theta[[i]]
    XTheta <- c(XTheta, XTheta.group) 
  }
  
  zAlpha = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
  # calculate the previous eta and associated probability
  eta.previous = as.vector(XTheta + cBeta + zAlpha)
  prob.previous = inv.logit(eta.previous)
  # adjust probabilities within tolerance levels
  prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
  prob.previous[prob.previous > model.options$prob.max] <-  model.options$prob.max
  loglikelihood.previous <- sum(log((1 - prob.previous[y==0]))) + sum(log(prob.previous[y==1]))    
  
  # create the covariance matrix for the intercept and theta coefficients for the splines - R0
  covarBeta <- diag(ncol(covariates)) * sigma.beta
  covarBetaInverse <- diag(ncol(covariates)) * 1/sigma.beta
  
  # build y-tilde
  yt = cBeta + (y - prob.previous) * (1 / (prob.previous * (1 - prob.previous)))
  # build the weight matrix
  weight <- Diagonal(x = prob.previous * (1-prob.previous))
  # build the covariance of the beta coefficients and make sure it is positive definite
  covar <- as.matrix(nearPD(solve(covarBetaInverse + crossprod(C, as.matrix(weight %*% C))))$mat)
  # build the mean
  mu <- covar  %*% (crossprod(C, as.matrix(weight %*% yt)))
  # draw the proposed coefficients for the fixed effects
  betaCovariates.star <- rmvnorm(1, mu, covar)
  
  # get proposal probabilities
  cBeta.star <- C %*% betaCovariates.star
  eta.star <- XTheta + cBeta.star + zAlpha
  prob.star = inv.logit(eta.star)
  # adjust probabilities within tolerance levels
  prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
  prob.star[prob.star > model.options$prob.max] <-  model.options$prob.max
  loglikelihood.star <- sum(log((1 - prob.star[y == 0]))) + sum(log(prob.star[y == 1]))    
  
  y.star <- cBeta.star + (y - prob.star) * (1 / (prob.star * (1 - prob.star)))
  weight.star <- Diagonal(x=prob.star * (1 - prob.star))
  covar.star <- as.matrix(nearPD(solve(covarBetaInverse + crossprod(C, as.matrix(weight.star %*% C))))$mat)
  mu.star <- covar.star %*% (crossprod(C, as.matrix(weight.star%*%y.star)))

  # Metropolis hastings step
  rho <- (loglikelihood.star - loglikelihood.previous + 
            log(dmvnorm(as.vector(betaCovariates.previous), as.vector(mu.star), covar.star)) - 
            log(dmvnorm(as.vector(betaCovariates.star), as.vector(mu), covar)) +
            (crossprod(as.vector(betaCovariates.previous)) - tcrossprod(betaCovariates.star))/ (2 * sigma.beta))
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
updateFixedEffectsCovariates.gaussian <- function(model.options, outcomes, covariates, 
                                                  X, Theta, Z, alpha, betaCovariates.previous,  
                                                  sigma.error) {
  
  # grab the relevant model options
  sigma.beta <- model.options$sigma.beta
  
  # build X * beta for each group and combine into a single vector
  XTheta = vector()
  for(i in 1:length(X)) {
    XTheta.group = X[[i]] %*% Theta[[i]]
    XTheta <- c(XTheta, XTheta.group) 
  }
  
  # calculate the residuals
  residuals <- outcomes - XTheta - Z$intercept * alpha$intercept - Z$slope * alpha$slope
  residuals = residuals[,outcomes.var]
  # get the proposed mean/variance of the fixed effects associated with covariates
  covariates.matrix = as.matrix(covariates)
  covarBeta <- diag(ncol(covariates)) * sigma.beta
  covarBetaInverse <- diag(ncol(covariates)) * 1/sigma.beta
  
  proposedCovariance <- ginv(covarBetaInverse + 
                               (1/sigma.error) * t(covariates.matrix) %*% covariates.matrix)
  proposedMean <- proposedCovariance %*% 
    ((1/sigma.error) * t(covariates.matrix) %*% as.matrix(residuals))
  
  # Scale Cov to adjust acceptance rate
  proposedCovariance <- proposedCovariance * 
    ifelse(!is.null(model.options$accept.rate.adjust), 
           model.options$accept.rate.adjust, 1)
  
  # ensure the covariance is positive definite
  proposedCovariance <- as.matrix(nearPD(proposedCovariance)$mat) 
  
  # draw a proposed set of coefficients for the covariates
  betaCovariates.star <- rmvnorm(1, proposedMean, proposedCovariance)
  
  # calculate the residuals
  resid.star <- residuals - covariates.matrix %*% t(betaCovariates.star)
  resid.previous <- residuals - covariates.matrix %*% (as.matrix(betaCovariates.previous))
  
  # calculate the acceptance probability
  rho<-sum(log(dnorm(as.vector(resid.star), 0, sqrt(sigma.error)))) + 
    log(dmvnorm(betaCovariates.star, rep(0, ncol(covariates)), covarBetaInverse)) + 
    log(dmvnorm(betaCovariates.previous, proposedMean, proposedCovariance)) - 
    sum(log(dnorm(resid.star, 0, sqrt(sigma.error)))) - 
    log(dmvnorm(betaCovariates.previous, rep(0, ncol(covariates)), covarBetaInverse)) - 
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
updateRandomEffects.binary <- function(numSubjects, numObservations, firstObsPerSubject,
                                       subjectsPerGroup, ids, outcomes, times.observation, 
                                       covariates, X, Theta, 
                                       Z, alpha, betaCovariates,  
                                       sigma.randomIntercept, sigma.randomSlope,
                                       sigma.randomInterceptSlope) {

  # Update B1 and B2 using MH step
  # get previous log likelihood
  y <- as.matrix(outcomes)
  C = as.matrix(covariates)
  if (!is.null(covariates)) {
    cBeta = as.vector(C %*% betaCovariates)
  }
  
  # build X * beta for each group and combine into a single vector
  XTheta = vector()
  for(i in 1:length(X)) {
    XTheta.group = X[[i]] %*% Theta[[i]]
    XTheta <- c(XTheta, XTheta.group) 
  }

  zAlpha.previous = Z[,1] * alpha[,1] + Z[,2] * alpha[,2]
  # calculate the previous eta and associated probability
  eta.previous = as.vector(XTheta + ifelse(!is.null(covariates), cBeta, 0) + zAlpha.previous)
  prob.previous = inv.logit(eta.previous)
  # adjust probabilities within tolerance levels
  prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
  prob.previous[prob.previous > model.options$prob.max] <-  model.options$prob.max
  loglikelihood.previous <- rep(0, length(y))    
  loglikelihood.previous[y==0] <- log(1 - prob.previous[y==0])  
  loglikelihood.previous[y==1] <- log(prob.previous[y==1])
  
  ### Get the proposal log likelihood by drawing a candidate from a symmetric random walk proposal
  # get the random intercepts, one per subject
  alpha.onePerSubject = alpha[firstObsPerSubject,]
  noise <- rmvt(nrow(alpha.onePerSubject),sigma=diag(2),3) # TODO: why df=3??? or b+rmvnorm(n,sigma=Sigma)
  # get the proposed random effects
  alpha.star.intercept <- alpha.onePerSubject[,1] + noise[,1]
  alpha.star.slope <- alpha.onePerSubject[,2] + noise[,2]
  # calculate the proposal log likelihood
  zAlpha.star = Z[,1] * rep(alpha.star.intercept, numObservations) + Z[,2] * rep(alpha.star.slope, numObservations)
  # calculate the previous eta and associated probability
  eta.star = as.vector(XTheta + ifelse(!is.null(covariates), cBeta, 0) + zAlpha.star)
  prob.star = inv.logit(eta.star)
  # adjust probabilities within tolerance levels
  prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
  prob.star[prob.star > model.options$prob.max] <-  model.options$prob.max
  loglikelihood.star <- rep(0, length(y))    
  loglikelihood.star[y==0] <- log(1 - prob.star[y==0])  
  loglikelihood.star[y==1] <- log(prob.star[y==1])
  
  sigma.alpha = matrix(c(sigma.randomIntercept, 
                         sigma.randomInterceptSlope, 
                         sigma.randomInterceptSlope, 
                         sigma.randomSlope), nrow=2, byrow=TRUE)
  alpha.star.onePerSubject = cbind(alpha.star.intercept, alpha.star.intercept)
  
  ratio <- (
    tapply(loglikelihood.star, ids, sum) +
    dmvnorm(alpha.star.onePerSubject, c(0,0), sigma.alpha, log=TRUE) - 
      (tapply(loglikelihood.previous, ids, sum) + 
         dmvnorm(alpha.onePerSubject, c(0,0), sigma.alpha, log=TRUE))
  )
  
  update <- 1*(log(runif(numSubjects)) < ratio)
  randomIntercepts <- alpha.onePerSubject[,1]
  randomIntercepts[update == 1] <- alpha.star.intercept[update==1]
  randomSlopes <- alpha.onePerSubject[,2]
  randomSlopes[update == 1] <- alpha.star.slope[update==1]
  
  # build the new random effects matrix
  alpha.total = data.frame(intercept=randomIntercepts, slope=randomSlopes)
  alpha = alpha.total[rep(seq_len(nrow(alpha.total)), numObservations),]
  # expand back out to complete alpha matrix
  return (alpha)
  
}



#' Update the random effects
#' 
#' @param 
#' @return 
#' @examples
#'
updateRandomEffects.gaussian <- function(numSubjects, numObservations, firstObsPerSubject,
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
updateCovarianceParameters.binary <- function(model.options, numSubjects, 
                                              firstObsPerSubject, alpha) {
  
  
  # sample from an inverse wishart to update the covariance of the random effects
  perSubjectAlpha = alpha[firstObsPerSubject,]
  sigma.alpha <- riwish(model.options$sigma.randomEffects.df + numSubjects, 
                        model.options$sigma.randomEffects.scale + crossprod(as.matrix(perSubjectAlpha)))
  
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
updateCovarianceParameters.gaussian <- function(model.options, totalObservations, numSubjects, 
                                                firstObsPerSubject,
                                                outcomes, covariates, X, Theta, 
                                                Z, alpha, betaCovariates,
                                                sigma.error) {
  
  # sample from an inverse wishart to update the covariance of the random effects
  perSubjectAlpha = alpha[firstObsPerSubject,]
  sigma.alpha <- riwish(model.options$sigma.randomEffects.df + numSubjects, 
                        model.options$sigma.randomEffects.scale + crossprod(as.matrix(perSubjectAlpha)))
  
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
  shape <- model.options$sigma.error.shape.tau + (totalObservations / 2) 
  rate <- model.options$sigma.error.rate.tau + (crossprod(residuals) / 2)
  sigma.error <- 1 / rgamma(1, shape, rate)
  
  return (list(sigma.error = sigma.error,
               sigma.randomIntercept = sigma.alpha[1,1],
               sigma.randomSlope = sigma.alpha[2,2],
               sigma.randomInterceptSlope = sigma.alpha[1,2]))
  
}

#'
#' Calculate the marginal slope at the given iteration
#'
#' @param knotsByGroup list of knots by group
#' @param ThetaByGroup list of spline coefficients by group
#' @param subjectsPerGroup list of number of subjects by group
#' @param times.dropout dropout times
#' 
#' @importFrom splines ns
#' @importFrom gtools rdirichlet
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

#'
#' Calculate the estimated slope at the given dropout times
#' 
#' @param dropoutEstimationTimes list of times at which to estimate the slope
#' @param knotsByGroup list of knots by group
#' @param ThetaByGroup list of spline coefficients by group
#' @param subjectsPerGroup list of number of subjects by group
#' @param times.observation observation times
#' @param times.dropout dropout times
#' 
#' @importFrom splines ns
#'
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
#' calculate the acceptance probability for a given update step in the Bayesian
#' spline model run.
#' 
#' @param fit the bayes.splines.fit object
#' @param action the update action.  Valid values are "knot.add", "knot.remove",
#' "knot.move", "fixedEffects", "fixedEffectsCovariates"
#'
#' @export prob.acceptance
#' 
prob.acceptance <- function(fit, action) {
  iterations = fit$iterations
  if (action == "knot.add") {
    total_accepts = sum(sapply(iterations, function(x) { return(as.numeric(x$accepted$knot.add)) }))
  } else if (action == "knot.remove") {
    total_accepts = sum(sapply(iterations, function(x) { return(as.numeric(x$accepted$knot.remove)) }))
  } else if (action == "knot.move") {
    total_accepts = sum(sapply(iterations, function(x) { return(as.numeric(x$accepted$knot.move)) }))
    
  } else if (action == "fixedEffects") {
    total_accepts = sum(sapply(iterations, function(x) { return(as.numeric(x$accepted$fixedEffects)) }))
  } else if (action == "fixedEffectsCovariates") {
    total_accepts = sum(sapply(iterations, function(x) { return(as.numeric(x$accepted$fixedEffectsCovariates)) }))
  }
  
  
  return(total_accepts/length(iterations))
}


#' Summarize a Bayesian spline model run
#'
#' @param fit bayes.splines.fit object from an informativeDropout run.
#'
#' @export summary.bayes.splines.fit
#' 
summary.bayes.splines.fit <- function(fit) {
  if (length(fit$iterations) <= 0) {
    stop("No results found in Bayesian spline model fit")
  }
  dist = fit$dist
  model.options = fit$model.options
  iterations = fit$iterations
  groups = fit$groups
  covariates.var = fit$covariates.var
  
  result.summary = list()
  
  for (group.index in 1:length(groups)) {
    group = groups[group.index]
    
    ## number of knots
    knots_sample = unlist(lapply(iterations, function(x) { 
      return(length(x$knots[[group.index]]))
    }))
    knots = data.frame(
      mean=mean(knots_sample),
      median=median(knots_sample),
      ci_lower=quantile(knots_sample, probs=0.025),
      ci_upper=quantile(knots_sample, probs=0.975)
    )
    row.names(knots) = c("knots")
    
    ## marginal slope
    marginal_slope_sample = unlist(lapply(iterations, function(x) { 
      slope = (x$slope.marginal[[group.index]])
      return (ifelse(is.null(slope), 0, slope))
    }))
    marginal_slope = data.frame(
      mean=mean(marginal_slope_sample),
      median=median(marginal_slope_sample),
      ci_lower=quantile(marginal_slope_sample, probs=0.025),
      ci_upper=quantile(marginal_slope_sample, probs=0.975)
    )
    row.names(marginal_slope) = c("marginal slope")
    
    ## dropout time specific slopes
    slopes_by_dropout_time = data.frame(
      time=model.options$dropout.estimationTimes,
      mean=numeric(length(model.options$dropout.estimationTimes)),
      median=numeric(length(model.options$dropout.estimationTimes)),
      ci_lower=numeric(length(model.options$dropout.estimationTimes)),
      ci_upper=numeric(length(model.options$dropout.estimationTimes))
    )
    for (time.index in 1:(length(model.options$dropout.estimationTimes))) {
      slope_sample <- unlist(lapply(iterations, function(x) { 
        if (is.null(x$slope.dropoutSpecific)) {
          return(0)
        } 
        return(x$slope.dropoutSpecific[[group.index]][time.index])
      }))
      slopes_by_dropout_time$mean[time.index] = mean(slope_sample)
      slopes_by_dropout_time$median[time.index] = median(slope_sample)
      slopes_by_dropout_time$ci_lower[time.index] = quantile(slope_sample, na.rm=TRUE, probs=0.025)
      slopes_by_dropout_time$ci_upper[time.index] = quantile(slope_sample, na.rm=TRUE, probs=0.975)
    }  
    row.names(slopes_by_dropout_time) = NULL
    
    # build the summary list
    result.summary[[paste("group", group, sep="_")]] = list(
      knots = knots,
      marginal_slope = marginal_slope,
      slopes_by_dropout_time = slopes_by_dropout_time
    )
  }
  
  ## covariance parameters
  param_list = c("sigma.residual","sigma.randomIntercept",
                 "sigma.randomSlope", "sigma.randomInterceptSlope")
  if (dist == "gaussian") {
    param_list = c(param_list, "sigma.error")
  }
  covariance_parameters = data.frame(
    mean=numeric(length(param_list)),
    median=numeric(length(param_list)),
    ci_lower=numeric(length(param_list)),
    ci_upper=numeric(length(param_list))
  )
  row = 1
  for (param in param_list) {
    param_sample <- unlist(lapply(iterations, function(x) { 
      return(x[[param]])
    }))
    covariance_parameters$mean[row] = mean(param_sample)
    covariance_parameters$median[row] = median(param_sample)
    covariance_parameters$ci_lower[row] = quantile(param_sample, probs=0.025)
    covariance_parameters$ci_upper[row] = quantile(param_sample, probs=0.975)
    row = row + 1
  }
  row.names(covariance_parameters) = param_list
  result.summary$covariance_parameters = covariance_parameters
  
  ## summarize covariate effects
  if (!is.null(covariates.var)) {
    num_covariates = length(covariates.var)
    covariate_effects = data.frame(mean=numeric(num_covariates), median=numeric(num_covariates),
                                   ci_lower=numeric(num_covariates), ci_upper=numeric(num_covariates))
    for(i in 1:length(covariates.var)) {
      beta_covariate_sample = unlist(lapply(iterations, function(x) { 
        return(x$betas.covariates[i])
      }))
      covariate_effects$mean[i] = mean(beta_covariate_sample)
      covariate_effects$median[i] = median(beta_covariate_sample)
      covariate_effects$ci_lower[i] = quantile(beta_covariate_sample, probs=0.025)
      covariate_effects$ci_upper[i] = quantile(beta_covariate_sample, probs=0.975)
    }
    row.names(covariate_effects) <- covariates.var
    result.summary$covariate_effects = covariate_effects
  }
  
  if (dist == "gaussian") {
    ## residual covariance parameters
    sigma_error_sample = unlist(lapply(iterations, function(x) { 
      return(x$sigma.error)
    }))
    sigma_error_result = data.frame(
      mean = mean(sigma_error_sample),
      median = median(sigma_error_sample),
      ci_lower = quantile(sigma_error_sample, probs=0.025),
      ci_upper = quantile(sigma_error_sample, probs=0.975)
    )
    row.names(sigma_error_result) = "sigma.error"
    
    result.summary$sigma.error = sigma_error_result
  }
  
  
  # calculate acceptance probability
  acceptance_probabilities = data.frame(probability=c(
    prob.acceptance(fit, "knot.add"),
    prob.acceptance(fit, "knot.remove"),
    prob.acceptance(fit, "knot.move"),
    prob.acceptance(fit, "fixedEffects"),
    prob.acceptance(fit, "fixedEffectsCovariates")
  ))
  row.names(acceptance_probabilities) = c(
    "knot.add", "knot.remove", "knot.move", "fixedEffects", "fixedEffectsCovariates"
    )
  result.summary$acceptance_probabilities = acceptance_probabilities
  
  return(result.summary)
  
}




#'
#' Fit a varying coefficient model for longitudinal studies with
#' informative dropout.  Uses a Bayesian approach with a spline fit
#' to model the relationship between dropout times and slope
#'
#' @param data the data set
#' @param ids.var column of the data set containing the participant id
#' @param outcomes.var column of the data set containing the outcome variable
#' @param groups.var column of the data set indicating treatment group
#' @param covariates.var list of columns of the data set containing covariates
#' @param times.dropout.var column of the data set containing dropout times
#' @param times.observation.var column of the data set containing observation times
#' @param dist distribution of the outcome (either "gaussian" or "binary")
#' @param model.options the bayes.spines.model.options object (@@seealso bayes.spines.model.options)
#'  
#' @importFrom splines ns
#' @importFrom abind abind
#' @importFrom gtools logit
#'
#' @export informativeDropout.bayes.splines
#' 
informativeDropout.bayes.splines <- function(data, ids.var, outcomes.var, groups.var, 
                                             covariates.var, times.dropout.var, times.observation.var, 
                                             dist, model.options) {

  # create some reasonable defaults for the start positions and candidate positions
  # if not specified
  if (is.null(model.options$knots.positions.start) || is.null(model.options$knots.positions.candidate)) {
    dropout.min = min(data[,times.dropout.var])
    dropout.max = max(data[,times.dropout.var])
    if (is.null(model.options$knots.positions.candidate)) {
      model.options$knots.positions.candidate = seq(dropout.min, dropout.max, 1)
    }
    if (is.null(model.options$knots.positions.start)) {
      model.options$knots.positions.start = sample(model.options$knots.positions.candidate, 
                                                   model.options$knots.min)
    }
  }
  
  if (dist == "gaussian") {
    if (is.null(model.options$sigma.error) || is.na(model.options$sigma.error) || model.options$sigma.error <= 0) {
      stop("Prior options error :: for Gaussian outcomes sigma.error must be greater than 0")
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
    knots.boundary = range(model.options$knots.positions.start)
    knots.interior = model.options$knots.positions.start[-c(1,length(model.options$knots.positions.start))] 
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
      eta.null = logit(sum(outcomes[,outcomes.var])/length(outcomes[,outcomes.var]))
    }
  }
  
  
  # initialize the first model iteration
  # TODO: mcmc options specify starting values
  iterations.saved = ceiling((model.options$iterations - model.options$burnin) / model.options$thin)
  modelIterationList <- vector(mode = "list", length = iterations.saved)
  
  model.previous = 
    bayes.splines.iteration(
      knots=lapply(1:length(groupList), function(i) { return (model.options$knots.positions.start); }), 
      Theta=Theta.init, betas.covariates=betaCovariate.init,
      sigma.error = model.options$sigma.error,
      sigma.residual = model.options$sigma.residual,
      sigma.randomIntercept = model.options$sigma.randomIntercept,
      sigma.randomSlope = model.options$sigma.randomSlope,
      sigma.randomInterceptSlope = model.options$sigma.randomInterceptSlope,
      sigma.error.shape = model.options$sigma.error.shape.tau,
      sigma.error.rate = model.options$sigma.error.rate.tau
    )
  
  #
  # Run the reversible jump MCMC
  #
  iterations.saved.next = 1
  for (i in 1:model.options$iterations) {
    
    if (i %% model.options$print == 0) {
      print(paste("Bayesian spline model iteration = ", i, sep=""))
    }
    
    # make a copy of the previous iteration which will be modified as we move through the iteration
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
      if ((u < model.options$knots.prob.birth && 
           length(model.current$knots[[group.index]]) < model.options$knots.max) || 
          (length(model.current$knots[[group.index]]) <= model.options$knots.min)) {
        if (dist == "gaussian") {
          # add a knot
          result = addKnot.gaussian(model.options, model.current$knots[[group.index]], 
                                    group.outcomes, group.times.dropout, group.times.observation, 
                                    group.covariates, X.previous=X[[group.index]], 
                                    model.current$Theta[[group.index]], group.Z, group.alpha, 
                                    model.current$betas.covariates,  
                                    sigma.error=model.current$sigma.error)
        } else {
          # binary case
          result = addKnot.binary(model.options, model.current$knots[[group.index]], group.outcomes, 
                                  group.times.dropout, group.times.observation, 
                                  group.covariates, X[[group.index]], model.current$Theta[[group.index]],
                                  group.Z, group.alpha, model.current$betas.covariates)
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
          result = removeKnot.gaussian(model.options, model.current$knots[[group.index]], 
                                       group.outcomes, group.times.dropout, group.times.observation, 
                                       group.covariates, X[[group.index]], model.current$Theta[[group.index]],
                                       group.Z, group.alpha, model.current$betas.covariates, 
                                       model.current$sigma.error)
        } else {
          # remove a knot
          result = removeKnot.binary(model.options, model.current$knots[[group.index]], 
                                     group.outcomes, group.times.dropout, 
                                     group.times.observation, group.covariates,
                                     X[[group.index]], model.current$Theta[[group.index]], 
                                     group.Z, group.alpha, model.current$betas.covariates)
        }
        # update the model iteration
        X[[group.index]] = result$X
        model.current$knots[[group.index]] = result$knots
        model.current$Theta[[group.index]] = result$Theta
        model.current$proposed$knot.remove = TRUE
        model.current$accepted$knot.remove = result$accepted
      }

      # Move knots
      if (dist == "gaussian") {
        result = moveKnot.gaussian(model.options, model.current$knots[[group.index]], 
                                   group.outcomes, group.times.dropout, group.times.observation, 
                                   group.covariates,
                                   X[[group.index]], model.current$Theta[[group.index]], 
                                   group.Z, group.alpha, model.current$betas.covariates, 
                                   model.current$sigma.error)
      } else {
        result = moveKnot.binary(model.options, model.current$knots[[group.index]], 
                                 group.outcomes, group.times.dropout, group.times.observation, 
                                 group.covariates,
                                 X[[group.index]], model.current$Theta[[group.index]], 
                                 group.Z, group.alpha, model.current$betas.covariates)
      }
      X[[group.index]] = result$X
      model.current$knots[[group.index]] = result$knots
      model.current$proposed$knot.move = result$proposed
      model.current$accepted$knot.move = result$accepted
      
      # update fixed effects (includes coefficients for covariates and time varying slopes)
      if (dist == "gaussian") {
        result = updateFixedEffects.gaussian(model.options, model.current$knots[[group.index]], 
                                             group.outcomes, group.times.dropout, group.times.observation, 
                                             group.covariates, X[[group.index]], model.current$Theta[[group.index]], 
                                             group.Z, group.alpha, model.current$betas.covariates, 
                                             model.current$sigma.error)
      } else {
        result = updateFixedEffects.binary(model.options, model.current$knots[[group.index]], 
                                           group.outcomes, group.times.dropout, group.times.observation, 
                                           group.covariates, X[[group.index]], model.current$Theta[[group.index]], 
                                           group.Z, group.alpha, model.current$betas.covariates)
      }
      model.current$Theta[[group.index]] = result$Theta
      model.current$proposed$fixedEffects = TRUE
      model.current$accepted$fixedEffects = result$accepted
      
    }  
    
    # update fixed effects associated with covariates
    if (!is.null(covariates)) {
      if (dist == "gaussian") {
        result = updateFixedEffectsCovariates.gaussian(model.options, outcomes, covariates, 
                                                       X, model.current$Theta, 
                                                       Z, alpha, model.current$betas.covariates,
                                                       model.current$sigma.error)
      } else {
        # binary outcome
        result = updateFixedEffectsCovariates.binary(model.options, outcomes, covariates, 
                                                     X, model.current$Theta, 
                                                     Z, alpha, model.current$betas.covariates)
      }
      model.current$betas.covariates = result$betaCovariates
      model.current$proposed$fixedEffectsCovariates = TRUE
      model.current$accepted$fixedEffectsCovariates = result$accepted
    }

    # update random effects
    if (dist == "gaussian") {
      alpha = 
        updateRandomEffects.gaussian(numSubjects, numObservations, firstObsPerSubject,
                                     subjectsPerGroup, data[,ids.var], outcomes, 
                                     data[,times.observation.var], covariates, 
                                     X, model.current$Theta, Z, alpha, 
                                     model.current$betas.covariates,
                                     model.current$sigma.randomIntercept, 
                                     model.current$sigma.randomSlope,
                                     model.current$sigma.randomInterceptSlope,
                                     model.current$sigma.error)
    } else {
      alpha = 
        updateRandomEffects.binary(numSubjects, numObservations, firstObsPerSubject,
                                   subjectsPerGroup, data[,ids.var], outcomes, 
                                   data[,times.observation.var], covariates, 
                                   X, model.current$Theta, Z, alpha, 
                                   model.current$betas.covariates,
                                   model.current$sigma.randomIntercept, 
                                   model.current$sigma.randomSlope,
                                   model.current$sigma.randomInterceptSlope)
    }
    # update variance components
    if (dist == "gaussian") {
      result = updateCovarianceParameters.gaussian(model.options, nrow(data), 
                                                   numSubjects, firstObsPerSubject,
                                                   outcomes, covariates, X, model.current$Theta,
                                                   Z, alpha, model.current$betas.covariates,
                                                   sigma.error)
    } else {
      result = updateCovarianceParameters.binary(model.options, numSubjects, firstObsPerSubject, alpha)
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
      dropoutEstimationTimes=model.options$dropout.estimationTimes, 
      knotsByGroup = model.current$knots, 
      ThetaByGroup = model.current$Theta, 
      subjectsPerGroup=subjectsPerGroup,
      times.observation=data[,times.observation.var],
      times.dropout=data[firstObsPerSubject,times.dropout.var]
    )
    model.current$slope.dropoutSpecific = result
    
    # save the current iteration
    model.previous = model.current
    if ((i %% model.options$thin == 0 && i > model.options$burnin &&
         iterations.saved.next <= length(modelIterationList)) || 
        iterations.saved.next == length(modelIterationList)) {
      modelIterationList[[iterations.saved.next]] = model.current
      iterations.saved.next = iterations.saved.next + 1
    }
  }
  
  # return the estimates, with distributions, and the model results from each iteration
  return(bayes.splines.fit(model.options, dist, groupList, covariates.var, modelIterationList))
  
}


