#
# Supporting functions for the RJMCMC loop
#
#
#

#'
#' Calculate the weighted least squares solution for non-Gaussian outcomes
#' 
#' 
irls.binary <- function(y, x, covariates, eta.wls, B1, B2, nobs, niter) { 
  prev = NULL
  for (i in 1:niter) {
    # get the probability from the linear predictor
    pt <- inv.logit(eta.wls)
    # truncate at min/max values for computational stability
    pt[pt<0.00001]<-0.00001
    pt[pt>0.99999]<-0.99999
    
    # calculate the residuals
    ywls <- etawls-rep(B1,nobs)-t*rep(B2,nobs) + (y-pt)*(1/pt)*(1/(1-pt))
    weight<-Diagonal(x=pt*(1-pt))
    mnt<-solve(as.matrix(nearPD(crossprod(x,weight%*%x))$mat))%*%(crossprod(x,weight%*%ywls))
    etawls<-as.vector(x%*%mnt)+rep(B1,nobs)+t*rep(B2,nobs)
  }
  return(as.vector(mnt))
} 

#' A class containing the current state of the model
#' during a RJMCMC run
#'
rjmcmc.iteration <- function(knots=NULL, Theta=NULL, betaCovariates=NULL,
                             sigma.error = 1, sigma.spline = 1,
                             sigma.randomIntercept = 1,
                             sigma.randomSlope = 1, sigma.randomInterceptSlope = 0,
                             shape.tau = 0.001, rate.tau = 0.001) {
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
    betaCovariates = betaCovariates,
    
    # residual error variance
    sigma.error = sigma.error,
    # variance of spline coefficients
    sigma.spline = sigma.spline,
    # variance / covariance of the random effects
    sigma.randomIntercept = sigma.randomIntercept,
    sigma.randomSlope = sigma.randomSlope,
    sigma.randomInterceptSlope = sigma.randomInterceptSlope,
    
    # hyperparameters describing the gamma distribution of the variance components
    shape.tau = shape.tau,
    rate.tau = rate.tau,
    
    # indicates which changes to the model were proposed at this iteration
    proposed=list(knot.add=FALSE, knot.remove=FALSE, knot.move=FALSE,
                  fixedEffects=FALSE, fixedEffectsCovariates=FALSE),
    # incidates which changes were accepted at this iteration
    accepted=list(knot.add=FALSE, knot.remove=FALSE, knot.move=FALSE,
                  fixedEffects=FALSE, fixedEffectsCovariates=FALSE)
  )
  
  class(mi) <- append(class(mi), "rjmcmc.iteration")
  return(mi)
  
}

#'
#'Convenience function to get the name of a group indicator column
#'
groupIndicatorColumn <- function(group) {
  return (paste("group.", group, sep=""))
}

#' Add a new knot and use a Metropolis-Hastings step
#' to determine if we keep the changes to the model
#' 
#' @param 
#' @param 
#' @return 
#' @examples
#' 
addKnot <- function(dist, knots.previous, knots.options, 
                    outcomes, times.dropout, times.observation, covariates,
                    X.previous, Theta.previous, 
                    Z, alpha, betaCovariates, sigma.residual, 
                    sigma.error, sigma.beta, lambda.numKnots) {
  
  # add a knot by randomly selecting a candidate knot
  candidatesPositions = knots.options$candidatePositions[! knots.options$candidatePositions %in% knots.previous]
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
  
  if (dist == 'binary') {
    eta.wls = eta0 + zAlpha + ifelse(!is.null(covariates), cBeta, 0)
    Theta.LSXprev <- wls(y, X.previous, eta.wls, B1, B2, nobs,1)
    Theta.LSXstar <-wls(y, Xstar, eta.wls, B1, B2, nobs,1)
  } else {
    # get the residuals
    yls <- as.vector(y - zAlpha)
    if (!is.null(covariates)) {
      yls <- (yls - cBeta)
    } 
    # Calculate least squares estimates for coefficients and differences between LS and current coefficients
    Theta.LSXprev <- ginv(crossprod(X.previous))%*%(crossprod(X.previous,yls))
    Theta.LSXstar <- ginv(crossprod(X.star))%*%(crossprod(X.star,yls))
  }
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
  probBirth <- ifelse(length(knots.previous) == knots.options$min, 1, knots.options$birthProbability)
  probDeath <- ifelse(length(knots.previous) == knots.options$max - 1, 1, 1 - knots.options$birthProbability)
  
  # Calculate residuals for likelihood ratio
  if (dist == 'gaussian') {
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
  } else if (dist == 'binary') {
    # calculate the residuals -- TODO
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
  }

  


  
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
removeKnot <- function(dist, knots.previous, knots.options, 
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
  probBirth <- ifelse(length(knots.star) == knots.options$min, 
                      1, knots.options$birthProbability)
  probDeath <- ifelse(length(knots.previous) == knots.options$max, 
                      1, 1 - knots.options$birthProbability)
  
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
moveKnot <- function(dist, knots.previous, knots.stepSize, knots.candidatePositions,
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
updateFixedEffects <- function(dist, knots.previous, knots.options, 
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
updateFixedEffectsCovariates <- function(dist, outcomes, covariates, 
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
updateRandomEffects <- function(dist, numSubjects, numObservations, firstObsPerSubject,
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
updateCovarianceParameters <- function(dist, totalObservations, numSubjects, 
                                       firstObsPerSubject,
                                       outcomes, covariates, X, Theta, 
                                       Z, alpha, betaCovariates,
                                       sigma.error,
                                       prior.options) {
  
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
  
  # sample from an inverse wishart to update the covariance of the random effects
  perSubjectAlpha = alpha[firstObsPerSubject,]
  sigma.alpha <- riwish(prior.options$sigmaError.df + numSubjects, 
                        prior.options$sigmaError.scaleMatrix + crossprod(as.matrix(perSubjectAlpha)))
  
  sigma.randomIntercept = sigma.alpha[1,1]
  sigma.randomSlope = sigma.alpha[2,2]
  sigma.randomInterceptSlope = sigma.alpha[1,2]
  
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
  } else {
    family = ifelse(dist == 'poisson', poisson(link = "log"), binomial(link = "logit"))
    fit.beta <- glm(formula, family=family, data=data.covar)
    return (as.vector(coef(fit.beta))[-1])
  }
}
