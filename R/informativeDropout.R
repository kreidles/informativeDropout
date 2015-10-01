#####################################################################
# 
#  Package informativeDropout implements Bayesian and Frequentist
#  approaches for fitting varying coefficient models in longitudinal
#  studies with informative dropout
#
#  Copyright (C) 2014 University of Colorado Denver.
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#####################################################################

library(matrixcalc)
library(splines)
#'
#' A class describing a bayesian spline fit.  Includes the
#' parameter estimates and iteration summary
#' 
#'

#'
#' A class descirbing a bayesian dirichlet process fit.
#'
#'
#'

#'
#'A class describing a mixed model fit 
#'


#' A class containing the current state of the model
#' during a RJMCMC run
#'
rjmcmc.iteration <- function(knots=NULL, Theta=NULL, betaCovariate=NULL,
                             sigma.error = 1, sigma.randomIntercept = 1,
                             sigma.randomSlope = 1, sigma.randomInterceptSlope = 0,
                             lambda.numKnots = 1, sigma.beta = 1,
                             shape.tau = 0.001, rate.tau = 0.001) {
  mi = list(
    knots = knots,
    # list of per group intercept and spline coefficients
    Theta = Theta,
    # regression coefficients associated with the shared covariates 
    # (includes group offsets)
    betaCovariates = betaCovariates,
    
    # residual error variance
    sigma.error = sigma.error,
    # variance / covariance of the random effects
    sigma.randomIntercept = sigma.randomIntercept,
    sigma.randomSlope = sigma.randomSlope,
    sigma.randomInterceptSlope = sigma.randomInterceptSlope,
    
    # Poisson parameters for the number of knots, D, where D ~ Poisson(lambda)
    lambda.numKnots = lambda.numKnots,
    
    # variance of the regression coefficients (assumes equal variance of each coefficient)
    sigma.beta = sigma.beta,
    
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
addKnot <- function(knots.previous, knots.options, 
                    outcomes, times.dropout, times.observation, covariates,
                    X.previous, Theta.previous, 
                    Z, alpha, betaCovariates,  
                    sigma.error, sigma.beta, lambda.numKnots) {
  
  # add a knot by randomly selecting a candidate knot
  index = sample(1:length(knots.options$candidatePositions), 1)
  newKnot.value = knots.options$candidatePositions[index]
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
  cBeta = covariates %*% betaCovariates
  zAlpha = Z %*% alpha
  
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
  residual <- rnorm(1, 0, sigma.error)
  beta.newKnot <- Theta.LSXstar[newKnot.position+1] + residual
  beta.other <- Theta.LSXstar[-(newKnot.position+1)] + Theta.LSresid
  # insert the new beta value in the correct position
  if (newKnot.position == 1) {
    Theta.star = c(beta.other[1], beta.newKnot, beta.other[2:length(beta.other)])
  } else if (newKnot.position == length(knots)) {
    Theta.star = c(beta.other, beta.newKnot)
  } else {
    Theta.star = c(beta.other[1:newKnot.position], beta.newKnot, beta.other[(newKnot.position+1):length(beta.other)])
  }
  
  # Calculate residuals for likelihood ratio
  if (!is.null(covariates)) {
    LRresid.star <- as.vector(y - X.star %*% Theta.star - cBeta - zAlpha)
    LRresid.prev <- as.vector(y - X.previous %*% Theta.previous - cBeta - zAlpha)
  } else {
    LRresid.star <-as.vector(y - X.star %*% Theta.star - zAlpha)
    LRresid.prev <- as.vector(y - X.previous %*% Theta.previous - zAlpha)
  }
  
  # Calculate birth and death probabilities                                                        
  probBirth <- ifelse(length(knots.previous) == knots.options$min, 1, knot.options$birthProbability)
  probDeath <- ifelse(length(knots.previous) == knots.options$max - 1, 1, 1 - knot.options$birthProbability)
  
  #Calculate Acceptance Probability                                                        
  rho <- (log(lambda.numKnots) - 
            log(length(knots)) + 
            log(probDeath) - log(probBirth) + 
            log(sqrt(sigma.error)) - 
            0.5 * log(sigma.beta) +
            (residual^2 / (2 * sigma.error * sigma.error)) + 
            ((crossprod(Theta) - crossprod(Theta.star)) / (2 * sigma.beta)) +
            ((crossprod(LRresid.prev) - crossprod(LRresid.star)) / (2 * sigma.beta))
  )
  
  if (rho > log(runif(1))) {
    return (list(X=X.star, knots=knots.star, Theta=Theta.star, accepted=true))
  } else {
    return (list(X=X.previous, knots=knots.previous, Theta=Theta.previous, accepted=false))          
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
removeKnot <- function(knots.previous, knots.options, 
                       outcomes, times.dropout, times.observation, covariates,
                       X.previous, Theta.previous, 
                       Z, alpha, betaCovariates,  
                       sigma.error, sigma.beta, lambda.numKnots) {
  
  # randomly remove an existing knot
  index = sample(1:length(knots.previous), 1)
  knots.star <- knots.previous[-index]
  
  knots.boundary = range(knots.star)
  knots.interior = knots.star[-c(1,length(knots.star))] 
  
  if (length(knots.star) > 1 ) {  
    X.star<-cbind(
      rep(1, length(times.dropout)),
      ns(times.dropout, knots=knotstar.i, Boundary.knots=knotstar.b, intercept=T) * times.observation
    )
  } else {
    X.star<-cbind(rep(1, nrow(data)), data[,times])
  }
  
  # Calculate residuals
  y = as.matrix(outcomes)
  cBeta = covariates %*% betaCovariates
  zAlpha = Z %*% alpha
  
  yls <- as.vector(y - zAlpha)
  if (!is.null(covariates)) {
    yls <- (yls - cBeta)
  } 
  
  # Calculate least squares estimates for coefficients and differences between LS and current coefficients
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
  probBirth <- ifelse(length(modelIteration.previous$knots) == knots.options$min, 1, knots.options$birthProbability)
  probDeath <- ifelse(length(modelIteration.previous$knots) == knots.options$max - 1, 1, 1 - knots.options$birthProbability)
  
  #Calculate Acceptance Probability                                                        
  rho <- (log(lambda.numKnots) - 
            log(length(knots.previous)) + 
            log(probDeath) - log(probBirth) + 
            log(sqrt(sigma.error)) - 
            0.5 * log(sigma.beta) +
            (residual^2 / (2 * sigma.error * sigma.error)) + 
            ((crossprod(Theta.previous) - crossprod(Theta.star)) / (2 * sigma.beta)) +
            ((crossprod(LRresid.prev) - crossprod(LRresid.star)) / (2 * sigma.beta))
  )
  
  if (1/rho > log(runif(1))) {
    return (list(X=X.star, knots=knots.star, Theta=Theta.star, accepted=true))
  } else {
    return (list(X=X.previous, knots=knots.previous, Theta=Theta.previous, accepted=false))         
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
moveKnot <- function(knots.previous, knots.stepSize, knots.candidatePositions,
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
  potentialLocations <- knots.candidatePositions[!(knots.candidatePositions %in% knots.previous)]
  # here we only allow movement within some small window
  potentialLocations <- potentialLocations[potentialLocations > (index - knots.stepSize) & 
                                             potentialLocations < (index + knots.stepSize)] 
  
  if (length(potentialLocations > 0)) {
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
    cBeta = covariates %*% betaCovariates
    zAlpha = Z %*% alpha
    
    if (!is.null(covariates)) {
      LRresid.star <- as.vector(y - X.star %*% Theta.star - cBeta - zAlpha)
      LRresid.prev <- as.vector(y - X.previous %*% Theta.previous - cBeta - zAlpha)
    } else {
      LRresid.star <-as.vector(y - X.star %*% Theta.star - zAlpha)
      LRresid.prev <- as.vector(y - X.previous %*% Theta.previous - zAlpha)
    }
    
    # Calculate the acceptance probability
    rho <- (log(length(potentialknots)) - log(length(knots.star)) + 
      (crossprod(LRresid.prev)-crossprod(LRresid.star))/2/sigma.error)
    if (rho > log(runif(1))) {
      return (knots=knots.star, X=X.star, proposed=true, accepted=true)
    } else {
      return (knots=knots.previous, X=X.previous, proposed=true, accepted=false)
    }
  } else {
    # nowhere to move to
    return (knots=knots.previous, X=X.previous, proposed=false, accepted=false)
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
updateFixedEffects <- function(knots.previous, knots.options, 
                               outcomes, times.dropout, times.observation, covariates,
                               X.previous, Theta.previous, 
                               Z, alpha, betaCovariates,  
                               sigma.error, sigma.beta, lambda.numKnots) {
  
  # create the covariance matrix for the intercept and theta coefficients for the splines - R0
  covarIntTheta <- diag(rep(sigma.beta, (length(knots.previous)+1))) 
  covarIntThetaInverse <- diag(rep(1/sigma.beta, (length(knots.previous)+1))) 
  # Calculate residuals 
  y = as.matrix(outcomes)
  cBeta = covariates %*% betaCovariates
  zAlpha = Z %*% alpha

  if (!is.null(covariates)) {
    LRresid <- as.vector(y - X.previous %*% Theta.previous - cBeta - zAlpha)
  } else {
    LRresid <- as.vector(y - X.previous %*% Theta.previous - zAlpha)
  }
  
  # calculate the covariance and mean of the proposal distribution
  proposedCovariance <- ginv(covarIntThetaInverse + (1/sigma.beta) * crossprod(X.previous))
  proposedMean <- proposedCovariance %*% ((1/sigma.beta) * crossprod(X.previous, LRresid))
  
  # Scale Cov to adjust acceptance rate
  proposedCovariance <- mcmc.options$fixedEffectAcceptRateAdjust * proposedCovariance 
  
  # ensure the covariance is positive definite
  proposedCovariance <- as.matrix(nearPD(proposedCovariance)$mat) 
  # draw a proposed set of coefficients
  Theta.star <- rmvnorm(1, proposedMean, proposedCovariance) 
  
  # Calculate residuals for likelihood ratio
  resid.star <- LRresid - crossprod(X.previous, Theta.star)
  resid.prev <- LRresid - crossprod(X.previous, Theta.previous)
  
  # Calculate acceptance probability
  rho <- (log(dmvnorm(Theta.previous, proposedMean, proposedCovariance)) - 
            log(dmvnorm(Theta.star, proposedMean, proposedCovariance)) + 
            (crossprod(Theta.previous)-tcrossprod(Theta.star))/(2 * sigma.beta) + 
            (crossprod(resid.prev)-crossprod(resid.star))/(2 * sigma.beta))
  
  if (rho > log(runif(1))) {
    return (Theta=Theta.star, accepted=TRUE)
  } else {
    return (Theta=Theta.previous, accepted=FALSE)
  }
}

#'
#' Update the regression coefficients related to common
#' covariates and group effects
#'
updateFixedEffectsCovariates <- function(outcomes, covariates, 
                                         X, Theta, Z, alpha, betaCovariates.previous,  
                                         sigma.error, sigma.beta) {
  
  # build the residuals
  
  residuals = vector(0)
  covariatesByGroup = NULL
  for(i in 1:length(X)) {
    residuals.group = outcomes[[i]] - Z[[i]] %*% alpha[[i]] - X[[i]] %*% Theta[[i]]
    residuals <- c(residuals, residuals.group) 
    
    if (is.null(covariatesByGroup)) {
      covariatesByGroup = covariates[[i]]
    } else {
      covariatesByGroup <- rbind(covariatesByGroup, covariates[[i]])
    }
  }
  
  covarIntThetaInverse <- diag(rep(1/sigma.beta, ncol(covariatesByGroup))) 
  
  proposedCovariance <- ginv(covarIntThetaInverse + (1/sigma.beta) * t(covariatesByGroup) %*% covariatesByGroup)
  proposedMean <- proposedCovariance %*% ((1/sigma.beta) * t(covariatesByGroup) %*% residuals)
  
  # Scale Cov to adjust acceptance rate
  proposedCovariance <- proposedCovariance * ifelse(!is.null(mcmc.options$fixedEffectAcceptRateAdjust), 
                                                    mcmc.options$fixedEffectAcceptRateAdjust, 1)
  
  # ensure the covariance is positive definite
  proposedCovariance <- as.matrix(nearPD(proposedCovariance)$mat) 
  
  # draw a proposed set of coefficients for the covariates
  betaCovariates.star <- rmvnorm(1, proposedMean, proposedCovariance)
  
  # calculate the residuals
  resid.star <- residuals - covariatesByGroup %*% t(betaCovariates.star)
  resid.previous <- residuals - covariatesByGroup %*% t(betaCovariates.previous)
  
  # calculate the acceptance probability
  rho<-sum(log(dnorm(resid.star, 0, sqrt(sigma.beta)))) + 
    log(dmvnorm(betaCovariates.star, rep(0, ncol(covariatesByGroup)), covarIntThetaInverse)) + 
    log(dmvnorm(betaCovariates.previous, proposedMean, proposedCovariance)) - 
    sum(log(dnorm(residuals, 0, sqrt(sigma.error)))) - 
    log(dmvnorm(betaCovariates.previous, rep(0, ncol(covariatesByGroup)), covarIntThetaInverse)) - 
    log(dmvnorm(betaCovariates.star, proposedMean, proposedCovariance))
  
  if (rho > log(runif(1))) {
    return (betaCovariates=betaCovariates.star, accepted=TRUE)
  } else {
    return (betaCovariates=betaCovariates.previous, accepted=FALSE)
  }
}


#' Update the random effects
#' 
#' @param 
#' @return 
#' @examples
#'
updateRandomEffects <- function(numSubjects, numObservations, subjectsPerGroup,
                                ids, outcomes, times.observation, 
                                covariates, X, Theta, 
                                Z, alpha, betaCovariates,  
                                sigma.randomIntercept, sigma.randomSlope,
                                sigma.randomInterceptSlope, sigma.error) {
  
  
  # calculate rho
  rho = (sigma.randomInterceptSlope / sqrt(sigma.randomIntercept * sigma.randomSlope))
  
  # some convenience variables
  tau.error = 1 / sigma.error
  tau.randomIntercept = 1 / sigma.randomIntercept
  tau.randomSlope = 1 / sigma.randomSlope
  tau.randomInterceptSlope = 1 / sigma.randomInterceptSlope
  
  residuals = vector(0)
  idsByGroup = NULL
  alphaByGroup = NULL
  timesByGroup = NULL
  for(i in 1:length(X)) {
    residuals.group = (outcomes[[i]] - 
                         Z[[i]][,2] %*% alpha[[i]][,2] - 
                         X[[i]] %*% Theta[[i]])
    if (!is.null(covariates)) {
      residuals.group = residuals.group - covariates[[i]] %*% betaCovariates
    }
    residuals <- c(residuals, residuals.group) 
    
    if (is.null(idsByGroup)) {
      idsByGroup = ids[[i]]
    } else {
      idsByGroup <- rbind(idsByGroup, ids[[i]])
    }
    
    if (is.null(alphaByGroup)) {
      alphaByGroup = alpha[[i]]
    } else {
      alphaByGroup <- rbind(alphaByGroup, alpha[[i]])
    }
    
    if (is.null(timesByGroup)) {
      timesByGroup = times.observation[[i]]
    } else {
      timesByGroup <- rbind(timesByGroup, alpha[[i]])
    }
  }
  
  # get the conditional distribution of the random intercept, given the random slope
  variance.randomIntercept <- 1 / (tau.error * numObservations + tau.randomIntercept * (1 - rho*rho))
  mean.randomIntercept <- (variance.randomIntercept * 
                             (tau.error * tapply(residuals,idsByGroup,sum) + 
                                alphaByGroup[,2] * (rho / (1 - rho*rho)) * 
                                (1 / sqrt(sigma.randomIntercept * sigma.randomSlope))))
  # draw a new random effect sample
  randomIntercepts <-rnorm(numSubjects, mean.randomIntercept, sqrt(variance.randomIntercept))
  
  # get the conditional distribution of the random slope -- FIX RESIDUALS!
  # update the residuals
  residuals = vector(0)
  for(i in 1:length(X)) {
    residuals.group = (outcomes[[i]] - 
                         Z[[i]][,1] %*% alpha[[i]][,1] - 
                         X[[i]] %*% Theta[[i]])
    if (!is.null(covariates)) {
      residuals.group = residuals.group - covariates[[i]] %*% betaCovariates
    }
    residuals <- c(residuals, residuals.group) 
  }
  
  variance.randomSlope <- 1 / (tau.error * tapply(timesByGroup^2, idsByGroup, sum) + tau.randomSlope * (1 - rho*rho))
  mean.randomSlope <- (variance.randomSlope * 
                         (tau.error * tapply(timesByGroup*residuals,patid,sum) + 
                            modelIteration$alpha.intercept * (rho / (1 - rho*rho)) * 
                            (1 / sqrt(sigma.randomIntercept * sigma.randomSlope))))
  
  
  randomSlopes <-rnorm(numSubjects, mean.randomSlope, sqrt(variance.randomSlope))
  
  # build the new random effects matrix
  alpha.total = cbind(randomIntercepts, randomSlopes)
  
  # split up the new random effects into groups again
  start = 1
  alpha = vector(mode = 'list', length(subjectsPerGroup))
  for(i in 1:length(subjectsPerGroup)) {
    alpha[[i]] = alpha.total[start:(start+subjectsPerGroup[[i]])]
    start = start + groupRows[[i]] + 1
  }

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
updateCovarianceParameters <- function(totalObservations, numSubjects,
                                       outcomes, covariates, X, Theta, 
                                       Z, alpha, betaCovariates,
                                       sigma.error,
                                       prior.options) {
  # update the residuals
  residuals = vector(0)
  alphaByGroup = NULL
  for(i in 1:length(X)) {
    residuals.group = (outcomes[[i]] - 
                         Z[[i]] %*% alpha[[i]] - 
                         X[[i]] %*% Theta[[i]])
    if (!is.null(covariates)) {
      residuals.group = residuals.group - covariates[[i]] %*% betaCovariates
    }
    residuals <- c(residuals, residuals.group) 
    
    if (is.null(alphaByGroup)) {
      alphaByGroup = alpha[[i]]
    } else {
      alphaByGroup <- rbind(alphaByGroup, alpha[[i]])
    }
  }

  # sample from an inverse Gamma to update sigma.error
  shape <- prior.options$shape.tau + (N / 2) 
  rate <- prior.options$rate.tau + (crossprod(residuals) / 2)
  sigma.error <- 1 / rgamma(1, shape, rate)
  
  # sample from an inverse wishart to update the covariance of the random effects
  sigma.alpha <- riwish(prior.options$sigmaError.df + numSubjects, 
                        prior.options$sigmaError.scaleMatrix + crossprod(as.matrix(alphaByGroup)))
  
  
  sigma.randomIntercept = sigma.alpha[1,1]
  sigma.randomSlope = sigma.alpha[2,2]
  sigma.randomInterceptSlope = sigma.alpha[1,2]
  
  return (list(sigma.randomIntercept = sigma.alpha[1,1],
               sigma.randomSlope = sigma.alpha[2,2],
               sigma.randomInterceptSlope = sigma.alpha[1,2]))
}

#
#
#

#'
#' Fit a varying coefficient model for longitudinal studies with
#' informative dropout.  Uses a Bayesian approach with a spline fit
#' to model the relationship between dropout times and slope
#'
#' 
#'
#'
#'
informativeDropout.bayes.splines <- function(data, ids.var, outcomes, groups, covariates, 
                                             times.dropout, times.observation, dist,
                                             knots.options, mcmc.options, prior.options) {
  
  # validate the knot options.
  if (knots.options$birthProbability <= 0 || knots.options$birthProbability >= 1) {
    stop("Knot options error :: Invalid birth probability. Please specified a value between 0 and 1.")
  }
  if (is.na(knots.options$min) || knots.options$min <= 0) {
    stop("Knot options error :: The minimum number of knots must be 1 or greater.")
  }
  if (is.na(knots.options$max) || knots.options$max <= knots.options$min) {
    stop("Knot options error :: The maximum number of knots must be greater than the minimum number of knots.")
  }
  # create some reasonable defaults for the start positions and candidate positions
  # if not specified
  if (is.null(knots.options$startPositions) || is.null(knots.options$candidatePositions)) {
    dropout.min = min(data[,times.dropout])
    dropout.max = max(data[,times.dropout])
    if (is.null(knots.options$candidatePositions)) {
      knots.options$candidatePositions = seq(dropout.min, dropout.max, 1)
    }
    if (is.null(knots.options$startPositions)) {
      knots.options$startPositions = sample(knots.options$candidatePositions, knots.options$min)
    }
  }
  
  # validate the mcmc options
  if (is.na(mcmc.options$iterations) || mcmc.options$iterations <= 1) {
    stop("RJMCMC options error :: invalid number of iterations")
  }
  if (is.na(mcmc.options$burnIn) || mcmc.options$burnIn >= mcmc.options$iterations) {
    stop("RJMCMC options error :: the burn in period must be less than the total number of iterations")
  }
  
  # validate the prior options
  if (is.na(prior.options$shape.tau) || prior.options$shape.tau <= 0) {
    stop("Prior options error :: shape.tau must be greater than 0")
  }
  if (is.na(prior.options$rate.tau) || prior.options$rate.tau <= 0) {
    stop("Prior options error :: rate.tau must be greater than 0")
  }
  if (is.na(prior.options$sigmaError.df) || prior.options$sigmaError.df <= 0) {
    stop("Prior options error :: sigmaError.df must be greater than 0")
  }
  if (is.null(prior.options$sigmaError.scaleMatrix) || 
      !is.positive.definite(prior.options$sigmaError.scaleMatrix)) {
    stop("Prior options error :: sigmaError.scaleMatrix must be a positive definite matrix")
  }
  
  ## TODO: finish validation of other input parameters
  
  # get the list of treatment groups
  groupList <- unique(data[,groups])
  
  # cache a few values we will use repeatedly
  cache = list(
    # number of subjects
    numSubjects = length(unique(data[,ids.var])),
    # number of observations per subject
    numObservations = as.vector(table(data[,ids.var])),
    # number of subjects per group
    subjectsPerGroup = lapply(groupList, function(group) { 
      # get the covariates
      return(length(unique(data[data[,groups] == group, ids.var])))
    }),
    # participant ids by group
    ids = lapply(groupList, function(group) { 
      # get the covariates
      return(data[data[,groups] == group, ids.var])
    }),
    # outcomes ordered by group
    outcomes = lapply(groupList, function(group) { 
      # get the covariates
      return(data[data[,groups] == group, outcomes])
    }),
    # observation times ordered by group
    times.observation = lapply(groupList, function(group) { 
      # get the covariates
      return(data[data[,groups] == group, times.observation])
    }),
    # data split into lists by group
    covariates = lapply(groupList, function(group) { 
      # get the covariates
      covar = data[data[,groups] == group, covariates]
      # add columns for the group beta coefficients
      for(nonreferenceGroup in groupList[2:length(groupList)]) {
        covar[[groupIndicatorColumn(nonreferenceGroup)]] = 
          as.numeric(data[data[,groups] == group,groups] == nonreferenceGroup)
      }
      return (covar)
    }),
    
    Z = lapply(groupList, function(group) { 
      # get the times
      groupTimes = data[data[,groups] == group, times.observation]
      # add a column of 1's for the intercept
      return (cbind(rep(1,length(groupTimes)), groupTimes))
    }),
    
    # random effects
    alpha = lapply(groupList, function(group) { 
      groupIds = data[data[,groups] == group, ids.var]
      return (matrix(rep(0, 2 * length(unique(groupIds)))))
    }),
    
    # X matrix from the previous iteration, split by group
    X = lapply(groupList, function(group) { 
      groupTimes = data[data[,groups] == group,times.observation]
      groupDropout = data[data[,groups] == group,times.dropout]
      knots.boundary = range(knots.options$startPositions)
      knots.interior = knots.options$startPositions[-c(1,length(knots.options$startPositions))] 
      return(as.matrix(cbind(
        rep(1,length(groupTimes)),
        ns(groupDropout, knots=knots.interior, Boundary.knots=knots.boundary, intercept=T) * groupTimes
      )))
    })
  )
  
  # get estimates for the initial regression coefficients
  Xfull = NULL
  CovarFull = NULL
  OutcomesFull = NULL
  for(i in 1:length(cache$X)) {
    if (is.null(Xfull)) {
      Xfull = cache$X[[i]]
    } else {
      Xfull <- rbind(Xfull, cache$X[[i]])
    }
    
    if (is.null(CovarFull)) {
      CovarFull = cache$covariates[[i]]
    } else {
      CovarFull <- rbind(CovarFull, cache$covariates[[i]])
    }
    
    if (is.null(OutcomesFull)) {
      OutcomesFull = cache$outcomes[[i]]
    } else {
      OutcomesFull <- c(OutcomesFull, cache$outcomes[[i]])
    }
  }
  OutcomesFull <- as.data.frame(OutcomesFull)
  names(OutcomesFull) <- outcomes
  Xfull = as.data.frame(Xfull)
  names(Xfull) <- sapply(0:(ncol(Xfull)-1), function(i) { return(paste("theta", i, sep=''))})
  data.theta = cbind(OutcomesFull, Xfull)
  fit.Theta <- lm(as.formula(paste(c(paste(outcomes, "~"), paste(names(Xfull)[2:length(names(Xfull))], 
                                                                 collapse=" + ")), collapse=" ")), 
                  data=data.theta)
  
  
  # initialize the first model iteration
  # TODO: mcmc options specify starting values
  modelIterationList <- vector(mode = "list", length = mcmc.options$iterations)
  modelIterationList[[1]] = rjmcmc.iteration(knots=knots.options$startPositions, 
                                             Theta=NULL, betaCovariate=NULL,
                                             sigma.error = 1, sigma.randomIntercept = 1,
                                             sigma.randomSlope = 1, sigma.randomInterceptSlope = 0,
                                             lambda.numKnots = 1, sigma.beta = 1,
                                             shape.tau = 0.001, rate.tau = 0.001)
  
  #
  # Run the reversible jump MCMC
  #
  for (i in 2:mcmc.options$iterations) {
    
    model.previous = modelIterationList[i-1]
    # make a copy which will be modified as we move through the iteration
    model.current = model.previous
    
    for (group.index in 1:length(groupList)) {
      group = groupList[group.index]
      # get the subset of data for this group
      groupData = data[data[,groups] == group]
      # randomly decide to add/remove a knot
      if (runif(1) < knots.options$birthProbability) {
        # add a knot
        result = addKnot(model.current$knots, knots.options, 
                         groupData[,outcomes], groupData[,times.dropout], 
                         groupData[,times.observation], cache$covariates[[group.index]],
                         cache$X[[group.index]], model.current$Theta[[group.index]],
                         cache$Z[[group.index]], cache$alpha[[group.index]], model.current$betaCovariates,  
                         model.current$sigma.error, model.current$sigma.beta, model.current$lambda)
        # update the model iteration
        model.current$X[[group.index]] = result$X
        model.current$knots = result$knots
        model.current$Theta[[group.index]] = result$Theta
        model.current$proposed$knot.add = TRUE
        model.current$accepted$knot.add = result$accepted
        
      } else {
        # remove a knot
        result = removeKnot(model.current$knots, knots.options, 
                            groupData[,outcomes], groupData[,times.dropout], 
                            groupData[,times.observation], cache$covariates[[group.index]],
                            cache$X[[group.index]], model.current$Theta[[group.index]],
                            cache$Z[[group.index]], cache$alpha[[group.index]], model.current$betaCovariates,  
                            model.current$sigma.error, model.current$sigma.beta, model.current$lambda)  
        
        # update the model iteration
        model.current$X[[group.index]] = result$X
        model.current$knots = result$knots
        model.current$Theta[[group.index]] = result$Theta
        model.current$proposed$knot.remove = TRUE
        model.current$accepted$knot.remove = result$accepted
      }
      
      # Move knots
      result = moveKnot(model.current$knots, knots.options$stepSize[[group.index]], 
                        knots.options$candidatePositions,
                        groupData[,outcomes], groupData[,times.dropout], 
                        groupData[,times.observation], cache$covariates[[group.index]],
                        cache$X[[group.index]], model.current$Theta[[group.index]],
                        cache$Z[[group.index]], cache$alpha[[group.index]], model.current$betaCovariates,  
                        model.current$sigma.error)
      model.current$X[[group.index]] = result$X
      model.current$knots = result$knots
      model.current$proposed$knot.move = result$proposed
      model.current$accepted$knot.move = result$accepted
      
      # update fixed effects (includes coefficients for covariates and time varying slopes)
      result = updateFixedEffects(model.current$knots, knots.options, 
                                  groupData[,outcomes], groupData[,times.dropout], 
                                  groupData[,times.observation], cache$covariates[[group.index]],
                                  cache$X[[group.index]], model.current$Theta[[group.index]],
                                  cache$Z[[group.index]], cache$alpha[[group.index]], model.current$betaCovariates,  
                                  model.current$sigma.error, model.current$sigma.beta, model.current$lambda)
      model.current$Theta[[group.index]] = result$Theta
      model.current$proposed$fixedEffects = TRUE
      model.current$accepted$fixedEffects = result$accepted
      
    }  
    
    # update fixed effects associated with covariates
    result = updateFixedEffectsCovariates(cache$outcomes, cache$covariates, 
                                          cache$X, model.current$Theta, 
                                          cache$Z, cache$alpha, model.current$betaCovariates,  
                                          model.current$sigma.error, model.current$sigma.beta)
    model.current$betaCovariates = result$betaCovariates
    model.current$proposed$fixedEffectsCovariates = TRUE
    model.current$accepted$fixedEffectsCovariates = result$accepted

    
    # update random effects
    cache$alpha = updateRandomEffects(cache$numSubjects, cache$numObservations, cache$subjectsPerGroup,
                                      cache$ids, cache$outcomes, cache$times.observations,
                                      cache$covariates, cache$X, model.current$Theta, 
                                      cache$Z, cache$alpha, model.current$betaCovariates,
                                      model.current$sigma.randomIntercept, 
                                      model.current$sigma.randomSlope,
                                      model.current$sigma.randomInterceptSlope,
                                      model.current$sigma.error)
    
    # update variance components
    result = updateVarianceComponents(modelIteration.current, modelIteration.previous, knots.options)
    model.current$sigma.randomIntercept = result$sigma.randomIntercept
    model.current$sigma.randomSlope = result$sigma.randomSlope
    model.current$sigma.randomInterceptSlope = result$sigma.randomInterceptSlope
    
    # save the current iteration
    modelIterationList[[i]] = model.current
  }
  
  # calculate the final estimates as the mean across the different iterations
  
  
  # return the estimates, with distributions, and the model results from each iteration
  return (model.fit)
  
}

#' Fit a varying coefficient model for longitudinal studies with
#' informative dropout. 
#' 
#' @param 
#' @param 
#' @return 
#' @examples
#' 
informativeDropout <- function(data, id.var, outcomes.var, group.var, covariates.var, 
                               times.dropout.var, times.observation.var,
                               method="bayes.splines", dist="normal",
                               knots.options=list(birthProbability=0.5, min=3, max=10, 
                                                  startPositions=NULL, candidatePositions=NULL), 
                               mcmc.options=list(iterations=100000, burnIn=50000),
                               prior.options=list(shape.tau = 0.001, rate.tau = 0.001, 
                                                  sigmaError.df = 3, sigmaError.scaleMatrix = diag(2))) {
  
  if (method == 'bayes.splines') {
    # model the relationship between dropout time and slope using natural splines
    return (informativeDropout.bayes.splines(data, id.var, outcomes.var, group.var, covariates.var, 
                                             times.dropout.var, times.observation.var, dist, 
                                             knots.options, mcmc.options, prior.options))
  } else if (method == 'bayes.dirichlet') {
    # account for informative dropout using a dirichlet process 
    return (informativeDropout.bayes.dirichlet(data, id.var, outcomes.var, group.var, covariates.var, 
                                               times.dropout.var, times.observation.var, dist, prior.options))
  } else if (method == 'mixed') {
    # fit a mixed model which models the relationship between dropout time and slope using natural splines
    return (informativeDropout.mixed(data, id.var, outcomes.var, group.var, covariates.var, 
                                     times.dropout.var, times.observation.var, dist, dist))
  }
  
}

test.example <- function() {
  data <- read.csv("test.csv")
  data$day = data$years * 365
  
  informativeDropout(data, "WIHSID", "logcd4", "hard", c("AGEATBL", "minority"), 
                     "drop", "day", "bayes.splines", "normal",
                     knots.options=list(birthProbability=0.5, min=3, max=10, 
                                        startPositions=c(330,550,1060), candidatePositions=seq(10,max(data$day),10)), 
                     mcmc.options=list(iterations=20, burnIn=10),
                     prior.options=list(shape.tau = 0.001, rate.tau = 0.001, 
                                        sigmaError.df = 3, sigmaError.scaleMatrix = diag(2)))
}


# for debugging
data = data
data = data
ids.var = "WIHSID"
outcomes = "logcd4"
groups = "hard"
covariates = c("AGEATBL", "minority")
times.dropout = "drop"
times.observation = "day"
dist = "normal"
knots.options=list(birthProbability=0.5, min=3, max=10, 
                   startPositions=c(330,550,1060), candidatePositions=seq(10,max(data$day),10)) 
mcmc.options=list(iterations=20, burnIn=10)
prior.options=list(shape.tau = 0.001, rate.tau = 0.001, 
                   sigmaError.df = 3, sigmaError.scaleMatrix = diag(2))



