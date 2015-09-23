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

#' A class containing the current state of the model
#' during a RJMCMC run
#'
modelIteration <- function(knots=NULL, X=NULL, ) {
  mi = list(
    knots = knots,
    beta0=NULL,
    Theta=NULL,
    knotAdded = FALSE,
    betaCovariate = NULL,
    proposed=list(knot.add=FALSE, knot.remove=FALSE, knot.move=FALSE,
                  fixedEffects=FALSE, randomEffects=FALSE, varianceComponents=FALSE),
    accepted=list(knot.add=FALSE, knot.remove=FALSE, knot.move=FALSE,
                    fixedEffects=FALSE, randomEffects=FALSE, varianceComponents=FALSE)
    
  )
  
  class(mi) <- append(class(mi), "modelIteration")
  return(mi)
  
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
addKnot <- function(data, group, outcomes, treatment, covariates, dropoutTimes, 
                    times, modelIteration.previous, X.prev, knots.options,
                    prior.options) {
  
  #   ### remove debugging
  #   data=dat[data$hard==1,]
  #   group="1"
  #   outcomes="logcd4"
  #   treatment="hard"
  #   covariates=NA
  #   dropoutTimes="drop"
  #   times = "years"
  #   
  #   knots.options=list(
  #     candidatePositions = candidates.g1,
  #     
  #     min=3, max=10
  #   )
  #   
  #   modelIteration.previous = list(
  #     knots = currentknots.g1,
  #     alpha = c(0, 1),
  #     beta.covar = NA,
  #     Theta = list("1"=(Theta.LSXprev + rnorm(1)))
  #   )
  #   
  #   prior.options = list(sigmaError=1.25^2)
  #   
  #   knots.boundary = range(modelIteration.previous$knots)
  #   knots.interior = modelIteration.previous$knots[-c(1,length(modelIteration.previous$knots))] 
  #   
  #   Xprev = cbind(
  #     rep(1,nrow(data)),
  #     ns(data[,dropoutTimes], knots=knots.interior, Boundary.knots=knots.boundary, intercept=T) * data[,times]
  #   )
  #   ### END DEBUG
  #   
  
  # add a knot by randomly selecting a candidate knot
  index = sample(1:length(knots.options$candidatePositions), 1)
  newKnot.value = knots.options$candidatePositions[index]
  knots <- sort(c(modelIteration.previous$knots, newKnot.value))
  # get the interior and boundary knots, and grab the position of the knot that
  # was just added
  knots.boundary = range(knots)
  knots.interior = knots[-c(1,length(knots))] 
  newKnot.position = which(knots == newKnot.value)
  
  # Calculate spline transformation of dropout time and create the proposed X matrix
  Xstar <- cbind(
    rep(1,nrow(data)),
    ns(data[,dropoutTimes], knots=knots.interior, Boundary.knots=knots.boundary, intercept=T) * data[,times]
  )
  
  # Calculate y-random effects for least squares calculations
  y = as.matrix(data[,outcomes])
  Z = cbind(rep(1,nrow(data)), data[,times])
  if (!is.na(covariates)) {
    yls <-as.vector(y - Z %*% modelIteration.previous$alpha - 
                      as.matrix(data[,covariates]) %*% modelIteration.previous$betaC)
  } else {
    yls <-as.vector(y - Z %*% modelIteration.previous$alpha)
  }
  
  
  # Calculate least squares estimates for coefficients and differences between LS and current coefficients
  Theta.LSXprev <- ginv(crossprod(Xprev))%*%(crossprod(Xprev,yls))
  Theta.LSXstar <- ginv(crossprod(Xstar))%*%(crossprod(Xstar,yls))
  Theta.LSresid <- modelIteration.previous$Theta[[group]] - Theta.LSXprev
  #Draw a residual for the added coefficient and calculate coefficient transformation
  residual <- rnorm(1, 0, sd.resid)
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
  
  
  #Calculate residuals for likelihood ratio
  if (!is.na(covariates)) {
    LRresid.star <- as.vector(y - Xstar %*% Theta.star - 
                                as.matrix(data[,covariates]) %*% modelIteration.previous$betaC -
                                Z %*% modelIteration.previous$alpha
    )
    LRresid.prev <- as.vector(y - XPrev %*% modelIteration.previous$Theta - 
                                as.matrix(data[,covariates]) %*% modelIteration.previous$betaC -
                                Z %*% modelIteration.previous$alpha
    )
  } else {
    LRresid.star <-as.vector(y - Xstar %*% Theta.star - Z %*% modelIteration.previous$alpha)
    LRresid.prev <- as.vector(y - XPrev %*% modelIteration.previous$Theta - Z %*% modelIteration.previous$alpha)
  }
  
  # Calculate birth and death probabilities                                                        
  probBirth <- ifelse(length(modelIteration.previous$knots) == knots.options$min, 1, knots.options$birthProbability)
  probDeath <- ifelse(length(modelIteration.previous$knots) == knots.options$max - 1, 1, 1 - knots.options$birthProbability)
  
  #Calculate Acceptance Probability                                                        
  rho <- (log(modelIteration.previous$prior.lambda) - 
            log(length(modelIteration.previous$knots)) + 
            log(probDeath) - log(probBirth) + 
            log(sqrt(modelIteration.previous$prior.sigmaE)) - 
            0.5 * log(modelIteration.previous$sigmaBeta) +
            (residual^2 / (2 * modelIteration.previous$prior.sigmaE * modelIteration.previous$prior.sigmaE)) + 
            ((crossprod(modelIteration.previous$Theta) - crossprod(Theta.star)) / (2 * modelIteration.previous$sigmaBeta)) +
            ((crossprod(LRresid.prev) - crossprod(LRresid.star)) / (2 * modelIteration.previous$sigmaBeta))
  )
  
  if (rho > log(runif(1))) {
    return (Xstar, knots, Theta.star, true)
  } else {
    return (Xprev, modelIteraton.previous$knots, modelIteration.previous$Theta[[group]], false)          
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
removeKnot <- function(data, group, outcomes, treatment, covariates, dropoutTimes, 
                       times, modelIteration.previous, X.prev, knots.options,
                       prior.options) {
  # randomly remove an existing knot
  index = sample(1:length(modelIteration.previous$knots), 1)
  knots <- modelIteration.previous$knots[-index]
  
  knots.boundary = range(knots)
  knots.interior = knots[-c(1,length(knots))] 
  
  if (length(knots) > 1 ) {  
    Xstar<-cbind(
      rep(1, nrow(data)),
      ns(data[,dropoutTimes], knots=knotstar.i, Boundary.knots=knotstar.b, intercept=T)*t[group==1]
    )
  } else {
    Xstar<-cbind(rep(1, nrow(data)), data[,times])
  }
  
  # Calculate y-random effects for least squares calculations
  y = as.matrix(data[,outcomes])
  Z = cbind(rep(1,nrow(data)), data[,times])
  if (!is.na(covariates)) {
    yls <-as.vector(y - Z %*% modelIteration.previous$alpha - 
                      as.matrix(data[,covariates]) %*% modelIteration.previous$betaC)
  } else {
    yls <-as.vector(y - Z %*% modelIteration.previous$alpha)
  }
  
  # Calculate least squares estimates for coefficients and differences between LS and current coefficients
  Theta.LSXprev <- ginv(crossprod(Xprev))%*%(crossprod(Xprev,yls))
  Theta.LSXstar <- ginv(crossprod(Xstar))%*%(crossprod(Xstar,yls))
  Theta.LSresid <- modelIteration.previous$Theta[[group]] - Theta.LSXprev
  # update the coefficients
  residual.deletedKnot <- Theta.LSresid[index]
  Theta.star <- Theta.LSXstar + Theta.LSresid[-index]
  
  
  #Calculate residuals for likelihood ratio
  if (!is.na(covariates)) {
    LRresid.star <- as.vector(y - Xstar %*% Theta.star - 
                                as.matrix(data[,covariates]) %*% modelIteration.previous$betaC -
                                Z %*% modelIteration.previous$alpha
    )
    LRresid.prev <- as.vector(y - XPrev %*% modelIteration.previous$Theta - 
                                as.matrix(data[,covariates]) %*% modelIteration.previous$betaC -
                                Z %*% modelIteration.previous$alpha
    )
  } else {
    LRresid.star <-as.vector(y - Xstar %*% Theta.star - Z %*% modelIteration.previous$alpha)
    LRresid.prev <- as.vector(y - XPrev %*% modelIteration.previous$Theta - Z %*% modelIteration.previous$alpha)
  }
  
  # Calculate birth and death probabilities                                                        
  probBirth <- ifelse(length(modelIteration.previous$knots) == knots.options$min, 1, knots.options$birthProbability)
  probDeath <- ifelse(length(modelIteration.previous$knots) == knots.options$max - 1, 1, 1 - knots.options$birthProbability)
  
  #Calculate Acceptance Probability                                                        
  rho <- (log(modelIteration.previous$prior.lambda) - 
            log(length(modelIteration.previous$knots)) + 
            log(probDeath) - log(probBirth) + 
            log(sqrt(modelIteration.previous$prior.sigmaE)) - 
            0.5 * log(modelIteration.previous$sigmaBeta) +
            (residual^2 / (2 * modelIteration.previous$prior.sigmaE * modelIteration.previous$prior.sigmaE)) + 
            ((crossprod(modelIteration.previous$Theta) - crossprod(Theta.star)) / (2 * modelIteration.previous$sigmaBeta)) +
            ((crossprod(LRresid.prev) - crossprod(LRresid.star)) / (2 * modelIteration.previous$sigmaBeta))
  )
  
  if (1/rho > log(runif(1))) {
    return (Xstar, knots, Theta.star, true)
  } else {
    return (Xprev, modelIteraton.previous$knots, modelIteration.previous$Theta[[group]], false)          
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
moveKnot <- function(data, group, outcomes, treatment, covariates, dropoutTimes, 
                     times, modelIteration.current, X.prev, knots.options,
                     prior.options) {
  
  knots <- modelIteration.current$knots[[group]]
  
  #Pick a knot to move 
  knotToMove <- sample(knots, 1) 
  # get index of knot to be moved
  index <- which(knots == knotToMove) 
  # get the knots that are staying in the same place
  knotsToKeep <- knots[-index] 
  
  # find a new location from the potential knot locations
  potentialLocations <- knots.options$candidatePositions[!(knots.options$candidatePositions %in% knots)]
  # here we only allow movement within some small window
  potentialLocations <- potentialLocations[potentialLocations > (index - knots.options$stepSize[[group]]) & 
                                             potentialLocations < (index + knots.options$stepSize[[group]])] 
  
  
  if (length(potentialLocations > 0)) {
    # pick a new location
    p <- sample(potentialLocations,1)
    # get the new knots
    knotstar <- sort(c(knotsToKeep, p))
    
    if (length(knotstar) > 1) {
      knotstar.b<-range(knotstar)
      knotstar.i<-knotstar[!(knotstar %in% knotstar.b)]
      # Calculate spline transformation of dropout time and proposed X matrix and residuals
      Xstar<-cbind(
        rep(1, nrow(data)),
        ns(data[,dropoutTimes], knots=knotstar.i, Boundary.knots=knotstar.b, intercept=T)*t[group==1]
      )
    } else {
      Xstar<-cbind(rep(1, nrow(data)), data[,times])
    }
    
    #Calculate residuals for likelihood ratio
    if (!is.na(covariates)) {
      LRresid.star <- as.vector(y - Xstar %*% modelIteration.current$Theta - 
                                  as.matrix(data[,covariates]) %*% modelIteration.current$betaC -
                                  Z %*% modelIteration.current$alpha
      )
      LRresid.prev <- as.vector(y - XPrev %*% modelIteration.current$Theta - 
                                  as.matrix(data[,covariates]) %*% modelIteration.current$betaC -
                                  Z %*% modelIteration.current$alpha
      )
    } else {
      LRresid.star <-as.vector(y - Xstar %*% Theta.star - Z %*% modelIteration.current$alpha)
      LRresid.prev <- as.vector(y - XPrev %*% modelIteration.current$Theta - Z %*% modelIteration.current$alpha)
    }
    
    #Calculate the acceptance probability
    rho <- log(length(potentialknots)) - log(length(pknotstar)) + 
      (crossprod(LRresid.prev)-crossprod(LRresid.star))/2/modelIteration.current$prior.sigmaE
    if (rho > log(runif(1))) {
      return (knotstar, Xstar, true, true)
    } else {
      return (knots, Xprev, true, false)
    }
  } else {
    # nowhere to move to
    return (knots, Xprev, false, false)
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
updateFixedEffects <- function(data, group, outcomes, treatment, covariates, dropoutTimes, 
                               times, modelIteration.current, X.prev, knots.options,
                               prior.options) {
  
  knots <- modelIteration.current$knots[[group]]
  sigmaBeta <- modelIteration.current$prior.sigmaBeta
  Xprev <- X.prev[[group]]
  ThetaPrev <- modelIteration.current$Theta[[group]]
  
  # create the covariance matrix for the intercept and theta coefficients for the splines - R0
  covarIntTheta <- diag(rep(sigmaBeta, (length(knots)+1))) 
  covarIntThetaInverse <- diag(rep(1/sigmaBeta, (length(knots)+1))) 
  # Calculate residuals 
  if (!is.na(covariates)) {
    LRresid <- as.vector(y - X.prev %*% modelIteration.current$Theta - 
                                as.matrix(data[,covariates]) %*% modelIteration.current$betaC -
                                Z %*% modelIteration.current$alpha
    )
  } else {
    LRresid <-as.vector(y - X.prev %*% Theta.star - Z %*% modelIteration.current$alpha)
  }
  
  # calculate the covariance and mean of the proposal distribution
  proposedCovariance <- ginv(covarIntThetaInverse + (1/sigmaBeta) * crossprod(Xprev))
  proposedMean <- proposedCovariance %*% ((1/sigmaBeta) * crossprod(Xprev, LRresid))

  # Scale Cov to adjust acceptance rate
  proposedCovariance <- mcmc.options$fixedEffectAcceptRateAdjust * proposedCovariance 
  
  # ensure the covariance is positive definite
  proposedCovariance <- as.matrix(nearPD(proposedCovariance)$mat) 
  # draw a proposed set of coefficients
  ThetaStar <- rmvnorm(1, proposedMean, proposedCovariance) 
  
  # Calculate residuals for likelihood ratio
  residStar <- LRresid - crossprod(Xprev, ThetaStar)
  residPrev <- LRresid - crossprod(Xprev, ThetaPrev)

  # Calculate acceptance probability
  rho <- log(dmvnorm(ThetaPrev, proposedMean, proposedCovariance)) - 
    log(dmvnorm(ThetaStar, proposedMean, proposedCovariance)) + 
    (crossprod(ThetaPrev)-tcrossprod(ThetaStar))/(2 * sigmabeta) + 
    (crossprod(residPrev)-crossprod(residStar))/(2 * sigmaBeta)
  
  if (rho > log(runif(1))) {
    return (ThetaStar, TRUE, TRUE)
  } else {
    return (ThetaPrev, TRUE, FALSE)
  }
}

#'
#'
#'
#'
#'
#'
#'
updateFixedEffectsCovariates <- function(data, outcomes, treatment, covariates, dropoutTimes, 
                                         times, modelIteration.current, X.prev, knots.options,
                                         prior.options) {
  sigmaBeta <- modelIteration.current$prior.sigmaBeta
  
  # subtract the fixed effects for each group
  groupList <- unique(data[,treatment])
  
  residuals = vector(0)
  covariatesByGroup = NULL
  for(group in groupList) {
    dataGroup <- data[data[,treatment] == group,]
    yGroup <- dataGroup[,treatment]
    covarGroup <- dataGroup[,covariates]
    Z = cbind(rep(1,nrow(dataGroup)), dataGroup[,times])
    
    residualsGroup <- yGroup - Z %*% modelIteration.previous$alpha - X.prev[[group]] %*% ThetaPrev
    residuals <- c(residuals, residualsGroup) 
    if (is.null(covarByGroup)) {
      covariatesByGroup = covarGroup
    } else {
      covariatesByGroup <- rbind(covariatesByGroup, covarGroup)
    }
  }
  
  covarIntThetaInverse <- diag(rep(1/sigmaBeta, length(covariates))) 
  
  proposedCovariance <- ginv(covarIntThetaInverse + (1/sigmaBeta) * t(covariatesByGroup) %*% covariatesByGroup)
  proposedMean <- proposedCovariance %*% ((1/sigmaBeta) * t(covariatesByGroup) %*% residuals)
  
  # Scale Cov to adjust acceptance rate
  proposedCovariance <- mcmc.options$fixedEffectAcceptRateAdjust * proposedCovariance 
  
  # ensure the covariance is positive definite
  proposedCovariance <- as.matrix(nearPD(proposedCovariance)$mat) 

  # draw a proposed set of coefficients for the covariates
  betaCovariateStar <- rmvnorm(1, proposedMean, proposedCovariance)
  
  # calculate the residuals
  residualStar <- residuals - covariatesByGroup %*% t(betaCovariateStar)
  residualPrev <- residuals - covariatesByGroup %*% t(modelIteration.current$betaCovariate)
  
  # calculate the acceptance probability
  rho<-sum(log(dnorm(residualStar, 0, sqrt(sigmaBeta)))) + 
    log(dmvnorm(betaCovariateStar, rep(0, length(covariates)), covarIntThetaInverse)) + 
    log(dmvnorm(modelIteration.current$betaCovariate, proposedMean, proposedCovariance)) - 
    sum(log(dnorm(yt, 0, sqrt(sig[i-1])))) - 
    log(dmvnorm(modelIteration.current$betaCovariate, rep(0, length(covariates)), covarIntThetaInverse)) - 
    log(dmvnorm(betaCovariateStar, proposedMean, proposedCovariance))
  
  if (rho > log(runif(1))) {
    return (betaCovariateStar, TRUE, TRUE)
  } else {
    return (modelIteration.current$betaCovariate, TRUE, FALSE)
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
updateRandomEffects <- function(data, group, outcomes, treatment, covariates, dropoutTimes, 
                                times, modelIteration.previous, X.prev, knots.options,
                                prior.options) {
  # Update b2|b1
  
  vb<-1/(tau21+tau*tapply(t^2,patid,sum))
  mb0<-rhob*sqrt(sigmab2/sigmab1)*B1     # Conditional Prior Mean of b2|b1
  yt<-y-rep(B1,nobs)-bt[i]-C%*%covar[i,]
  yt[group==1]<-yt[group==1]-Xt.g1%*%as.vector(ct.g1[[i]])
  yt[group==2]<-yt[group==2]-hardeffect.g1[i]*hard[group==2]-Xt.g2%*%as.vector(ct.g2[[i]])
  mb<-vb*(tau21*mb0+tau*tapply(t*(yt),patid,sum)) 
  B2<-rnorm(nsub,mb,sqrt(vb))
  
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
updateCovarianceParameters <- function(data, group, outcomes, treatment, covariates, dropoutTimes, 
                                       times, modelIteration.previous, X.prev, knots.options,
                                       prior.options) {
  # Update tau  
  # get the current residual
  if (!is.na(covariates)) {
    yt <-as.vector(y - Z %*% modelIteration.previous$alpha - 
                      as.matrix(data[,covariates]) %*% modelIteration.previous$betaC)
  } else {
    yt <-as.vector(y - Z %*% modelIteration.previous$alpha)
  }
  
  
  yt<-y-rep(B1,nobs)-t*rep(B2,nobs)-bt[i]-C%*%covar[i,]
  yt[group==1]<-yt[group==1]-Xt.g1%*%as.vector(ct.g1[[i]])
  yt[group==2]<-yt[group==2]-Xt.g2%*%as.vector(ct.g2[[i]])-hardeffect.g1[i]*hard[group==2]
  g<-g0+crossprod(yt)/2
  taus[i]<-tau<-rgamma(1,d0+N/2,g)
  sig[i]<-1/taus[i]
  
  # Update Sigma.b  
  bmat<-cbind(B1,B2)            
  Sigma.b<-riwish(nu0+nsub,c0+crossprod(bmat))
  sigmab1<-Sigma.b[1,1]                    # Marginal variance of b1
  sigmab2<-Sigma.b[2,2]                    # Marginal variance of b2
  rhob<-Sigma.b[1,2]/sqrt(sigmab1*sigmab2) # Corr(b1,b2)
  tau12<-1/(sigmab1*(1-rhob^2))            # Conditional precision of b1|b2
  tau21<-1/(sigmab2*(1-rhob^2))            # Conditional precision of b2|b1
  Sigma.bs[i,]<-c(Sigma.b)   
}

#
#
#

#'
#' Internal function which uses the Bayesian approach to fit
#' varying coefficient models for Gaussian outcomes 
#'
#'
#'
#'
#'
_informativeDropout.bayesian.gaussian <- function(data, treatment, covariates, dropoutTimes, 
                                                  times, knots.options, mcmc.options) {
  
  
  groupList <- unique(data$treatment)
  
  modelIterationList <- list(mcmc.options$iterations)
  
  knotPositions <- knots.options$startPositions
  
  for (i in 2:mcmc.options$iterations) {
    
    modelIteration.previous = modelIterationList[i-1]
    
    for (group in groupList) {
      groupData = data[data[,treatment] == group]
      # randomly decide to add/remove a knot
      if (runif(1) < knots.options$birthProbability && ) {
        modelIteration.current = addKnot(modelIteration.previous, knots.options)
      } else {
        modelIteration.current = removeKnot(modelIteration.previous, knots.options)                 
      }
      
      # Move knots
      modelIteration.current = moveKnots(modelIteration.current, modelIteration.previous, knots.options)
      
      # update fixed effects (includes coefficients for covariates and time varying slopes)
      modelIteration.current = updateFixedEffects(modelIteration.current, modelIteration.previous, knots.options)
      
      
    }  
    
    # update fixed effects associated with covariates
    modelIteration.current = updateFixedEffectsCovariates(modelIteration.current, modelIteration.previous, knots.options)
    
    # update random effects
    modelIteration.current = updateRandomEffects(modelIteration.current, modelIteration.previous, knots.options)
    
    # update variance components
    modelIteration.current = updateVarianceComponents(modelIteration.current, modelIteration.previous, knots.options)
  }
  
  # calculate the final estimates as the mean across the different iterations
  
  
  # return the estimates, with distributions, and the model results from each iteration
  return model.fit
  
}

#' Add together two numbers.
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
informativeDropout <- function(data, treatment, covariates, dropoutTimes, times,
                               method="bayesian", dist="normal",
                               knots.options=list(birthProbability=0.5, min=3, max=10, 
                                                  startPositions=NULL, candidatePositions=NULL), 
                               mcmc.options=(iterations=100000,burnIn=50000)) {
  
  #
  
  # return an informativeDropout.fit class
  # split into either 
}


# functions
# - summary: show estimates
# - print
# - predict: get distribution for 


