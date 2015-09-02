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
    ThetaStar=NULL,
    knotAdded = FALSE,
    betaC = NULL,
    acceptance=list(knot.add=FALSE, knot.remove=FALSE, knot.move=FALSE,
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
  
  knots <- modelIteration.current$knots
  
  # If nknots>2, propose to move interior knot 
  if (length(knots) > 1) {
    
    # get boundaries and interior
    knots.boundary = range(knots)
    knots.interior = knots[-c(1,length(knots))] 
    
    #knots can only move to open bins between adjacent knots
    #Pick a knot to move: 
    knotToMove <- sample(knots, 1) 
    # get index of knot to be moved
    index <- which(knots == knotToMove) 
    # get the knots that are staying in the same place
    knotsToKeep <- sort(c(knots.interior, knots.boundary))[-index] 
    
    # find a new location from the potential knot locations
    potentialLocations <- knots.options$candidatePositions[!(knots.options$candidatePositions %in% knots)]
    # here we only allow movement within some small window
    potentialLocations <- potentialLocations[potentialLocations > (index - knots.options$stepSize[[group]] * 3) & 
                                               potentialLocations < (index + knots.options$stepSize[[group]] * 3)] 
    
    
    if (length(potentialLocations > 0)) {
      prop.move.1 <- prop.move.1 + 1
      if(length(potentialknots) == 1){
        p <- potentialknots
      }else{
        p <- sample(potentialknots,1)
      }
      
      pknotstar<-candidates.g1[!(candidates.g1 %in% c(keep,p))]
      pknotstar<-pknotstar[pknotstar>(p-step1*3) & pknotstar<(p+step1*3)]
      knotstar<-sort(c(keep, p))
      knotstar.b<-range(knotstar)
      knotstar.i<-knotstar[!(knotstar %in% knotstar.b)]
      #Calculate spline transformation of dropout time and proposed X matrix and residuals
      Xstar<-ns(u[group==1], knots=knotstar.i, Boundary.knots=knotstar.b, intercept=T)*t[group==1]
      ystar<-as.vector(y[group==1]-bt[i-1]-Xstar%*%ct.g1[[i]]-rep(B1[group.u==1],nobs[group.u==1])-t[group==1]*rep(B2[group.u==1],nobs[group.u==1])-C[group==1,]%*%covar[i-1,])
      yt<-as.vector(y[group==1]-bt[i-1]-Xt.g1%*%ct.g1[[i]]-rep(B1[group.u==1],nobs[group.u==1])-t[group==1]*rep(B2[group.u==1],nobs[group.u==1])-C[group==1,]%*%covar[i-1,])
      #Calculate the acceptance probability
      rho<-log(length(potentialknots))-log(length(pknotstar))+(crossprod(yt)-crossprod(ystar))/2/sig[i-1]
      if(rho>log(runif(1))){interior.g1[[i]]<-knotstar.i
                            boundary.g1[i,]<-knotstar.b
                            Xt.g1<-Xstar
                            accept.move.1<-accept.move.1+1
      }
    } else {
      # nowhere to move to
      return (knots, false, false)
    }
  } else if (length(knots) == 1) {
    
    # find a new location from the potential knot locations
    #     potentialknots<-candidates.g1[(candidates.g1 != boundary.g1[i,1])]
    potentialLocations <- knots.options$candidatePositions[knots.options$candidatePositions != knotToMove]
    # here we only allow movement within some small window
    potentialLocations <- potentialLocations[potentialLocations > (index - knots.options$stepSize[[group]]) & 
                                               potentialLocations < (index + knots.options$stepSize[[group]])] 
    # get the new location
    p <- sample(potentialLocations, 1)  
    
    potentialLocations.star <- knots.options$candidatePositions[(knots.options$candidatePositions != p)]
    potentialLocations.star <- potentialLocations.star[potentialLocations.star > (index - knots.options$stepSize[[group]]) & 
                                                         potentialLocations.star < (index + knots.options$stepSize[[group]])] 
    
    rho <- log(length(potentialLocations)) - log(length(potentialLocations.star))
    if (rho > log(runif(1))) {
      return (sort(knotsToKeep, p), true, true)
    } else {
      return (sort(knotsToKeep, p), true, false)
    }
  } else {
    return (NA, false, false)
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
updateFixedEffects <- function() {
  R0<-diag(rep(sig.beta, ((nknots.g1[i]+1)))) #Prior Variance
  yt<-y[group==1]-rep(B1[group.u==1],nobs[group.u==1])-t[group==1]*rep(B2[group.u==1],nobs[group.u==1])-C[group==1,]%*%covar[i-1,]   # Y - random effects
  mnt<-ginv(ginv(R0)+(1/sig[i-1])*crossprod(cbind(one[group==1],Xt.g1)))%*%(crossprod(cbind(one[group==1],Xt.g1),yt))*(1/sig[i-1]) #WLSP Mean
  covt<-ginv(ginv(R0)+(1/sig[i-1])*crossprod(cbind(one[group==1],Xt.g1))) #WLSP Cov
  covt<-vs*covt # Scale Cov to adjust acceptance rate
  covt<-as.matrix(nearPD(covt)$mat) #Ensure Cov is PD
  pbeta<-rmvnorm(1, mnt, covt) #Draw proposal
  #Calculate residuals for likelihood ratio
  ystar<-yt-tcrossprod(cbind(one[group==1],Xt.g1),pbeta)
  yt<-yt-cbind(one[group==1],Xt.g1)%*%c(bt[i-1], ct.g1[[i]])
  #Calculate acceptance prob
  rho<-log(dmvnorm(c(bt[i-1], ct.g1[[i]]), mnt,covt))-log(dmvnorm(pbeta, mnt,covt))+(crossprod(c(bt[i-1], ct.g1[[i]]))-tcrossprod(pbeta))/2/sig.beta+(crossprod(yt)-crossprod(ystar))/2/sig[i-1]
  
  if(rho>log(runif(1))){bt[i]<-pbeta[1]
                        ct.g1[[i]]<-pbeta[2:(nknots.g1[i]+1)]
                        accept.betas1<-accept.betas1+1
  }else{bt[i]<-bt[i-1]
  }
  
}


updateFixedEffectsCovariates <- function() {
  R0<-diag(rep(100, (ncovar))) #prior variance
  yt<-y-rep(B1,nobs)-t*rep(B2,nobs)
  yt[group==1]<-yt[group==1]-bt[i]-Xt.g1%*%as.vector(ct.g1[[i]])
  yt[group==2]<-yt[group==2]-bt[i]-hardeffect.g1[i]-Xt.g2%*%as.vector(ct.g2[[i]])
  mnt<-ginv(ginv(R0)+(1/sig[i-1])*t(cbind(C))%*%cbind(C))%*%(t(cbind(C))%*%yt)*(1/sig[i-1])
  covt<-ginv(ginv(R0)+(1/sig[i-1])*t(cbind(C))%*%cbind(C))
  covt<-vs*covt
  covt<-as.matrix(nearPD(covt)$mat)
  pbeta<-rmvnorm(1, mnt, covt)
  ystar<-yt-cbind(C)%*%t(pbeta)
  yt<-yt-cbind(C)%*%c(covar[i-1,])
  rho<-sum(log(dnorm(ystar, 0, sqrt(sig[i-1]))))+log(dmvnorm(pbeta, rep(0, (ncovar)),R0))+log(dmvnorm(c(covar[i-1,]), mnt,covt)) - sum(log(dnorm(yt, 0, sqrt(sig[i-1]))))-log(dmvnorm(c(covar[i-1,]), rep(0, (ncovar)),R0))-log(dmvnorm(pbeta, mnt,covt))
  if(rho>log(runif(1))){
    covar[i,]<-pbeta
    accept.covar<-accept.covar+1
  }else{
    covar[i,]<-covar[i-1,]}
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
updateRandomEffects <- function() {
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
updateCovarianceParameters <- function() {
  # Update tau  
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


