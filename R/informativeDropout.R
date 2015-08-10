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
modelInformation <- function(knots=NULL, method="bayesian", X=NULL, ) {
  mi = list(
    knots = knots,
    X=X,
    XStar=X,
    beta0=NULL,
    ThetaStar=NULL,
    knotAdded = FALSE,
    acceptance=list(knot.add=FALSE, knot.remove=FALSE, knot.move=FALSE,
                    fixedEffects=FALSE, randomEffects=FALSE, varianceComponents=FALSE)
    
  )
  
  class(mi) <- append(class(mi), "modelInformation")
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
addKnot <- function(data, modelInformation.previous, knots.options) {
  # add a knot by randomly selecting a candidate knot
  index = sample(1:length(knots.options$candidatePositions), 1)
  newKnot.value = knots.candidatePositions[index]
  knots <- sort(c(knotPositions, newKnot.value))
  # get the interior and boundary knots, and grab the position of the knot that
  # was just added
  knots.boundary = range(knots)
  knots.interior = knots[-c(1,length(knots))] 
  newKnot.position = which(knotstar == newKnot.value)
  
  # Calculate spline transformation of dropout time and create the proposed X matrix
  Xstar <- ns(data$dropoutTimes, knots=knots.interior, Boundary.knots=knots.boundary, intercept=T)*data$times

  #Calc y-random effects for least squares calculations
  yls<-as.vector(y[group==1]-bt[i-1]-rep(B1[group.u==1],nobs[group.u==1])-t[group==1]*rep(B2[group.u==1],nobs[group.u==1])-C[group==1,]%*%covar[i-1,])
  #Calculate least squares estimates for coefficients and differences between LS and current coefficients
  LS<-ginv(crossprod(Xt.g1))%*%(crossprod(Xt.g1,yls))
  LSstar<-ginv(crossprod(Xstar))%*%(crossprod(Xstar,yls))
  LSresid<-ct.g1[[i-1]]-LS
  #Draw a residual for the added coefficient and calculate coefficient transformation
  cjresid<-rnorm(1, 0, sd.resid)
  cj<-LSstar[knotj]+cjresid
  coeff<-LSstar[-knotj]+LSresid
  if (knotj==1){cstar<-c(cj, coeff)}
  if (knotj==nstar){cstar<-c(coeff, cj)}
  if (knotj<nstar & knotj>1){cstar<-c(coeff[1:(knotj-1)], cj, coeff[knotj:nknots.g1[i-1]])}
  #Calculate residuals for likelihood ratio
  ystar<-as.vector(y[group==1]-bt[i-1]-Xstar%*%cstar-rep(B1[group.u==1],nobs[group.u==1])-t[group==1]*rep(B2[group.u==1],nobs[group.u==1])-C[group==1,]%*%covar[i-1,])
  yt<-as.vector(y[group==1]-bt[i-1]-Xt.g1%*%ct.g1[[i-1]]-rep(B1[group.u==1],nobs[group.u==1])-t[group==1]*rep(B2[group.u==1],nobs[group.u==1])-C[group==1,]%*%covar[i-1,])
  #Calculate birth and death probabilities                                                        
  bk<-ifelse(nknots.g1[i-1]==kmin, 1, pbirth)
  dk<-ifelse(nknots.g1[i-1]==(kmax-1),1,1-pbirth)
  #Calculate Acceptance Probability                                                        
  rho<-log(poi.lam)-log(nstar-1)+log(dk)-log(bk)+log(sd.resid)-0.5*log(sig.beta)+(cjresid^2/2/sd.resid/sd.resid)+((crossprod(ct.g1[[i-1]])-crossprod(cstar))/sig.beta/2)+(crossprod(yt)-crossprod(ystar))/2/sig[i-1]
  if(rho>log(runif(1))){nknots.g1[i]<-nstar
                        interior.g1[[i]]<-knotstar.i
                        boundary.g1[i,]<- knotstar.b
                        Xt.g1<-Xstar
                        accept.bd1<-accept.bd1+1
                        ct.g1[[i]]<-cstar
  }else{interior.g1[[i]]<-interior.g1[[i-1]]
        boundary.g1[i,]<-boundary.g1[i-1,]
        ct.g1[[i]]<-ct.g1[[i-1]]
        nknots.g1[i]<-nknots.g1[i-1]}                      
  
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
removeKnot <- function() {
  # randomly remove an existing knot
  index = sample(1:length(knotPositions), 1)
  knotPositions <- knotPositions[-index]
  
  knotsd<-sort(c(interior.g1[[i-1]],boundary.g1[i-1,]))
  d<-sample(knotsd, 1)
  
  nstar<-nknots.g1[i-1]-1
  knotstar<-sort(subset(knotsd, knotsd != d))
  
  if (nstar>1){
    knotstar.b<-range(knotstar)
    
    knotstar.i<-knotstar[!(knotstar %in% knotstar.b)]
    
    Xstar<-ns(u[group==1], knots=knotstar.i, Boundary.knots=knotstar.b, intercept=T)*t[group==1]}
  
  if (nstar==1){knotstar.b<-c(knotstar, NA)
                
                knotstar.i<-NA
                
                Xstar<-cbind(t[group==1])}
  
  
  
  
  knotj<-which(knotsd==d)
  
  
  yls<-as.vector(y[group==1]-bt[i-1]-rep(B1[group.u==1],nobs[group.u==1])-t[group==1]*rep(B2[group.u==1],nobs[group.u==1])-C[group==1,]%*%covar[i-1,])
  
  LS<-ginv(crossprod(Xt.g1))%*%(crossprod(Xt.g1,yls))
  LSstar<-ginv(crossprod(Xstar))%*%(crossprod(Xstar,yls))
  LSresid<-ct.g1[[i-1]]-LS
  cjresid<-LSresid[knotj]
  LSresid<-LSresid[-knotj]
  
  
  cstar<-LSstar+LSresid
  
  ystar<-as.vector(y[group==1]-bt[i-1]-Xstar%*%cstar-rep(B1[group.u==1],nobs[group.u==1])-t[group==1]*rep(B2[group.u==1],nobs[group.u==1])-C[group==1,]%*%covar[i-1,])
  yt<-as.vector(y[group==1]-bt[i-1]-Xt.g1%*%ct.g1[[i-1]]-rep(B1[group.u==1],nobs[group.u==1])-t[group==1]*rep(B2[group.u==1],nobs[group.u==1])-C[group==1,]%*%covar[i-1,])
  bk<-ifelse(nstar==kmin, 1, pbirth)
  dk<-ifelse(nknots.g1[i-1]==kmax,1,1-pbirth)
  rho<--log(poi.lam)+log(nknots.g1[i-1]-1)-log(dk)+log(bk)-log(sd.resid)+0.5*log(sig.beta)-(cjresid^2/2/sd.resid/sd.resid)+((crossprod(ct.g1[[i-1]])-crossprod(cstar))/sig.beta/2)+(crossprod(yt)-crossprod(ystar))/2/sig[i-1]
  
  if(rho>log(runif(1))){nknots.g1[i]<-nstar
                        interior.g1[[i]]<-knotstar.i
                        boundary.g1[i,]<-knotstar.b
                        Xt.g1<-Xstar
                        accept.bd1<-accept.bd1+1
                        ct.g1[[i]]<-cstar
  }else{interior.g1[[i]]<-interior.g1[[i-1]]
        boundary.g1[i,]<-boundary.g1[i-1,]
        ct.g1[[i]]<-ct.g1[[i-1]]
        nknots.g1[i]<-nknots.g1[i-1]}           
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
moveKnot <- function() {
  
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
_informativeDropout.bayesian.gaussian <- function(iterations, data, covariates, dropoutTimes, 
                                                  times, knots.options, mcmc.options) {
  

  modelInformation = list(iterations)
  
  knotPositions <- knots.startPositions
  
  for (i in 2:iterations) {
    for (group in groupList) {
      # randomly decide to add/remove a knot
      if (runif(1) < knot.birthProbability && ) {
        modelInformation[i] = addKnot(modelInformation[i-1], knots.options)
      } else {
        modelInformation[i] = removeKnot(modelInformation[i-1], knots.options)                 
      }
      
      # Move knots
      modelInformation[i] = moveKnots(modelInformation[i], modelInformation[i-1], knots.options)
      
      # update fixed effects (includes coefficients for covariates and time varying slopes)
      modelInformation[i] = updateFixedEffects(modelInformation[i], modelInformation[i-1], knots.options)
      
      # update random effects
      modelInformation[i] = updateRandomEffects(modelInformation[i], modelInformation[i-1], knots.options)
      
      # update variance components
      modelInformation[i] = updateVarianceComponents(modelInformation[i], modelInformation[i-1], knots.options)
    }    
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


