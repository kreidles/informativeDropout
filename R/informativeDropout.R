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
library(MASS)
library(Matrix)
library(nlme)
library(lme4)
library(mvtnorm)
library(MCMCpack)
library(gtools)
library(abind)

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




informativeDropout.bayes.dirichlet <- function(data, ids.var, outcomes.var, groups.var, 
                                               covariates.var, 
                                               times.dropout.var, times.observation.var, 
                                               dist, mcmc.options, prior.options, startValues) {
  
  numClusters <- startValues$numClusters
  clusterWeights <- rep(1/numClusters, numClusters)
  
  for (i in 2:mcmc.options$iterations) {
    for (group.index in groupList) {
      # Update the cluster indicator (c_i)
      weightProd[[group.index]] <- cumprod(1 - clusterWeights)
      pi.g0[1]<-v.g0[1]
      for (h in 2:H) pi.g0[h]<-v.g0[h]*cumv[h-1]
      for (h in 1:H) tmp2.g0[,h]<-pi.g0[h]*dmvnorm(betas.g0,mu.g0[h,],matrix(covar.g0, 3))
      p.g0<-tmp2.g0/apply(tmp2.g0,1,sum)
      
      C.g0[i,]<-c.g0<-rMultinom(p.g0,1)
      Pi.g0[i,]<-pi.g0
      for (h in 1:H) ns.g0[h]<-length(c.g0[c.g0==h])  # Must allow zeros for empty clusters
      
      
      # Update stick breaking weights (V_k)
      for (h in 1:(H-1)) v.g0[h]<-rbeta(1,1+ns.g0[h],alpha.g0+sum(ns.g0[(h+1):H]))
      V.g0[i,]<-v.g0
      
      ncluster.g0<-Ncluster.g0[i]<-length(ns.g0[ns.g0!=0])
      
      
      
      # sample from the stick breaking weights
      
      
      # sample the means for each cluster 
      
      
      # sample the covariance of the subject specific coefficients
      
      
      # Update the subject specific coefficients
      
      
      # Sample the hyperparameters for the baseline distribution
      
      
      # Sample the concentration parameter 
      
      
      # update the covariate effects
      
      
      # update the scale parameters if included in the model
      
      
      # calculate the marginal effects
    }

    

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
                                             covariates.var, 
                                             times.dropout.var, times.observation.var, 
                                             dist, knots.options, mcmc.options, prior.options,
                                             dropoutEstimationTimes) {
  
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
    knots.boundary = range(knots.options$startPositions)
    knots.interior = knots.options$startPositions[-c(1,length(knots.options$startPositions))] 
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
  
  # initialize the first model iteration
  # TODO: mcmc options specify starting values
  modelIterationList <- vector(mode = "list", length = mcmc.options$iterations)
  modelIterationList[[1]] = 
    rjmcmc.iteration(knots=lapply(1:length(groupList), function(i) { 
      return (knots.options$startPositions); }), 
      Theta=Theta.init, betaCovariate=betaCovariate.init,
      sigma.error = 1, sigma.spline = 1.25^2, sigma.randomIntercept = 1,
      sigma.randomSlope = 1, sigma.randomInterceptSlope = 0,
      shape.tau = 0.001, rate.tau = 0.001)
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
      if ((u < knots.options$birthProbability && 
           length(model.current$knots[[group.index]]) < knots.options$max) || 
          (length(model.current$knots[[group.index]]) <= knots.options$min)) {
        # add a knot
        result = addKnot(dist=dist,
                         knots.previous=model.current$knots[[group.index]], 
                         knots.options=knots.options, 
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
        # update the model iteration
        X[[group.index]] = result$X
        model.current$knots[[group.index]] = result$knots
        model.current$Theta[[group.index]] = result$Theta
        model.current$proposed$knot.add = TRUE
        model.current$accepted$knot.add = result$accepted
        
      } else {
        # remove a knot
        result = removeKnot(dist=dist, knots.previous=model.current$knots[[group.index]], 
                            knots.options=knots.options, 
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
      result = moveKnot(dist=dist, knots.previous=model.current$knots[[group.index]], 
                        knots.stepSize=knots.options$stepSize, 
                        knots.candidatePositions=knots.options$candidatePositions,
                        outcomes=group.outcomes, 
                        times.dropout=group.times.dropout, 
                        times.observation=group.times.observation, 
                        covariates=group.covariates,
                        X.previous=X[[group.index]], 
                        Theta.previous=model.current$Theta[[group.index]],
                        Z=group.Z, alpha=group.alpha, 
                        betaCovariates=model.current$betaCovariates,  
                        sigma.error=model.current$sigma.error)
      X[[group.index]] = result$X
      model.current$knots[[group.index]] = result$knots
      model.current$proposed$knot.move = result$proposed
      model.current$accepted$knot.move = result$accepted
      
      # update fixed effects (includes coefficients for covariates and time varying slopes)
      print("FIXED")
      result = updateFixedEffects(dist=dist,
                                  knots.previous=model.current$knots[[group.index]], 
                                  knots.options=knots.options, 
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
      model.current$Theta[[group.index]] = result$Theta
      model.current$proposed$fixedEffects = TRUE
      model.current$accepted$fixedEffects = result$accepted
      
    }  
    
    # update fixed effects associated with covariates
    if (!is.null(covariates)) {
      result = updateFixedEffectsCovariates(dist=dist, 
                                            outcomes=outcomes, 
                                            covariates=covariates, 
                                            X=X, 
                                            Theta=model.current$Theta, 
                                            Z=Z, alpha=alpha, 
                                            betaCovariates.previous=model.current$betaCovariates,  
                                            sigma.error=model.current$sigma.error, 
                                            sigma.beta=prior.options$sigma.beta)
      model.current$betaCovariates = result$betaCovariates
      model.current$proposed$fixedEffectsCovariates = TRUE
      model.current$accepted$fixedEffectsCovariates = result$accepted
    }
    
    # update random effects
    alpha = 
      updateRandomEffects(dist=dist, 
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
    
    # update variance components
    result = updateCovarianceParameters(dist=dist,
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

#' Fit a varying coefficient model for longitudinal studies with
#' informative dropout. 
#' 
#' @param 
#' @param 
#' @return 
#' @examples
#' 
informativeDropout <- function(data, ids.var, outcomes.var, groups.var, covariates.var, 
                               times.dropout.var, times.observation.var,
                               method="bayes.splines", dist="normal",
                               knots.options, 
                               mcmc.options=list(iterations=100000, burnIn=50000),
                               prior.options, dropoutEstimationTimes) {
  
  if (method == 'bayes.splines') {
    # model the relationship between dropout time and slope using natural splines
    return (informativeDropout.bayes.splines(data, ids.var, outcomes.var, groups.var, covariates.var, 
                                             times.dropout.var, times.observation.var, dist, 
                                             knots.options, mcmc.options, prior.options,
                                             dropoutEstimationTimes))
    
  } else if (method == 'bayes.dirichlet') {
    # account for informative dropout using a dirichlet process 
    return (informativeDropout.bayes.dirichlet(data, ids.var, outcomes.var, groups.var, covariates.var, 
                                               times.dropout.var, times.observation.var, dist, prior.options))
  } else if (method == 'mixed') {
    # fit a mixed model which models the relationship between dropout time and slope using natural splines
    return (informativeDropout.mixed(data, ids.var, outcomes.var, groups.var, covariates.var, 
                                     times.dropout.var, times.observation.var, dist, dist))
  }
  
}

acceptanceProbability <- function(fit, action) {
  if (action == "knot.add") {
    total_accepts = sum(sapply(fit, function(x) { return(as.numeric(x$accepted$knot.add)) }))
  } else if (action == "knot.remove") {
    total_accepts = sum(sapply(fit, function(x) { return(as.numeric(x$accepted$knot.remove)) }))
  } else if (action == "knot.move") {
    total_accepts = sum(sapply(fit, function(x) { return(as.numeric(x$accepted$knot.move)) }))
    
  } else if (action == "fixedEffects") {
    total_accepts = sum(sapply(fit, function(x) { return(as.numeric(x$accepted$fixedEffects)) }))
  } else if (action == "fixedEffectsCovariates") {
    total_accepts = sum(sapply(fit, function(x) { return(as.numeric(x$accepted$fixedEffectsCovariates)) }))
  }
  
  
  return(total_accepts/length(fit))
}

test.example <- function() {
  data <- read.csv("test.csv")
  data$day = data$years * 365
  
  # for debugging
  ids.var = "WIHSID"
  outcomes.var = "logcd4"
  groups.var = "hard"
  covariates.var = c("AGEATBL", "minority")
  times.dropout.var = "drop"
  times.observation.var = "day"
  method="bayes.splines"
  dist = "gaussian"
  knots.options=list(birthProbability=0.5, min=3, max=10, stepSize=3,
                     startPositions=c(330,550,1060), candidatePositions=seq(10,max(data$day),10)) 
  mcmc.options=list(iterations=20, burnIn=10, sigma.residual=1)
  prior.options=list(shape.tau = 0.001, rate.tau = 0.001, lambda.numKnots = 1,
                     sigma.beta = 1, sigmaError.df = 3, sigmaError.scaleMatrix = diag(2))
  
  set.seed(1066)
  result = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
                              times.dropout.var, times.observation.var, 
                              method, dist,
                              knots.options = knots.options, 
                              mcmc.options = mcmc.options,
                              prior.options = prior.options)
}

test.sim <- function() {
  data <- read.table("../../Rnsv/code/sim_sml_1.dat")
  #data$day = data$years * 365
  
  names(data) <- c("patid", "alpha", "drptm", "b1", "b2",
                   "b2ui", "b2uii", "b2uiii", "t", "e", "yi", "yii", "yiii")
  data$group <- rep(1,nrow(data))
  # for debugging
  ids.var = "patid"
  outcomes.var = "yiii"
  groups.var = "group"
  covariates.var = NULL
  times.dropout.var = "drptm"
  times.observation.var = "t"
  method="bayes.splines"
  dist = "gaussian"
  knots.options=list(birthProbability=0.2, min=1, max=10, stepSize=0.1,
                     startPositions=c(0,7/30,0.5, 23/30,1), candidatePositions=seq(0,1,0.1/3)) 
  mcmc.options=list(iterations=40000, burnIn=10, sigma.residual=1.25^2)
  prior.options=list(shape.tau = 0.001, rate.tau = 0.001, lambda.numKnots = 5,
                     sigma.beta = 25, sigmaError.df = 3, 
                     sigmaError.scaleMatrix = diag(2))
  start.options=list()
  dropoutEstimationTimes = seq(0,1,1/15)
  
  set.seed(1066)
  result = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
                              times.dropout.var, times.observation.var, 
                              method, dist,
                              knots.options = knots.options, 
                              mcmc.options = mcmc.options,
                              prior.options = prior.options,
                              dropoutEstimationTimes = dropoutEstimationTimes)
  
  
  acceptanceProbability(result, "knot.add")
  acceptanceProbability(result, "knot.remove")
  acceptanceProbability(result, "knot.move")
  acceptanceProbability(result, "fixedEffects")
  #acceptanceProbability(result, "fixedEffectsCovariates")
  
  
  nknots = unlist(lapply(result, function(x) { return(length(x$knots[[1]])) } ))
  ts.plot(nknots)
  summary(nknots)
  
  slopes = unlist(lapply(result, function(x) { return(x$slope.marginal[[1]]) }))
  summary(slopes)
  ts.plot(slopes)
  
  sum.sigma.error = unlist(lapply(result, function(x) { return(x$sigma.error) }))
  summary(sum.sigma.error)
  ts.plot(sum.sigma.error)
  
  sum.sigma.randomIntercept = unlist(lapply(result, function(x) { return(x$sigma.randomIntercept) }))
  summary(sum.sigma.randomIntercept)
  ts.plot(sum.sigma.randomIntercept)
  
  sum.sigma.randomSlope = unlist(lapply(result, function(x) { return(x$sigma.randomSlope) }))
  summary(sum.sigma.randomSlope)
  ts.plot(sum.sigma.randomSlope)
  
  sum.sigma.randomInterceptSlope = unlist(lapply(result, function(x) { return(x$sigma.randomInterceptSlope) }))
  summary(sum.sigma.randomInterceptSlope)
  ts.plot(sum.sigma.randomInterceptSlope)
  
  
  dropout.slopes1 = unlist(lapply(result, function(x) { return(x$slope.dropoutSpecific[[1]][1])}))
  summary(dropout.slopes1)
  ts.plot(dropout.slopes1)
  
  dropout.slopes2 = unlist(lapply(result, function(x) { return(x$slope.dropoutSpecific[[1]][2])}))
  summary(dropout.slopes2)
  ts.plot(dropout.slopes2)
  
  dropout.slopes3 = unlist(lapply(result, function(x) { return(x$slope.dropoutSpecific[[1]][3])}))
  summary(dropout.slopes3)
  ts.plot(dropout.slopes3)
}






