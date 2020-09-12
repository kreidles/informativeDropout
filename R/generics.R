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

#
# Register some S3 generic methods for summarizing the fit object
#

#' trace plot for a parameter
#'
#' @param fit the model fit object
#' @export plot.trace
#' 
plot.trace <- function (fit, ...) {
  UseMethod("plot.trace", fit)
}

#' density plot for a parameter
#'  
#' @param fit the model fit object
#' @export
#' 
plot.density <- function (fit, ...) {
  UseMethod("plot.density", fit)
}

#' plot the slope by dropout time
#'
#' @param fit the model fit object
#' @export  
plot.slopeByDropout <- function (fit, ...) {
  UseMethod("plot.slopeByDropout", fit)
}

#' perform sensitivity analysis on the slope results
#'
#' @param fit the model fit object
#' @export 
sensitivity <- function(fit, ...) {
  UseMethod("sensitivity", fit)
}

#'
#' Single iteration of weighted least squares
#' 
#' @param y vector of outcomes
#' @param X data frame of predictors
#' @param eta.wls current value of the linear predictor
#' @param model.options the model options for the current run
#' @param eta.null the value of the linear predictor for a null model - logit(p event)
#' 
wls.binary <- function(y, X, eta.wls, model.options, eta.null) { 
  
  # get the probability from the linear predictor
  prob <- inv.logit(eta.wls)
  # truncate at min/max values for computational stability
  prob[prob < model.options$prob.min] <- model.options$prob.min
  prob[prob > model.options$prob.max] <- model.options$prob.max
  
  # calculate the weight matrix and estimates 
  var = (prob * (1 - prob))
  y.wls <- eta.null + (y - prob) * (1 / var)
  weight <- Diagonal(x = var)
  result <- solve(as.matrix(nearPD(crossprod(X, as.matrix(weight %*% X)))$mat)) %*% (crossprod(X, as.matrix(weight %*% y.wls)))
  return(as.vector(result))
} 

#' Fit a simple linear model to get initial values for regression
#' coefficients associated with covariates
#' 
#' @param dist the distribution of the outcome ("gaussian" or "binary") 
#' @param covariates data frame containing covariate values
#' @param times vector of observation times
#' @param outcomes vector of outcomes
#' 
#'
getInitialEstimatesCovariates <- function(dist, covariates, times, outcomes) {
  if (is.null(covariates) || ncol(covariates) == 0) {
    return (NULL)
  }
  
  data.covar = cbind(outcomes, times,covariates)
  formula = as.formula(paste(c(paste(names(outcomes), "~"), 
                               paste(names(data.covar[,-1]), collapse=" + ")), collapse=" "))
  if (dist == 'gaussian') {
    fit.beta <- lm(formula, data=data.covar)
    return (as.vector(coef(fit.beta))[-c(1,2)])
  } else if (dist == 'binary') {
    fit.beta <- glm(formula, family=binomial, data=data.covar)
    return (as.vector(coef(fit.beta))[-c(1,2)])
  } else {
    stop("unsupported distribution")
  }
}

#' Fit a simple linear model to get initial values for regression
#' coefficients associated with splines
#' 
#' @param dist the distribution of the outcome ("gaussian" or "binary") 
#' @param groupList list of groups
#' @param groups vector of group assignments
#' @param X data frame containing predictor variables (assumes a column of ones for an intercept)
#' @param covariates data frame containing covariate values
#' @param outcomes vector of outcomes
#' 
getInitialEstimatesTheta <- function(dist, groupList, groups, X, covariates, outcomes) {
  data.theta = cbind(outcomes, X)
  if(!is.null(covariates)){data.theta = cbind(data.theta, covariates)}
  formula = as.formula(paste(c(paste(names(outcomes), "~"), 
                               paste(names(X)[2:length(names(X))], 
                                     collapse=" + ")), 
                             collapse=" "))
  if (dist == 'gaussian') {
    return (lapply(1:length(groupList), function(i) {
      group = groupList[i]
      formula = as.formula(paste(c(paste(names(outcomes), "~"), 
                                   paste(c(names(X)[2:length(names(X))],names(covariates)), 
                                         collapse=" + ")), 
                                 collapse=" "))
      fit.Theta <- lm(formula, data=data.theta[groups==group,])
      
      return (as.vector(coef(fit.Theta)[1:length(names(X))]))
    }))
  } else {
    # binomial
    return (lapply(1:length(groupList), function(i) {
      group = groupList[i]
      formula = as.formula(paste(c(paste(names(outcomes), "~"), 
                                   paste(c(names(X)[2:length(names(X))],names(covariates)), 
                                         collapse=" + ")), 
                                 collapse=" "))
      fit.Theta <- glm(formula, family='binomial',data=data.theta[groups==group,])
      
      return (as.vector(coef(fit.Theta)[1:length(names(X))])) 
    }))
  }
}

