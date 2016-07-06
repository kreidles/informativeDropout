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
  # delta factor is multiplier on slope after dropout
  # y = int + drop * slope + (time - dropout) * delta * slope
  # arguments: min/max times, list of estimation times, 
  # multiple deltas
  # vector of covariate values by time
  UseMethod("sensitivity.slope", fit)
}

#'
#' Single iteration of weighted least squares
#' 
#' 
wls.binary <- function(y, X, eta.wls, model.options) { 
  
  # get the probability from the linear predictor
  prob <- inv.logit(eta.wls)
  # truncate at min/max values for computational stability
  prob[prob < model.options$prob.min] <- model.options$prob.min
  prob[prob > model.options$prob.max] <- model.options$prob.max
  
  # calculate the weight matrix and estimates 
  var = (prob * (1 - prob))
  y.wls <- model.options$eta.null + (y - prob) * (1 / var)
  weight <- Diagonal(x = var)
  result <- solve(as.matrix(nearPD(crossprod(X, as.matrix(weight %*% X)))$mat)) %*% (crossprod(X, as.matrix(weight %*% y.wls)))
  return(as.vector(result))
} 

#' Fit a simple linear model to get initial values for regression
#' coefficients associated with covariates
#' 
#' @param dist the distribution of the outcome ("gaussian" or "binary") 
#' @param covariates date frame containing covariate values
#' @param outcomes vector of outcomes
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

#' Fit a simple linear model to get initial values for regression
#' coefficients associated with splines
#' 
#' @param dist the distribution of the outcome ("gaussian" or "binary") 
#' @param covariates date frame containing covariate values
#' @param outcomes vector of outcomes
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
    # binomial
    fit.Theta <- glm(formula, family=binomial, data=data.theta)
    return (lapply(1:length(groupList), function(i) {
      return (as.vector(coef(fit.Theta)))
    }))
  }
}

