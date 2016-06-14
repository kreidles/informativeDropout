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
library(Hmisc)

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
  result <- solve(as.matrix(nearPD(crossprod(X, weight %*% X))$mat)) %*% (crossprod(X, weight %*% y.wls))
  return(as.vector(result))
} 


#
# Register some S3 generic methods for summarizing the fit object
#

# trace plot for a parameter
plot.trace <- function (x, ...) {
  UseMethod("plot.trace", x)
}
# density plot for a parameter
plot.density <- function (x, ...) {
  UseMethod("plot.density", x)
}
# plot the slope by dropout time
plot.slopeByDropout <- function (x, ...) {
  UseMethod("plot.slopeByDropout", x)
}
# perform sensitivity analysis on the slope results
sensitivity <- function(x, ...) {
  # delta factor is multiplier on slope after dropout
  # y = int + drop * slope + (time - dropout) * delta * slope
  # arguments: min/max times, list of estimation times, 
  # multiple deltas
  # vector of covariate values by time
  UseMethod("sensitivity.slope")
}


#' Fit a varying coefficient model for longitudinal studies with
#' informative dropout. 
#' 
#' @param 
#' @param 
#' @return 
#' @examples
#' 
#' @export informativeDropout
#' 
informativeDropout <- function(data, ids.var, outcomes.var, groups.var, covariates.var, 
                               times.dropout.var, times.observation.var,
                               method="bayes.splines", dist="normal",
                               model.options) {
  
  if (method == 'bayes.splines') {
    # model the relationship between dropout time and slope using natural splines
    return (informativeDropout.bayes.splines(data, ids.var, outcomes.var, groups.var, covariates.var, 
                                             times.dropout.var, times.observation.var, dist, model.options))
    
  } else if (method == 'dirichlet') {
    # account for informative dropout using a dirichlet process 
    return (informativeDropout.bayes.dirichlet(data, ids.var, outcomes.var, groups.var, covariates.var, 
                                               times.dropout.var, times.observation.var, dist, model.options))
  } else if (method == 'mixed.splines') {
    # fit a mixed model which models the relationship between dropout time and slope using natural splines
    return (informativeDropout.mixed(data, ids.var, outcomes.var, groups.var, covariates.var, 
                                     times.dropout.var, times.observation.var, dist, 
                                     model.options))
  }
  
}



