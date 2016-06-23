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

#' @include informativeDropout.R

#####################################
# Bayesian spline model examples
#####################################


#'
#'
example.bayes_splines_gaussian_2group_covar <- function() {
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
  
  
#   sigma.beta=NULL, sigma.residual=NULL,
#   sigma.error=0,
#   sigma.error.shape.tau=NULL, sigma.error.rate.tau=NULL,
#   lambda.numKnots=NULL,
#   sigma.randomIntercept = NULL,
#   sigma.randomSlope = NULL,
#   sigma.randomInterceptSlope = NULL,
#   sigma.randomEffects.df = NULL,
#   sigma.randomEffects.scale = NULL,
#   eta.null=NULL
  
  model.options <- bayes.splines.model.options(iterations=100, burnin=10, thin=1, print=1,
                                               knots.prob.birth=0.5, knots.min=3, knots.max=10, knots.stepSize=3,
                                               knots.positions.start=c(330,550,1060), 
                                               knots.positions.candidate=seq(10,max(data$day),10),
                                               dropout.estimationTimes=c(1,2,3),
                                               sigma.error.shape.tau=0.001, sigma.error.rate.tau=0.001,
                                               sigma.beta=1, lambda.numKnots=1,
                                               sigma.residual=1,
                                               sigma.randomEffects.df = 3,
                                               sigma.randomEffects.scale = diag(2))
  
  
  set.seed(1066)
  fit = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
                              times.dropout.var, times.observation.var, 
                              method, dist, model.options)  
}


example.bayes_splines_gaussian_1group_nocovar <- function() {
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
  
  model.options <- bayes.splines.model.options(iterations=100, burnin=10, thin=1,
                                               knots.prob.birth=0.2, knots.min=1, knots.max=10, 
                                               knots.stepSize=0.1,
                                               knots.positions.start=c(0,7/30,0.5, 23/30,1), 
                                               knots.positions.candidate=seq(0,1,0.1/3),
                                               dropout.estimationTimes=seq(0,1,1/15),
                                               sigma.error.shape.tau=0.001, sigma.error.rate.tau=0.001,
                                               sigma.beta=25, lambda.numKnots=5,
                                               sigma.residual = 1,
                                               sigma.error=1,
                                               sigma.randomIntercept = 1,
                                               sigma.randomSlope = 1,
                                               sigma.randomInterceptSlope = 0.001,
                                               sigma.randomEffects.df = 3,
                                               sigma.randomEffects.scale = diag(2),
                                               eta.null=NULL)
  model.options$sigma.randomIntercept
  
  
  set.seed(1066) 
  fit = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
                           times.dropout.var, times.observation.var, 
                           method, dist, model.options)
  
  summary(fit)
  
  
  
  nknots = unlist(lapply(fit$iterations, function(x) { return(length(x$knots[[1]])) } ))
  ts.plot(nknots)
  summary(nknots)
  
  slopes = unlist(lapply(fit$iterations, function(x) { return(x$slope.marginal[[1]]) }))
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





example.bayes_splines_binary_1group_nocovar <- function() {
  data <- read.table("../../Rnsv/code/sim_sml_1.dat")
  #data$day = data$years * 365
  
  names(data) <- c("patid", "alpha", "drptm", "b1", "b2",
                   "b2ui", "b2uii", "b2uiii", "t", "e", "yi", "yii", "yiii")
  data$group <- rep(1,nrow(data))
  data$yi_bin = (data$yi> 0)
  data$drop <- data$drptm + 1/15
  
  # for debugging
  ids.var = "patid"
  outcomes.var = "yi_bin"
  groups.var = "group"
  covariates.var = NULL
  times.dropout.var = "drop"
  times.observation.var = "t"
  method="bayes.splines"
  dist = "binary"
  
  model.options=bayes.splines.model.options(iterations=100, n.clusters=15, burnin=0,
                                        dropout.estimationTimes = seq(1/15,1,1/15),
                                        dp.concentration=1,
                                        dp.concentration.alpha=1,
                                        dp.concentration.beta=1,
                                        dp.cluster.sigma = diag(3),
                                        dp.cluster.sigma.nu0 = 5,
                                        dp.cluster.sigma.T0 = diag(3),
                                        dp.dist.mu0 = c(0,0,0),
                                        dp.dist.mu0.mb = c(0,0,0),
                                        dp.dist.mu0.Sb = diag(3),
                                        dp.dist.sigma0 = diag(3),
                                        dp.dist.sigma0.nub = 5,
                                        dp.dist.sigma0.Tb = diag(3),
                                        betas.covariates = NULL,
                                        betas.covariates.mu = NULL,
                                        betas.covariates.sigma = NULL,
                                        sigma.error.tau=0.01)
  
  
  set.seed(1066)
  fit = informativeDropout(data, ids.var, 
                           outcomes.var, groups.var,
                           covariates.var, 
                           times.dropout.var, times.observation.var,
                           method, dist, model.options)
  
  summary(fit)
  
}

example.bayes_splines_binary_1group_covar <- function() {
  data <- read.table("../../Rnsv/code/sim_sml_1.dat")
  #data$day = data$years * 365
  
  names(data) <- c("patid", "alpha", "drptm", "b1", "b2",
                   "b2ui", "b2uii", "b2uiii", "t", "e", "yi", "yii", "yiii")
  data$group <- rep(1,nrow(data))
  data$yi_bin = (data$yi> 0)
  data$gender = data$yi_bin
  data$drop <- data$drptm + 1/15
  
  data$gender[data$yi_bin == 1] <- as.numeric(runif(length(data$yi_bin[data$yi_bin == 1])) > 0.8)
  data$gender[data$yi_bin == 0] <- as.numeric(runif(length(data$yi_bin[data$yi_bin == 0])) > 0.15)
  
  # for debugging
  ids.var = "patid"
  outcomes.var = "yi_bin"
  groups.var = "group"
  covariates.var = "gender"
  times.dropout.var = "drop"
  times.observation.var = "t"
  method="bayes.splines"
  dist = "binary"
  
  model.options=bayes.splines.model.options(iterations=100, n.clusters=15, burnin=0,
                                        dropout.estimationTimes = seq(1/15,1,1/15),
                                        dp.concentration=1,
                                        dp.concentration.alpha=1,
                                        dp.concentration.beta=1,
                                        dp.cluster.sigma = diag(3),
                                        dp.cluster.sigma.nu0 = 5,
                                        dp.cluster.sigma.T0 = diag(3),
                                        dp.dist.mu0 = c(0,0,0),
                                        dp.dist.mu0.mb = c(0,0,0),
                                        dp.dist.mu0.Sb = diag(3),
                                        dp.dist.sigma0 = diag(3),
                                        dp.dist.sigma0.nub = 5,
                                        dp.dist.sigma0.Tb = diag(3),
                                        betas.covariates = NULL,
                                        betas.covariates.mu = 0,
                                        betas.covariates.sigma = matrix(0.7),
                                        sigma.error.tau = 0.01)
  
  
  set.seed(1066)
  fit = informativeDropout(data, ids.var, 
                           outcomes.var, groups.var,
                           covariates.var, 
                           times.dropout.var, times.observation.var,
                           method, dist, model.options)
  
  summary(fit)
}

