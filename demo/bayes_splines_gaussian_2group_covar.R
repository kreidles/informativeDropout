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

# load the data set
data(hiv)
data <- hiv
#data <- data[data$hard==0,]
# set the model options
model.options <- bayes.splines.model.options(
  iterations=500, burnin=0, thin=1, print=100,
  knots.prob.birth=0.2, knots.min=1, knots.max=10, knots.stepSize=30,
  knots.positions.start=list(c(142, 562, 982, 1402, 1822), c(139, 559, 979, 1399, 1819)),
  knots.positions.candidate=list(seq(142,1826,15),seq(139,1819,15)),
  dropout.estimationTimes=c(1,2,3),
  sigma.error=1,
  sigma.error.shape.tau=0.001, sigma.error.rate.tau=0.001,
  sigma.beta=25, lambda.numKnots=5,
  sigma.residual=1.25^2,
  sigma.randomIntercept = 1,
  sigma.randomSlope = 1,
  sigma.randomInterceptSlope = 0,
  sigma.randomEffects.df = 3,
  sigma.randomEffects.scale = diag(2))

# select columns to include in the model
ids.var = "WIHSID"
outcomes.var = "logcd4"
groups.var = "hard"
#groups.var=NULL
data$timebaselogcd4 = data$baselogcd4 * data$years
covariates.var = c("baselogcd4", "timebaselogcd4")
#covariates.var = NULL
times.dropout.var = "drop"
times.observation.var = "years"

# set the model fitting method
method="bayes.splines"
# select the distribution of the outcome
dist = "gaussian"

# set a random seed
set.seed(1066)

# fit the model
fit = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
                         times.dropout.var, times.observation.var, 
                         method, dist, model.options) 

# summarize the result
summary(fit)

