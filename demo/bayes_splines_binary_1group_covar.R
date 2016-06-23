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
data(sim_1)
data <- sim_1

# set the model options
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
                                             eta.null=0.01)

# set the columns to use in the model
ids.var = "patid"
outcomes.var = "yi_bin"
groups.var = "group"
covariates.var = "gender"
times.dropout.var = "drop"
times.observation.var = "t"

# set the model fitting method
method="bayes.splines"
# set the distribution of the outcome
dist = "binary"

# set a random seed
set.seed(1066)
# fit the model
fit = informativeDropout(data, ids.var, 
                         outcomes.var, groups.var,
                         covariates.var, 
                         times.dropout.var, times.observation.var,
                         method, dist, model.options)

# summarize the result
summary(fit)


