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
data$group <- rep(1,nrow(data))

# set the model options
model.options <- bayes.splines.model.options(iterations=100, burnin=10, thin=1,
                                             knots.prob.birth=0.2, knots.min=1, knots.max=10, 
                                             knots.stepSize=0.1,
                                             knots.positions.start=c(0,7/30,0.5, 23/30,1), 
                                             knots.positions.candidate=seq(0,1,0.1/3),
                                             dropout.estimationTimes=seq(0,1,1/15),
                                             sigma.error.shape.tau=0.001, sigma.error.rate.tau=0.001,
                                             sigma.beta=25, lambda.numKnots=5,
                                             sigma.residual = 1.25^2,
                                             sigma.error=1,
                                             sigma.randomIntercept = 1,
                                             sigma.randomSlope = 1,
                                             sigma.randomInterceptSlope = 0.001,
                                             sigma.randomEffects.df = 3,
                                             sigma.randomEffects.scale = diag(2),
                                             eta.null=NULL)

# set the columns to use in the model
ids.var = "patid"
outcomes.var = "yiii"
groups.var = "group"
covariates.var = NULL
times.dropout.var = "drptm"
times.observation.var = "t"

# set the model fitting method
method="bayes.splines"
# set the distribution of the outcome
dist = "gaussian"

# set a random seed
set.seed(1066) 

# fit the model
fit = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
                         times.dropout.var, times.observation.var, 
                         method, dist, model.options)
# summarize the result
summary(fit)



# nknots = unlist(lapply(fit$iterations, function(x) { return(length(x$knots[[1]])) } ))
# ts.plot(nknots)
# summary(nknots)
# 
# slopes = unlist(lapply(fit$iterations, function(x) { return(x$slope.marginal[[1]]) }))
# summary(slopes)
# ts.plot(slopes)
# 
# sum.sigma.error = unlist(lapply(result, function(x) { return(x$sigma.error) }))
# summary(sum.sigma.error)
# ts.plot(sum.sigma.error)
# 
# sum.sigma.randomIntercept = unlist(lapply(result, function(x) { return(x$sigma.randomIntercept) }))
# summary(sum.sigma.randomIntercept)
# ts.plot(sum.sigma.randomIntercept)
# 
# sum.sigma.randomSlope = unlist(lapply(result, function(x) { return(x$sigma.randomSlope) }))
# summary(sum.sigma.randomSlope)
# ts.plot(sum.sigma.randomSlope)
# 
# sum.sigma.randomInterceptSlope = unlist(lapply(result, function(x) { return(x$sigma.randomInterceptSlope) }))
# summary(sum.sigma.randomInterceptSlope)
# ts.plot(sum.sigma.randomInterceptSlope)
# 
# 
# dropout.slopes1 = unlist(lapply(result, function(x) { return(x$slope.dropoutSpecific[[1]][1])}))
# summary(dropout.slopes1)
# ts.plot(dropout.slopes1)
# 
# dropout.slopes2 = unlist(lapply(result, function(x) { return(x$slope.dropoutSpecific[[1]][2])}))
# summary(dropout.slopes2)
# ts.plot(dropout.slopes2)
# 
# dropout.slopes3 = unlist(lapply(result, function(x) { return(x$slope.dropoutSpecific[[1]][3])}))
# summary(dropout.slopes3)
# ts.plot(dropout.slopes3)
