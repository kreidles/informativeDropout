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

#####################################
# Demo of a dirichlet process model with a gaussian outcome
#####################################

# load the simulated data
data(sim_1)
data <- sim_1

# define the model options
model.options=dirichlet.model.options(iterations=100, n.clusters=15, burnin=0, print=10,
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
                                      sigma.error = 1,
                                      sigma.error.tau = 0.01,
                                      density.intercept.domain = seq(-4,4,0.01),
                                      density.slope.domain = seq(-4,4,0.01))

# Set the columns to be used in the model
ids.var = "patid"
outcomes.var = "yi"
groups.var = NULL 
covariates.var = NULL
times.dropout.var = "drop"
times.observation.var = "t"
method="dirichlet"
dist = "gaussian"

# set a random seed if desired
set.seed(1066)

# fit the model
fit = informativeDropout(data, ids.var, 
                         outcomes.var, groups.var,
                         covariates.var, 
                         times.dropout.var, times.observation.var,
                         method, dist, model.options)

# summarise the result
summary(fit)
