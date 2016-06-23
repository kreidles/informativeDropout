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
# Dirichlet process model examples
#####################################

# load the simulated Multicenter AIDS Cohort Study data
data(macs)

# create the model options
model.options=dirichlet.model.options(iterations=100, n.clusters=60, burnin=0, thin=1,
                                      print=10,
                                      dropout.offset=0,
                                      dropout.estimationTimes = seq(2,13,1),
                                      dp.concentration=1,
                                      dp.concentration.alpha=1,
                                      dp.concentration.beta=1,
                                      dp.cluster.sigma = diag(3),
                                      dp.cluster.sigma.nu0 = 5,
                                      dp.cluster.sigma.T0 = diag(3),
                                      dp.dist.mu0 = c(0,0,log(7)),
                                      dp.dist.mu0.mb = c(0,0,log(7)),
                                      dp.dist.mu0.Sb = diag(3),
                                      dp.dist.sigma0 = diag(3),
                                      dp.dist.sigma0.nub = 5,
                                      dp.dist.sigma0.Tb = diag(3),
                                      betas.covariates = c(1.0, 0.1),
                                      betas.covariates.mu = c(0,0),
                                      betas.covariates.sigma = 100*diag(2),
                                      sigma.error = 1,
                                      sigma.error.tau = 0.001)


# define the columns in the data set to be used in the model
data=macs
ids.var="CASEID"
outcomes.var="logcd4"
groups.var="hard"
covariates.var=c("logbase", "basebytime") 
times.dropout.var="dropouttime"
times.observation.var='time'
method="dirichlet"
dist='gaussian'

# set a random seed if desired
set.seed(1066)

# fit the Dirichlet process model
fit = informativeDropout(data, ids.var, 
                         outcomes.var, groups.var,
                         covariates.var, 
                         times.dropout.var, times.observation.var,
                         method, dist, model.options)

# print a summary
summary(fit)

