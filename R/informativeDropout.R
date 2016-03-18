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
                                               dist, model.options, prior.options) {
  
  # make sure we have the correct type of mcmc opts
  if (!is(model.options,"dirichlet.model.options")) {
    stop("Model options error :: options must be of type dirichlet.model.options")
  }
  
  # perform some additional validation on the priors
  if (!is(prior.options,"dirichlet.prior.options")) {
    stop("Prior options error :: options must be of type dirichlet.prior.options")
  }
  # prior for covariate effects
  if (!is.null(covariates.var)) {
    if (is.na(prior.options$betas.covariates.mu) || 
        length(prior.options$betas.covariates.mu) != length(covariates.var)) {
      stop("Prior options error :: invalid prior mean for fixed effects related to covariates")
    }
    if (is.na(prior.options$betas.covariates.R0) || 
        !is.matrix(prior.options$betas.covariates.R0) ||
        nrow(prior.options$betas.covariates.R0) != length(covariates.var) || 
        ncol(prior.options$betas.covariates.R0) != length(covariates.var) ||
        !is.positive.definite(prior.options$betas.covariates.R0)) {
      stop("Prior options error :: invalid prior variance for fixed effects related to covariates")
    }
  }
  
  if (dist == 'gaussian') {
    if (is.null(prior.options$sigma.error.tau)) {
      stop("Prior options error :: invalid tau (hyperparameter) for the error variance")
    }
  }
  
  
  ## reorder the data and cache information on the groups and subjects
  # number of clusters
  n.clusters <- model.options$n.clusters
  # reorder by group, subject, and observation time
  data <- data[order(data[,groups.var], data[,ids.var], data[,times.observation.var]),]
  rownames(data) <- NULL
  n.total = nrow(data)
  # get the list of treatment groups
  groupList <- unique(data[,groups.var])
  # number of subjects
  n.subjects = length(unique(data[,ids.var]))
  # number of subjects per group
  n.perGroup = lapply(groupList, function(group) { 
    # get the covariates
    return(length(unique(data[data[,groups.var] == group, ids.var])))
  })
  # number of observations per subject
  nobj.perGroup = lapply(groupList, function(group) {
    group.data = data[data[,groups.var] == group, ]
    return(sapply(unique(group.data[,ids.var]), function(id) {
      return(nrow(group.data[group.data[,ids.var] == id,]))
    }))
  })
  ## set starting values
  # mixing weights
  if (!is.null(model.options$start.weights.mixing)) {
    weights.mixing = lapply(groupList, function(group) {
      weights = model.options$start.weights.mixing
      return(weights)
    })
  } else {
    weights.mixing = lapply(groupList, function(group) {
      weights = rep(0, n.clusters)
      return(weights)
    })
  }
  
  # number of subjects in each cluster
  n.perCluster = lapply(groupList, function(group) {
    n.subjects = rep(0, n.clusters)
    return(n.subjects)
  }) 
  # Conditional weights -- pr(c_i=h|c_i not in l<h)
  weights.conditional = lapply(groupList, function(group) {
    weights = rep(1/n.clusters, n.clusters)
    weights[n.clusters] = 1 # Apply DP truncation to last cluster
    return(weights)
  })
  
  # initialize the concentration parameters for each group
  concentration = lapply(groupList, function(group) {
    return(prior.options$dp.concentration)
  })
  
  # initialize cluster specific means
  cluster.mu = lapply(groupList, function(group) {
    means = matrix(0,n.clusters,3)
    return(means)
  })
  # initialize the cluster specific variance
  cluster.sigma = lapply(groupList, function(group) {
    covar = prior.options$dp.cluster.sigma
    return(covar)
  })
  # initialize the slope estimates at each droptime
  start.slope.dropoutTimes = NULL
  if (!is.null(model.options$dropout.estimationTimes)) {
    start.slope.dropoutTimes = lapply(groupList, function(group) {
      return (matrix(0,length(model.options$dropout.estimationTimes), n.clusters))
    })
  }
  
  # initialize the regression coefficients for the random intercept
  # and slope.  We append the log of the dropout time
  betas = lapply(groupList, function(group) {
    N = length(unique(data[data[,groups.var] == group, ids.var]))
    logDropout = log(unique(data[data[,groups.var] == group,
                                 c(ids.var,times.dropout.var)])[,times.dropout.var] + 
                              model.options$dropout.offset)
    return(cbind(intercept=rep(0,N), slope=rep(0,N), logDropout))
  })
  # starting value for subject specific deviations from the cluster mean
  betas.deviations = lapply(groupList, function(group) {
    N = length(unique(data[data[,groups.var] == group, ids.var]))
    return(cbind(intercept.dev=rep(0,N), slope.dev=rep(0,N), logDropout.dev=rep(0,N)))
  })
  
  # covariate coefficients 
  start.betas.covariates = NULL
  if (!is.null(covariates.var)) {
    if (!is.null(prior.options$betas.covariates)) {
      start.betas.covariates = prior.options$betas.covariates
    } else {
      start.betas.covariates = getInitialEstimatesCovariates(dist, 
                                                             as.data.frame(data[,covariates.var]), 
                                                             data[,outcomes.var])
    }
  }
  
  # starting values for the mean and variance of the Dirichlet baseline distribution
  dp.dist.mu0 = lapply(groupList, function(group) {
    return(prior.options$dp.dist.mu0)
  })
  dp.dist.sigma0 = lapply(groupList, function(group) {
    return(prior.options$dp.dist.sigma0)
  })
  # starting values for the common cluster variance
  dp.cluster.sigma = lapply(groupList, function(group) {
    return(prior.options$dp.cluster.sigma)
  })
  
  # for Gaussian outcomes, the residual error
  if (dist == "gaussian") {
    sigma.error = prior.options$sigma.error
  }
  
  
  # initialize the first model iteration
  modelIterationList <- vector(mode = "list", length = model.options$iterations)
  modelIterationList[[1]] = dirichlet.iteration(
    weights.mixing=weights.mixing, 
    weights.conditional=weights.conditional,
    betas = betas, 
    betas.deviations = betas.deviations,
    betas.covariates = start.betas.covariates, 
    betas.covariates.mu = prior.options$betas.covariates.mu,
    betas.covariates.sigma = prior.options$betas.covariates.sigma,
    dp.dist.mu0 = dp.dist.mu0, dp.dist.sigma0 = dp.dist.sigma0,
    dp.cluster.sigma = dp.cluster.sigma,
    dp.concentration=concentration, 
    cluster.N = n.perCluster,
    cluster.mu=cluster.mu, 
    expected.intercept=NULL, expected.slope=NULL,
    slope.dropoutTimes = start.slope.dropoutTimes,
    sigma.error = sigma.error
  )
  
  for (i in 2:model.options$iterations) {
    
    print(paste("ITER = ", i, sep=""))
    
    model.previous = modelIterationList[[i-1]]
    # make a copy which will be modified as we move through the iteration
    model.current = model.previous
    
    for (group.index in 1:length(groupList)) {
      
      group = groupList[group.index]
      # get the subset of data for this group
      group.data = data[data[,groups.var] == group,]
      # convenience variables for the group specific information
      group.n = n.perGroup[[group.index]]
      group.nobs = nobj.perGroup[[group.index]]
      group.weights.conditional = model.current$weights.conditional[[group.index]]
      group.betas = model.current$betas[[group.index]]
      group.betas.deviations = model.current$betas.deviations[[group.index]]
      group.cluster.mu = model.current$cluster.mu[[group.index]]
      group.dp.cluster.sigma = model.current$dp.cluster.sigma[[group.index]]
      group.dp.dist.mu0 = model.current$dp.dist.mu0[[group.index]]
      group.dp.dist.sigma0 = model.current$dp.dist.sigma0[[group.index]]
      group.dp.cluster.sigma = model.current$dp.cluster.sigma[[group.index]]
      
      # calculate the overall probability of belonging to a cluster (length = #clusters)
      cumulative.weights.conditional <- cumprod(1 - group.weights.conditional)
      model.current$weights.mixing[[group.index]] = sapply(1:n.clusters, function(i) {
        if (i <= 1) {
          return(group.weights.conditional[i])
        } else {
          return(group.weights.conditional[i]*cumulative.weights.conditional[i-1])
        }
      })
      group.weights.mixing = model.current$weights.mixing[[group.index]]
      # calculate the probability that each individual belongs to a cluster (#subjects x #clusters)
      # overall cluster prob * normal density (mean by cluster, common covariance) 
      # at subject specific int,slope, log(u)
      clusterProbabilities = sapply(1:n.clusters, function(h) {
        return(group.weights.mixing[h] * dmvnorm(group.betas,
                                                 group.cluster.mu[h,],
                                                 group.dp.cluster.sigma))
      })
      # Build multinomial probabilities (Equation from step 1, Appendix A)
      prob = (clusterProbabilities / apply(clusterProbabilities,1,sum))
      # draw from multinomial to determine which cluster the individual belongs to
      # C is the cluster indicator
      group.cluster.assignments = rMultinom(prob, 1)
      model.current$cluster.assignments[[group.index]] = group.cluster.assignments
      
      # calculate/save the number of people in each cluster
      # Must allow zeros for empty clusters
      model.current$cluster.N[[group.index]] = sapply(1:n.clusters, function(h) {
        return(length(group.cluster.assignments[group.cluster.assignments==h]))
      })
      group.cluster.N = model.current$cluster.N[[group.index]]
      
      # Update the stick breaking weights, drawn from beta distributions (Step 2 appendix A)
      model.current$weights.conditional[[group.index]] = sapply(1:(n.clusters), function(h) {
        if (h < n.clusters) {
          return(rbeta(1, 1 + group.cluster.N[h], 
                       prior.options$dp.concentration + 
                         sum(group.cluster.N[(h+1):n.clusters])))
        } else {
          return(1)
        }
      })
      
      # number of clusters that have at least one subject (common to have empty clusters)
      numNonEmptyClusters = length(group.cluster.N[group.cluster.N != 0])
      
      # Alpha is concentration parameter of the dirichlet process (step 7)
      # sample latent beta distributed variable (Escobar and West 1995)
      eta <- rbeta(1, prior.options$dp.concentration + 1, group.n)
      # sample the concentration parameter from a mixture of Gamma distributions
      pi.eta <- ((prior.options$dp.concentration.alpha + numNonEmptyClusters - 1) /
                   (group.n * (prior.options$dp.concentration.beta - log(eta)) +
                      prior.options$dp.concentration.alpha + numNonEmptyClusters - 1))
      # now draw the new concentration parameter for the Dirichlet process
      model.current$dp.concentration[[group.index]] <- 
        ((pi.eta * rgamma(1, prior.options$dp.concentration.alpha + numNonEmptyClusters, 
                          prior.options$dp.concentration.beta - log(eta))) + 
           ((1 - pi.eta) * rgamma(1, prior.options$dp.concentration.alpha + numNonEmptyClusters - 1, 
                                  prior.options$dp.concentration.beta - log(eta))))
      
      # Update cluster means and covariances (density estimate)
      # assumes common covariance across clusters but potentially different means
      dp.dist.sigma0.inv = solve(group.dp.dist.sigma0)
      dp.cluster.sigma.inv = solve(group.dp.cluster.sigma)
      
      for (h in 1:n.clusters) {
        # calculate the means per cluster (Step 3)
        S <- solve(dp.dist.sigma0.inv + group.cluster.N[h] * dp.cluster.sigma.inv)
        if (group.cluster.N[h] == 1) {
          m <- S %*% (dp.dist.sigma0.inv %*% group.dp.dist.mu0 + 
                        dp.cluster.sigma.inv %*% ((group.betas[group.cluster.assignments==h,])) )
        } else {
          m <- S %*%(dp.dist.sigma0.inv %*% group.dp.dist.mu0 + 
                       dp.cluster.sigma.inv %*% (colSums(group.betas[group.cluster.assignments==h,])) )
        }
        means.currentCluster = rmvnorm(1, m, S)
        model.current$cluster.mu[[group.index]][h,] = means.currentCluster
        
        # All subjects within a cluster have same distribution of random effects
        # Step 5: draw intercept conditional on slope and dropout time
        #Update beta i given cluster mean and var and ui for each subject
        # B0|B1 and U
        # inverse of the 2x2 covariance of the slope and the log dropout time
        inv.slopeU = solve(group.dp.cluster.sigma[-1,-1])
        covar = matrix(group.dp.cluster.sigma[1, c(2,3)],1)
        betas.currentCluster = group.betas[group.cluster.assignments==h,,drop=FALSE]
        betas.currentCluster.rows = ifelse(group.cluster.N[h] > 0, nrow(betas.currentCluster), 0)
        subjectMeans.currentCluster = matrix(rep(means.currentCluster, betas.currentCluster.rows),
                                             ncol=3,byrow=TRUE)
        ## conditional prior mean/var
        prior.mean = as.numeric((means.currentCluster[1] + covar %*% inv.slopeU %*% 
                                   t(betas.currentCluster[,c(2,3),drop=FALSE] - subjectMeans.currentCluster[,c(2,3),drop=FALSE])))
        prior.var = as.numeric((group.dp.cluster.sigma[1,1] - covar %*% inv.slopeU %*% t(covar)))
        ## calculate the posterior mean and variance for the random intercept
        posterior.var = 1 / ((1 / prior.var) + (group.nobs[group.cluster.assignments==h] / model.current$sigma.error))
        # get the data for subjects in this cluster
        tempid = rep(group.cluster.assignments, group.nobs)
        data.currentCluster = group.data[tempid==h,]
        # calculate the posterior mean
        if (is.null(covariates.var)) {
          randomInts = data.currentCluster[,outcomes.var] - 
            data.currentCluster[, times.observation.var] * 
            rep(group.betas[group.cluster.assignments==h, 2],
                group.nobs[group.cluster.assignments==h])
        } else {
          randomInts = data.currentCluster[,outcomes.var] - 
            data.currentCluster[, times.observation.var] * 
            rep(group.betas[group.cluster.assignments==h, 2],
                group.nobs[group.cluster.assignments==h]) - 
            data.currentCluster[, covariates.var] %*% group.beta.covariates
        }
        posterior.mean = (posterior.var * 
                            ((prior.mean/prior.var) + 
                               ((1/sigma.error) * 
                                  tapply(randomInts, data.currentCluster[, ids.var],sum)))) 
        # draw the random intercepts
        group.betas[group.cluster.assignments == h, 1] = 
          rnorm(group.cluster.N[h], posterior.mean, sqrt(posterior.var))
        
        ######  B1|B0 and U (Step 5 continued, draw random slope given intercept and u)
        inv.intU = solve(group.dp.cluster.sigma[-2,-2])
        covar = matrix(c(group.dp.cluster.sigma[1,2], group.dp.cluster.sigma[2,3]),1)
        # conditional prior mean/var
        prior.mean = as.numeric((means.currentCluster[2] + covar %*% inv.intU %*% 
                                   t(betas.currentCluster[,c(1,3),drop=FALSE] - subjectMeans.currentCluster[,c(1,3),drop=FALSE])))
        prior.var = as.numeric((group.dp.cluster.sigma[2,2] - covar %*% inv.intU %*% t(covar)))
        ## calculate the posterior mean and variance for the random intercept
        posterior.var = 1 / ((1 / prior.var) + 
                               (tapply(data.currentCluster[, times.observation.var]^2, 
                                       data.currentCluster[, ids.var], sum) / 
                                  sigma.error))
        # calculate the posterior mean
        if (is.null(covariates.var)) {
          randomSlopes = (data.currentCluster[, times.observation.var]) * 
            (data.currentCluster[,outcomes.var] - 
               rep(group.betas[group.cluster.assignments==h, 1],
                   group.nobs[group.cluster.assignments==h]))
        } else {
          randomSlopes = (data.currentCluster[, times.observation.var]) * 
            (data.currentCluster[,outcomes.var] - 
               rep(group.betas[group.cluster.assignments==h, 1],
                   group.nobs[group.cluster.assignments==h]) - 
               data.currentCluster[, covariates.var] %*% group.beta.covariates)
        }
        posterior.mean = (posterior.var * 
                            ((prior.mean/prior.var) + 
                               ((1/sigma.error) * 
                                  tapply(randomSlopes, data.currentCluster[, ids.var],sum)))) 
        # draw the random intercepts
        group.betas[group.cluster.assignments == h, 2] = 
          rnorm(group.cluster.N[h], posterior.mean, sqrt(posterior.var))
        
        # calculating the subject specific deviations from the cluster means
        # used to simplify calculation of the covariance
        #calculate betas minus their means
        group.betas.deviations[group.cluster.assignments == h,] = (
          group.betas[group.cluster.assignments == h,] - 
            (matrix(rep(means.currentCluster, group.cluster.N[h]), ncol=3, byrow=TRUE))
        )
        
      } # END CLUSTER-SPECIFIC UPDATE LOOP
      
      # update the model iteration
      model.current$betas[[group.index]] = group.betas
      model.current$betas.deviations[[group.index]] = group.betas.deviations

      #### TODO: add flag to indicate if cluster covariance is updated per group
      ## or updated across groups
      # Common covariance for each cluster (Step 4)
      # Update Cluster Covariance
      model.current$dp.cluster.sigma[[group.index]] = riwish(
        prior.options$dp.cluster.sigma.nu0 + group.n, # <- check typo
        prior.options$dp.cluster.sigma.T0 + crossprod(group.betas.deviations)
      )
      
      ### Update the hyprparameters of the baseline distribution of 
      ### the Dirichlet process (Step 6)
      # update the mean
      Sb.inv = solve(prior.options$dp.dist.mu0.Sb)
      var = solve(Sb.inv + numNonEmptyClusters * dp.dist.sigma0.inv)
      if (numNonEmptyClusters == 1) {
        m <- var %*% (Sb.inv %*% prior.options$dp.dist.mu0.mb + 
                        dp.dist.sigma0.inv %*% ((model.current$cluster.mu[[group.index]][group.cluster.N > 0,])) )
      } else {
        m <- var %*% (prior.options$dp.dist.mu0.Sb %*% prior.options$dp.dist.mu0.mb + 
                        dp.dist.sigma0.inv %*% (colSums(model.current$cluster.mu[[group.index]][group.cluster.N > 0,])) )
      }
      model.current$dp.dist.mu0[[group.index]] = as.vector(rmvnorm(1,m,var))
      
      # update the variance of the baseline distribution
      model.current$dp.dist.sigma0[[group.index]] = riwish(
        prior.options$dp.dist.sigma0.nub + numNonEmptyClusters,
        prior.options$dp.dist.sigma0.Tb + 
          crossprod(model.current$cluster.mu[[group.index]][group.cluster.N > 0,] - 
                      matrix(rep(model.current$dp.dist.mu0[[group.index]],numNonEmptyClusters), 
                             nrow=numNonEmptyClusters, byrow=T))
      )
      
      
      ### calculate per-group summary statistics
      #Estimate Density of Slope - for plotting density of the random slope (output only)  
      domain <-seq(-4,4,0.01)
      sim.slope = vector()
      sim.int = vector()
      for (h in 1:n.clusters) {
        sim.int <-rbind(sim.int,
                        (clusterProbabilities[h] * 
                           dnorm(domain, model.current$cluster.mu[[group.index]][h,1], 
                                 sqrt(model.current$dp.dist.sigma0[[group.index]][1,1]))))
        sim.slope <-rbind(sim.slope,
                          (clusterProbabilities[h] * 
                             dnorm(domain, model.current$cluster.mu[[group.index]][h,2], 
                                   sqrt(model.current$dp.dist.sigma0[[group.index]][2,2]))))
      }
      model.current$density.intercept[[group.index]] <- colSums(sim.int)
      model.current$density.slope[[group.index]] <- colSums(sim.slope)
      
      # Equation 19
      #Estimate slope at each dropout time
      if (!is.null(model.options$dropout.estimationTimes)) {
        tmp = matrix(0,length(model.options$dropout.estimationTimes), n.clusters)
        for (h in 1:n.clusters) {
          tmp[,h] = (
            clusterProbabilities[h] * 
              dnorm(log(model.options$dropout.estimationTimes),
                    model.current$cluster.mu[[group.index]][h,3],
                    sqrt(model.current$dp.cluster.sigma[[group.index]][3,3]))
          ) 
        }
        
        p = tmp / apply(tmp,1,sum)
        for(u in 1:length(model.options$dropout.estimationTimes)) {
          covar = model.current$dp.cluster.sigma[[group.index]]
          dropoutTime = model.options$dropout.estimationTimes[u]
          model.current$slope.dropoutTimes[[group.index]][u] = (
            sum(p[u,] * 
                  (model.current$cluster.mu[[group.index]][,2] + 
                     covar[2,3] * 
                     (log(u) - model.current$cluster.mu[[group.index]][,2]) /
                     covar[3,3]))
          )
        }
      } else {
        model.current$slope.dropoutTimes = NULL
      }
      
      # Estimate E(B1), E(B0) - section 2.4 expectation of random effects (int, slope)
      model.current$expected.intercept[[group.index]] = 
        sum(clusterProbabilities * model.current$cluster.mu[[group.index]][,1])
      model.current$expected.slope[[group.index]] = 
        sum(clusterProbabilities * model.current$cluster.mu[[group.index]][,2])
      
      
    } # END GROUP-SPECIFIC UPDATE LOOP
    
    # combine the random effects for each group into complete arrays
    intercepts = vector()
    slopes = vector()
    for (group.index in 1:length(groupList)) {
      intercepts = c(intercepts, rep(model.current$betas[[group.index]][,1], 
                                     nobj.perGroup[[group.index]]))
      slopes = c(slopes, rep(model.current$betas[[group.index]][,2], 
                             nobj.perGroup[[group.index]]))
    }
    
    if (dist == "gaussian") {
      ## TODO: remove group specific variables, use betas (random effects)
      # Step 9 Update sigma.error with inverse gamma
      residual = data[, outcomes.var] - (intercepts + data[, times.observation.var] * slopes)
      if (!is.null(covariates.var)) {
        residual = residual - as.vector(as.matrix(data[,covariates.var]) %*% 
                                          model.current$betas.covariates)
      }
      g <- prior.options$sigma.error.tau + crossprod(residual)/2
      tau <- rgamma(1, prior.options$sigma.error.tau + n.total / 2, g)
      model.current$sigma.error = 1 / tau
    } 
    
    # update fixed effects associated with covariates
    if (!is.null(covariates.var)) {
      if (dist == "gaussian") {
        sigma.error.inv = 1/model.current$sigma.error
        prior.sigma.inv = solve(prior.options$betas.covariates.sigma)
        
        residuals = data[, outcomes.var] - intercepts - slopes * data[, times.observation.var]
        var = solve(prior.sigma.inv + crossprod(data[,covariates.var]) * sigma.error.inv)
        
        m <- var %*% (crossprod(data[, covariates.var], residuals) * sigma.error.inv)
        model.current$betas.covariates = rmvnorm(1,m,var)
        
      } else {
        # metropolis hastings
      }
    }
    
    # save the current iteration
    modelIterationList[[i]] = model.current
    
  } # END ITERATION LOOP
  
  return(modelIterationList)
} # END FUNCTION 




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
  
  # for binary case, initialize eta
  if (dist == 'binary') {
    eta.null = lapply(subjectsPerGroup, function(N) {
      return(logit(rep(sum(y)/N, N)))
    })
  }
  
  
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






