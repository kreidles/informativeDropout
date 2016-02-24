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
                                               dist, mcmc.options, prior.options, 
                                               startValues=NULL) {
  
  # validate the mcmc options
  if (is.na(mcmc.options$iterations) || mcmc.options$iterations <= 1) {
    stop("MCMC options error :: invalid number of iterations")
  }
  if (is.na(mcmc.options$burnIn) || mcmc.options$burnIn >= mcmc.options$iterations) {
    stop("MCMC options error :: the burn in period must be less than the total number of iterations")
  }
  if (is.na(mcmc.options$numClusters) || mcmc.options$numClusters <= 1) {
    stop("MCMC options error :: invalid number of clusters")
  }
  
  # validate the prior options
  # Dirichlet process details
  if (is.na(prior.options$dp.concentration) || prior.options$dp.concentration <= 0) {
    stop("Prior options error :: invalid concentration parameter for the Dirichlet process")
  }
  if (is.na(prior.options$dp.concentration.alpha) || prior.options$dp.concentration.alpha <= 0) {
    stop("Prior options error :: invalid gamma prior for the Dirichlet process concentration parameter")
  }
  if (is.na(prior.options$dp.concentration.beta) || prior.options$dp.concentration.beta <= 0) {
    stop("Prior options error :: invalid gamma prior for the Dirichlet process concentration parameter")
  }
  # priors for baseline distribution of the Dirichlet process
  # mean of the baseline distribution
  if (is.na(prior.options$dp.dist.mu0)) {
    stop("Prior options error :: invalid mean for the Dirichlet process baseline distribution")
  }
  # hyperparams for the mean of the baseline distribution
  if (is.na(prior.options$dp.dist.mu0.mb)) {
    stop("Prior options error :: invalid mean (hyperparameter) for the mean of the Dirichlet process baseline distribution")
  }
  if (is.na(prior.options$dp.dist.mu0.Sb) || 
      !is.matrix(prior.options$dp.dist.mu0.Sb) ||
      nrow(prior.options$dp.dist.mu0.Sb) != 3 || ncol(prior.options$dp.dist.mu0.Sb) != 3 ||
      !is.positive.definite(prior.options$dp.dist.mu0.Sb)) {
    stop("Prior options error :: invalid variance (hyperparameter) for the mean of the Dirichlet process baseline distribution")
  }
  
  # variance of the baseline distribution
  if (is.na(prior.options$dp.dist.sigma0)) {
    stop("Prior options error :: no variance for the Dirichlet process baseline distribution")
  }
  # hyperparameters of the variance of the baseline distribution
  if (is.na(prior.options$dp.dist.sigma0.nub)) {
    stop("Prior options error :: invalid df (hyperparameter) for the variance of the Dirichlet process baseline distribution")
  }
  if (is.na(prior.options$dp.dist.sigma0.Tb) || 
      !is.matrix(prior.options$dp.dist.sigma0.Tb) ||
      nrow(prior.options$dp.dist.sigma0.Tb) != 3 || ncol(prior.options$dp.dist.sigma0.Tb) != 3 ||
      !is.positive.definite(prior.options$dp.dist.sigma0.Tb)) {
    stop("Prior options error :: invalid scale matrix (hyperparameter) for the variance of the Dirichlet process baseline distribution")
  }
  
  # cluster specific variance for the Dirichlet process
  if (is.na(prior.options$dp.cluster.sigma) || 
      !is.matrix(prior.options$dp.cluster.sigma) ||
      nrow(prior.options$dp.cluster.sigma) != 3 || ncol(prior.options$dp.cluster.sigma) != 3 ||
      !is.positive.definite(prior.options$dp.cluster.sigma)) {
    stop("Prior options error :: invalid cluster-specific variance for the Dirichlet process")
  }
  # hyperparameters of the cluster specific variance
  if (is.na(prior.options$dp.cluster.sigma.nu0)) {
    stop("Prior options error :: invalid df (hyperparameter) for the cluster-specific variance of the Dirichlet process")
  }
  if (is.na(prior.options$dp.cluster.sigma.T0) || 
      !is.matrix(prior.options$dp.cluster.sigma.T0) ||
      nrow(prior.options$dp.cluster.sigma.T0) != 3 || ncol(prior.options$dp.cluster.sigma.T0) != 3 ||
      !is.positive.definite(prior.options$dp.cluster.sigma.T0)) {
    stop("Prior options error :: invalid scale matrix (hyperparameter) for the cluster-specific variance of the Dirichlet process")
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
    if (is.na(prior.options$sigma.error.tau1)) {
      stop("Prior options error :: invalid tau1 (hyperparameter) for the error variance")
    }
    if (is.na(prior.options$sigma.error.tau2)) {
      stop("Prior options error :: invalid tau2 (hyperparameter) for the error variance")
    }
  }
  
  ## reorder the data and cache information on the groups and subjects
  # number of clusters
  numClusters <- mcmc.options$numClusters
  # reorder by group, subject, and observation time
  data <- data[order(data[,groups.var], data[,ids.var], data[,times.observation.var]),]
  rownames(data) <- NULL
  # get the list of treatment groups
  groupList <- unique(data[,groups.var])
  # number of subjects
  numSubjects = length(unique(data[,ids.var]))
  # number of subjects per group
  subjectsPerGroup = lapply(groupList, function(group) { 
    # get the covariates
    return(length(unique(data[data[,groups.var] == group, ids.var])))
  })
  # number of observations per subject
  numObservationsPerGroup = lapply(groupList, function(group) {
    group.data = data[data[,groups.var] == group, ]
    return(sapply(unique(group.data[,ids.var]), function(id) {
      return(nrow(group.data[group.data[,ids.var] == id,]))
    }))
  })
  ## set starting values
  # mixing weights
  mixingWeights = lapply(groupList, function(group) {
    weights = rep(0, numClusters)
    return(weights)
  })
  # number of subjects in each cluster
  subjectsPerCluster = lapply(groupList, function(group) {
    numSubjects = rep(0, numClusters)
    return(numSubjects)
  }) 
  # Conditional weights -- pr(c_i=h|c_i not in l<h)
  conditionalWeights = lapply(groupList, function(group) {
    weights = rep(1/numClusters, numClusters)
    weights[numClusters] = 1 # Apply DP truncation to last cluster
    return(weights)
  })
  
  # initialize the concentration parameters for each group
  concentration = lapply(groupList, function(group) {
    return(prior.options$dp.concentration)
  })

  # initialize cluster specific means
  clusterMeans = lapply(groupList, function(group) {
    means = matrix(0,numClusters,3)
    return(means)
  })
  # initialize the cluster specific variance
  clusterCovariance = lapply(groupList, function(group) {
    covar = diag(3)
    return(covar)
  })
  
  # initialize the regression coefficients for the random intercept
  # and slope.  We append the log of the 
  betas = lapply(groupList, function(group) {
    N = length(unique(data[data[,groups.var] == group, ids.var]))
    logDropout = log(unique(data[data[,groups.var] == group, ids.var]))
    return(cbind(intercept=rep(0,N), slope=rep(0,N), logDropout))
  })
  # starting value for subject specific deviations from the cluster mean
  betas.deviations = lapply(groupList, function(group) {
    N = length(unique(data[data[,groups.var] == group, ids.var]))
    return(cbind(intercept.dev=rep(0,N), slope.dev=rep(0,N), logDropout.dev=rep(0,N)))
  })
  
  # covariate coefficients -- TODO
  betaCovariate.init = NULL
  if (!is.null(covariates.var)) {
    betaCovariate.init = getInitialEstimatesCovariates(dist, 
                                                       as.data.frame(data[,covariates.var]), 
                                                       data[,outcomes.var])
  }

  
  
  # what is this?
  # posterior probabilities
  posteriorClusterProb = lapply(subjectsPerGroup, function(nPerGroup) {
    return(matrix(0, nPerGroup, numClusters))
  })

  # starting values for the mean and variance of the Dirichlet baseline distribution
  dp.dist.mu0 = lapply(groupList, function(group) {
    return(c(0,0,0.5))
  })
  dp.dist.sigma0 = lapply(groupList, function(group) {
    covar = diag(3)
    return(covar)
  })
  # starting values for the common cluster variance
  dp.cluster.sigma = lapply(groupList, function(group) {
    covar = diag(3)
    return(covar)
  })
  
  # for Gaussian outcomes, the residual error
  sigma.error = lapply(groupList, function(group) {
    return(1)
  })
  
  # initialize the first model iteration
  modelIterationList <- vector(mode = "list", length = mcmc.options$iterations)
  modelIterationList[[1]] = dirichlet.iteration(
    mixingWeights=mixingWeights, 
    conditionalWeights=conditionalWeights,
    betas = betas, betas.covariates = betaCovariate.init,
    betas.deviations = betas.deviations,
    perClusterN=subjectsPerCluster,
    dp.dist.mu0 = dp.dist.mu0,
    dp.dist.sigma0 = dp.dist.sigma0,
    dp.cluster.sigma = dp.cluster.sigma,
    clusterMeans=clusterMeans, clusterCovariance=clusterCovariance, 
    expected.intercept=NULL, expected.slope=NULL,
    concentration=concentration, 
    sigma.error = sigma.error
  )

  for (i in 2:mcmc.options$iterations) {
    
    print(paste("ITER = ", i, sep=""))
    
    model.previous = modelIterationList[[i-1]]
    # make a copy which will be modified as we move through the iteration
    model.current = model.previous
    
    for (group.index in groupList) {
      
      group = groupList[group.index]
      # get the subset of data for this group
      group.data = data[data[,groups.var] == group,]
      # convenience variables for the group specific information
      group.N = subjectsPerGroup[[group.index]]
      group.nobs = numObservationsPerGroup[[group.index]]
      group.conditionalWeights = model.current$conditionalWeights[[group.index]]
      group.betas = model.current$betas[[group.index]]
      group.betas.deviations = model.current$betas.deviations[[group.index]]
      group.beta.covariates = model.current$beta.covariates[[group.index]]
      group.clusterMeans = model.current$clusterMeans[[group.index]]
      group.clusterCovariance = model.current$clusterCovariance[[group.index]]
      group.dp.dist.mu0 = model.current$dp.dist.mu0[[group.index]]
      group.dp.dist.sigma0 = model.current$dp.dist.sigma0[[group.index]]
      group.dp.cluster.sigma = model.current$dp.cluster.sigma[[group.index]]
      group.sigma.error = model.current$sigma.error[[group.index]]
      # calculate the overall probability of belonging to a cluster (length = #clusters)
      cumulative.conditionalWeights <- cumprod(1 - group.conditionalWeights)
      model.current$mixingWeights[[group.index]] = sapply(1:numClusters, function(i) {
        if (i <= 1) {
          return(group.conditionalWeights[i])
        } else {
          return(group.conditionalWeights[i]*cumulative.conditionalWeights[i-1])
        }
      })
      group.mixingWeights = model.current$mixingWeights[[group.index]]
      
      # calculate the probability that each individual belongs to a cluster (#subjects x #clusters)
      # overall cluster prob * normal density (mean by cluster, common covariance) 
      # at subject specific int,slope, log(u)
      clusterProbabilities = sapply(1:numClusters, function(h) {
        return(group.mixingWeights[h] * dmvnorm(group.betas,
                                                group.clusterMeans[h,],
                                                group.clusterCovariance))
      })
      # Build multinomial probabilities (Equation from step 1, Appendix A)
      prob = (clusterProbabilities / apply(clusterProbabilities,1,sum))
      # draw from multinomial to determine which cluster the individual belongs to
      # C is the cluster indicator
      group.clusterAssignments = rMultinom(prob, 1)
      model.current$clusterAssignments[[group.index]] = group.clusterAssignments
      
      # calculate/save the number of people in each cluster
      # Must allow zeros for empty clusters
      model.current$perClusterN[[group.index]] = sapply(1:numClusters, function(h) {
        return(length(group.clusterAssignments[group.clusterAssignments==h]))
      })
      perClusterN = model.current$perClusterN[[group.index]]
      
      # Update the stick breaking weights, drawn from beta distributions (Step 2 appendix A)
      model.current$conditionalWeights[[group.index]] = sapply(1:(numClusters), function(h) {
        if (h < numClusters) {
          return(rbeta(1, 1 + perClusterN[h], 
                       prior.options$dp.concentration + sum(perClusterN[(h+1):numClusters])))
        } else {
          return(rbeta(1, 1 + perClusterN[h], prior.options$dp.concentration))
        }
      })
      
      # number of clusters that have at least one subject (common to have empty clusters)
      numNonEmptyClusters = length(perClusterN[perClusterN != 0])

      # Alpha is concentration parameter of the dirichlet process (step 7)
      # sample latent beta distributed variable (Escobar and West 1995)
      eta <- rbeta(1, prior.options$dp.concentration + 1, group.N)
      # sample the concentration parameter from a mixture of Gamma distributions
      pi.eta <- ((prior.options$dp.concentration.alpha + numNonEmptyClusters - 1) /
                   (group.N * (prior.options$dp.concentration.beta - log(eta)) +
                      prior.options$dp.concentration.alpha + numNonEmptyClusters - 1))
      # now draw the new concentration parameter for the Dirichlet process
      model.current$concentration[[group.index]] <- 
        ((pi.eta * rgamma(1, prior.options$dp.concentration.alpha + numNonEmptyClusters, 
                         prior.options$dp.concentration.beta - log(eta))) + 
           ((1 - pi.eta) * rgamma(1, prior.options$dp.concentration.alpha + numNonEmptyClusters - 1, 
                                  prior.options$dp.concentration.beta - log(eta))))

      # Update cluster means and covariances (density estimate)
      # assumes common covariance across clusters but potentially different means
      dp.dist.sigma0.inv = solve(group.dp.dist.sigma0)
      dp.cluster.sigma.inv = solve(group.dp.cluster.sigma)
      
      for (h in 1:numClusters) {
        cat ("CLUSTER: ", h, "\n")
        # calculate the means per cluster (Step 3)
        S <- solve(dp.dist.sigma0.inv + perClusterN[h] * dp.cluster.sigma.inv)
        if (perClusterN[h] == 1) {
          m <- S %*% (dp.dist.sigma0.inv %*% group.dp.dist.mu0 + 
                        dp.cluster.sigma.inv %*% ((group.betas[group.clusterAssignments==h,])) )
        } else {
          m <- S %*%(dp.dist.sigma0.inv %*% group.dp.dist.mu0 + 
                       dp.cluster.sigma.inv %*% (colSums(group.betas[group.clusterAssignments==h,])) )
        }
        means.currentCluster = rmvnorm(1, m, S)
        model.current$clusterMeans[[group.index]][h,] = means.currentCluster

        # All subjects within a cluster have same distribution of random effects
        # Step 5: draw intercept conditional on slope and dropout time
        #Update beta i given cluster mean and var and ui for each subject
        # B0|B1 and U
        # inverse of the 2x2 covariance of the slope and the log dropout time
        inv.slopeU = solve(group.dp.cluster.sigma[-1,-1])
        covar = matrix(group.dp.cluster.sigma[1, c(2,3)],1)
        betas.currentCluster = group.betas[group.clusterAssignments==h,]
        subjectMeans.currentCluster = matrix(rep(means.currentCluster, 
                                                 nrow(betas.currentCluster)),
                                             ncol=3,byrow=TRUE)
        ## conditional prior mean/var
        prior.mean = as.numeric((means.currentCluster[1] + covar %*% inv.slopeU %*% 
                t(betas.currentCluster[,c(2,3)] - subjectMeans.currentCluster[,c(2,3)])))
        prior.var = as.numeric((group.dp.cluster.sigma[1,1] - covar %*% inv.slopeU %*% t(covar)))
        ## calculate the posterior mean and variance for the random intercept
        posterior.var = 1 / ((1 / prior.var) + (group.nobs[group.clusterAssignments==h] / group.sigma.error))
        # get the data for subjects in this cluster
        tempid = rep(group.clusterAssignments, group.nobs)
        data.currentCluster = group.data[tempid==h,]
        # calculate the posterior mean
        if (is.null(covariates.var)) {
          randomInts = data.currentCluster[,outcomes.var] - 
            data.currentCluster[, times.observation.var] * 
            rep(group.betas[group.clusterAssignments==h, 2],
                group.nobs[group.clusterAssignments==h])
        } else {
          randomInts = data.currentCluster[,outcomes.var] - 
            data.currentCluster[, times.observation.var] * 
            rep(group.betas[group.clusterAssignments==h, 2],
                group.nobs[group.clusterAssignments==h]) - 
            data.currentCluster[, covariates.var] %*% group.beta.covariates
        }
        posterior.mean = (posterior.var * 
                            ((prior.mean/prior.var) + 
                               ((1/group.sigma.error) * 
                               tapply(randomInts, data.currentCluster[, ids.var],sum)))) 
        # draw the random intercepts
        group.betas[group.clusterAssignments == h, 1] = 
          rnorm(perClusterN[h], posterior.mean, sqrt(posterior.var))
        
        ######  B1|B0 and U (Step 5 continued, draw random slope given intercept and u)
        inv.intU = solve(group.dp.cluster.sigma[-2,-2])
        covar = matrix(c(group.dp.cluster.sigma[1,2], group.dp.cluster.sigma[2,3]),1)
        # conditional prior mean/var
        prior.mean = as.numeric((means.currentCluster[2] + covar %*% inv.intU %*% 
                                   t(betas.currentCluster[,c(1,3)] - subjectMeans.currentCluster[,c(1,3)])))
        prior.var = as.numeric((group.dp.cluster.sigma[2,2] - covar %*% inv.intU %*% t(covar)))
        ## calculate the posterior mean and variance for the random intercept
        posterior.var = 1 / ((1 / prior.var) + 
                               (tapply(data.currentCluster[, times.observation.var]^2, 
                                       data.currentCluster[, ids.var], sum) / 
                                  group.sigma.error))
        # calculate the posterior mean
        if (is.null(covariates.var)) {
          randomSlopes = (data.currentCluster[, times.observation.var]) * 
            (data.currentCluster[,outcomes.var] - 
            rep(group.betas[group.clusterAssignments==h, 1],
                group.nobs[group.clusterAssignments==h]))
        } else {
          randomSlopes = (data.currentCluster[, times.observation.var]) * 
            (data.currentCluster[,outcomes.var] - 
               rep(group.betas[group.clusterAssignments==h, 1],
                   group.nobs[group.clusterAssignments==h]) - 
            data.currentCluster[, covariates.var] %*% group.beta.covariates)
        }
        posterior.mean = (posterior.var * 
                            ((prior.mean/prior.var) + 
                               ((1/group.sigma.error) * 
                                  tapply(randomSlopes, data.currentCluster[, ids.var],sum)))) 
        # draw the random intercepts
        group.betas[group.clusterAssignments == h, 2] = 
          rnorm(perClusterN[h], posterior.mean, sqrt(posterior.var))

        # calculating the subject specific deviations from the cluster means
        # used to simplify calculation of the covariance
        #calculate betas minus their means
        group.betas.deviations[group.clusterAssignments == h,] = (
          group.betas[group.clusterAssignments == h,] - 
            (matrix(rep(means.currentCluster, perClusterN[h]), ncol=3, byrow=TRUE))
        )

      } # END CLUSTER-SPECIFIC UPDATE LOOP
      
      # update the model iteration
      model.current$betas[[group.index]] = group.betas
      model.current$betas.deviations[[group.index]] = group.betas.deviations
      
      # Common covariance for each cluster (Step 4)
      # Update Cluster Covariance
      model.current$dp.cluster.sigma[[group.index]] = riwish(
        prior.options$dp.cluster.sigma.nu0 + nsub,
        prior.options$dp.cluster.sigma.T0 + crossprod(group.betas.deviations)
      )

      ### Update the hyprparameters of the baseline distribution of 
      ### the Dirichlet process (Step 6)
      # update the mean
      sigma0.inv = solve(dp.dist.sigma0)
      Sb.inv = solve(prior.options$dp.dist.mu0.Sb)
      var = solve(Sb.inv + numNonEmptyClusters * sigma0.inv)
      if (numNonEmptyClusters == 1) {
        m <- var %*% (Sb.inv %*% prior.options$dp.dist.mu0.mb + 
                        sigma0.inv %*% ((model.current$clusterMeans[perClusterN > 0,])) )
      } else {
        m <- var %*% (prior.options$dp.dist.mu0.Sb %*% mb0 + 
                        sigma0.inv %*% (colSums(model.current$clusterMeans[[group.index]][perClusterN > 0,])) )
      }
      model.current$dp.dist.mu0[[group.index]] = rmvnorm(1,m,var)
      
      # update the variance of the baseline distribution
      model.current$dp.dist.sigma0[[group.index]] = riwish(
        prior.options$dp.dist.sigma0.nub + numNonEmptyClusters,
        prior.options$dp.dist.sigma0.Tb + 
          crossprod(model.current$clusterMeans[[group.index]][perClusterN > 0,] - 
                      matrix(rep(model.current$dp.dist.mu0[[group.index]],ncluster), 
                             numNonEmptyClusters, byrow=T))
      )
      
      if (dist = "gaussian") {
        # Step 9 Update sigma.error with inverse gamma
        residual = (group.data[, outcomes.var] - 
                      (Z[,1] * alpha[,1] + Z[,2] * alpha[,2]) -
                      ifelse(is.null(covariates), 
                             as.vector(as.matrix(group.data[,covariates.var]) %*% 
                                         model.current$betas.covariates[[group.index]]), 
                             0))
        g <- prior.options$sigma.error.tau1 + crossprod(residual)/2
        tau <- rgamma(1, prior.options$sigma.error.tau1 + group.N / 2, g)
        model.current$sigma.error[[group.index]] = 1 / tau
      } else {
        
      }

      #Estimate Density of Slope - for plotting density of the random slope (output only)  
      tmp5<-tmp6<-NULL
      for(h in 1:H){tmp5<-rbind(tmp5,pi[h]*dnorm(grid, mu[h,2], sqrt(Sigma.b[2,2])))
      tmp6<-rbind(tmp6,pi[h]*dnorm(grid, mu[h,1], sqrt(Sigma.b[1,1])))
      }
      dens1[i,]<-colSums(tmp5)
      dens0[i,]<-colSums(tmp6)
      
      # Equation 19
      #Estimate slope at each dropout time
      for (h in 1:H) tmp4[,h]<-pi[h]*dnorm(log(seq(1/15,16/15,1/15)),mu[h,3],sqrt(covar[9]))
      p.u<-tmp4/apply(tmp4,1,sum)
      for(uu in 1:16){uspecific.zz[i,uu]<-sum(p.u[uu,]*(mu[,2]+covar[6]*(log(uu/15)-mu[,3])/covar[9]))}
      
      #Estimate E(B1), E(B0) - section 2.4 expectation of random effects (int, slope)
      expected.zz1[i]<-sum(pi*mu[,2])
      expected.zz0[i]<-sum(pi*mu[,1])
      
    } # END GROUP-SPECIFIC UPDATE LOOP
      
    
    # update fixed effects associated with covariates
    if (dist == "gaussian") {
      sigma.inv = solve(model.current$beta.covariates.sigma)
      prior.sigma.inv = solve(prior.options$betas.covariates.sigma)
      var = solve(prior.sigma.inv + numSubjects * sigma.inv)
      m <- var %*% (prior.sigma.inv %*% prior.options$betas.covariates.mu + 
                      sigma.inv %*% (colSums()) )
      model.current$betas.covariates = rmvnorm(1,m,var)
      
    } else {
      # metropolis hastings
    }

    

  } # END ITERATION LOOP
}
    



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






