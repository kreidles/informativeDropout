#
# Functions for fitting the Dirichlet process approach to controlling
# for non-informative dropout
#
#
require(matrixcalc)

#' Data stored with each iteration of the Dirichlet process model
#'
#' @param weights.mixing vector containing the probability of belonging to cluster k, 
#'      for k = 1 to the number of clusters
#' @param weights.conditional vector containing the probability of belonging to cluster k, given
#' that the subject was not in clusters 1 to k - 1
#' @param cluster.assignments current cluster assignments
#' @param betas A (k x 3) matrix of regression coefficients for the random intercept, slope,
#' and log of the dropout time for each cluster, with k = number of clusters
#' @param betas.covariates A (c x 1) vector of regression coefficients for covariates, 
#' with c = number of covariates
#' @param betas.covariates.mu A (c x 1) vector representing the mean of the distribution of 
#' regression coefficients related to covariates
#' @param betas.covariates.sigma A (c x c) vector representing the covariance of the distribution of 
#' regression coefficients related to covariates
#' @param betas.deviations An (N x k) matrix of subject specific deviations from the cluster means
#' @param dp.dist.mu0 A (3 x 1) vector of means for the baseline distribution of the Dirichlet process
#' @param dp.dist.sigma0 A (3 x 3) matrix representing the covariance of the baseline 
#' distribution of the Dirichlet process
#' @param dp.cluster.sigma A (3 x 3) matrix representing the covariance of the random intercept, slope,
#' and log dropout time for each cluster.  This covariance is the same for each cluster
#' @param dp.concentration the concentration parameter (alpha) for the Dirichlet process
#' @param cluster.N A (k x 1) vector indicating the number of subjects in cluster k, for k = 1 to the
#' number of clusters
#' @param cluster.mu A (k x 3) matrix of means for the random intercept, slope, and log dropout time
#' in cluster k, for k = 1 to the number of clusters
#' @param expected.intercept the expected value of the random intercept
#' @param expected.slope the expected value of the random slope
#' @param sigma.error For Gaussian outcomes, the residual error variance
#' 
#' @exportClass dirichlet.iteration
#' 
dirichlet.iteration <- function(weights.mixing=NULL, weights.conditional=NULL,
                                cluster.assignments = NULL,
                                betas = NULL, 
                                betas.deviations = NULL,
                                betas.covariates = NULL, 
                                betas.covariates.mu = NULL,
                                betas.covariates.sigma = NULL,
                                dp.dist.mu0 = NULL, dp.dist.sigma0 = NULL,
                                dp.cluster.sigma = NULL,
                                dp.concentration=NULL, 
                                cluster.N = NULL,
                                cluster.mu=NULL, 
                                sigma.error = NULL,
                                expected.intercept=NULL, expected.slope=NULL,
                                slope.dropoutTimes=NULL,
                                density.intercept = NULL,
                                density.slope = NULL) {
  mi = list(
    # mixing weights
    weights.mixing = weights.mixing,
    # posterior probabilities that subject i belongs to cluster h
    weights.conditional = weights.conditional,
    # current cluster assignments by group
    cluster.assignments = cluster.assignments,
    # regression coefficients for random intercept and slope,
    # concatenated with the log of the dropout time for convenience
    betas = betas,
    # regression coefficients associated with the shared covariates 
    # (includes group offsets)
    betas.covariates = betas.covariates,
    # mean of distribution of covariate regression coefficients
    betas.covariates.mu = betas.covariates.mu,
    # covariance of distribution of covariate regression coefficients
    betas.covariates.sigma = betas.covariates.sigma,
    # subject specific deviations from the cluster mean
    betas.deviations = betas.deviations,
    # mean and variance of the baseline distribution of the Dirichlet process
    dp.dist.mu0 = dp.dist.mu0,
    dp.dist.sigma0 = dp.dist.sigma0,
    # common covariance for cluster random effects
    dp.cluster.sigma = dp.cluster.sigma,
    # concentration parameter
    dp.concentration = dp.concentration,
    # number of subjects per cluster
    cluster.N = cluster.N,
    # cluster specific means
    cluster.mu = cluster.mu,
    # residual model error
    sigma.error = sigma.error,
    
    ## summary statistics calculated at each iteration
    # expected values of the random intercept and slope
    expected.intercept = expected.intercept,
    expected.slope = expected.slope,
    # estimated slopes at specified dropout times
    slope.dropoutTimes = slope.dropoutTimes,
    # densities of the random effects
    density.intercept = density.intercept,
    density.slope = density.slope

  )
  
  class(mi) <- append(class(mi), "dirichlet.iteration")
  return(mi)
  
}

#' Model options for the Dirichlet Process model.  Includes starting values, priors and
#' simulation/MCMC parameters.
#'
#' @param iterations number of iterations for the MCMC simulation
#' @param burnin burn in period for the simulation, i.e. the number of 
#' iterations to throw away at the beginning of the simulation
#' @param thin thinning interval, i.e. if thin=n, only keep every nth iteration
#' @param numClusters number of clusters for the Dirichlet Process stick breaking model
#' @param dp.concentration prior value for the concentration parameter of the Dirichlet process
#' @param dp.concentration.alpha Shape parameter for the hyperprior Gamma distribution of the
#' concentration parameter of the Dirichlet process
#' @param dp.concentration.beta Rate parameter for the hyperprior Gamma distribution of the
#' concentration parameter of the Dirichlet process
#' @param dp.cluster.sigma prior for the common cluster covariance
#' @param dp.cluster.sigma.nu0 Degrees of freedom for the hyperprior inverse Wishart distribution of the
#' common cluster covariance
#' @param dp.cluster.sigma.T0 Scale matrix for the hyperprior inverse Wishart distribution of the
#' common cluster covariance
#' @param dp.dist.mu0 prior mean for the baseline distribution of the Dirichlet process
#' @param dp.dist.mu0.mb Mean for the hyperprior Gaussian distribution of the mean of
#' the baseline distribution of the Dirichlet process
#' @param dp.dist.mu0.Sb Covariance for the hyperprior Gaussian distribution of the mean of
#' the baseline distribution of the Dirichlet process
#' @param dp.dist.sigma0 prior covariance for the baseline distribution of the Dirichlet process
#' @param dp.dist.sigma0.nub Degrees of freedom for the hyperprior inverse Wishart distribution of 
#' the covariance of the baseline distribution of the Dirichlet process
#' @param dp.dist.sigma0.Tb Scale matrix for the hyperprior inverse Wishart distribution of 
#' the covariance of the baseline distribution of the Dirichlet process
#' @param betas.covariates prior value for covariate regression coefficients
#' @param betas.covariates.mu Prior mean for the covariate regression coefficients
#' @param betas.covariates.sigma Prior covariance for the covariate regression coefficients
#' @param sigma.error Prior for the residual error (Gaussian outcomes only)
#' 
#' @exportClass dirichlet.model.options
#' 
dirichlet.model.options <- function(iterations=10000, burnin=500, thin=NULL,
                                    start.weights.mixing=NULL, n.clusters=15,
                                    dropout.estimationTimes=NULL, dropout.offset=1/15,
                                    dp.concentration=NULL,
                                    dp.concentration.alpha=NULL,
                                    dp.concentration.beta=NULL,
                                    dp.cluster.sigma = NULL,
                                    dp.cluster.sigma.nu0 = NULL,
                                    dp.cluster.sigma.T0 = NULL,
                                    dp.dist.mu0 = NULL,
                                    dp.dist.mu0.mb = NULL,
                                    dp.dist.mu0.Sb = NULL,
                                    dp.dist.sigma0 = NULL,
                                    dp.dist.sigma0.nub = NULL,
                                    dp.dist.sigma0.Tb = NULL,
                                    betas.covariates = NULL,
                                    betas.covariates.mu = NULL,
                                    betas.covariates.sigma = NULL,
                                    sigma.error = NULL,
                                    sigma.error.tau = NULL) {
  
  # validate the mcmc options
  if (is.na(iterations) || iterations <= 1) {
    stop("MCMC options error :: invalid number of iterations")
  }
  if (is.na(burnin) || burnin >= iterations) {
    stop("MCMC options error :: the burn in period must be less than the total number of iterations")
  }
  if (is.na(n.clusters) || n.clusters <= 1) {
    stop("Start options error :: invalid number of clusters")
  }
  if (!is.null(start.weights.mixing) && length(start.weights.mixing) != n.clusters) {
    stop("Start options error :: length of mixing weight start values must equal the number of clusters")
  }
  
  # validate the prior options
  # Dirichlet process details
  if (is.null(dp.concentration) || dp.concentration <= 0) {
    stop("Prior options error :: invalid concentration parameter for the Dirichlet process")
  }
  if (is.null(dp.concentration.alpha) || dp.concentration.alpha <= 0) {
    stop("Prior options error :: invalid gamma prior for the Dirichlet process concentration parameter")
  }
  if (is.null(dp.concentration.beta) || dp.concentration.beta <= 0) {
    stop("Prior options error :: invalid gamma prior for the Dirichlet process concentration parameter")
  }
  # priors for baseline distribution of the Dirichlet process
  # mean of the baseline distribution
  if (is.null(dp.dist.mu0)) {
    stop("Prior options error :: invalid mean for the Dirichlet process baseline distribution")
  }
  # hyperparams for the mean of the baseline distribution
  if (is.null(dp.dist.mu0.mb)) {
    stop("Prior options error :: invalid mean (hyperparameter) for the mean of the Dirichlet process baseline distribution")
  }
  if (is.null(dp.dist.mu0.Sb) || 
      !is.matrix(dp.dist.mu0.Sb) ||
      nrow(dp.dist.mu0.Sb) != 3 || ncol(dp.dist.mu0.Sb) != 3 ||
      !is.positive.definite(dp.dist.mu0.Sb)) {
    stop("Prior options error :: invalid variance (hyperparameter) for the mean of the Dirichlet process baseline distribution")
  }
  
  # variance of the baseline distribution
  if (is.null(dp.dist.sigma0)) {
    stop("Prior options error :: no variance for the Dirichlet process baseline distribution")
  }
  # hyperparameters of the variance of the baseline distribution
  if (is.null(dp.dist.sigma0.nub)) {
    stop("Prior options error :: invalid df (hyperparameter) for the variance of the Dirichlet process baseline distribution")
  }
  if (is.null(dp.dist.sigma0.Tb) || 
      !is.matrix(dp.dist.sigma0.Tb) ||
      nrow(dp.dist.sigma0.Tb) != 3 || ncol(dp.dist.sigma0.Tb) != 3 ||
      !is.positive.definite(dp.dist.sigma0.Tb)) {
    stop("Prior options error :: invalid scale matrix (hyperparameter) for the variance of the Dirichlet process baseline distribution")
  }
  
  # cluster specific variance for the Dirichlet process
  if (is.null(dp.cluster.sigma) || 
      !is.matrix(dp.cluster.sigma) ||
      nrow(dp.cluster.sigma) != 3 || ncol(dp.cluster.sigma) != 3 ||
      !is.positive.definite(dp.cluster.sigma)) {
    stop("Prior options error :: invalid cluster-specific variance for the Dirichlet process")
  }
  # hyperparameters of the cluster specific variance
  if (is.null(dp.cluster.sigma.nu0)) {
    stop("Prior options error :: invalid df (hyperparameter) for the cluster-specific variance of the Dirichlet process")
  }
  if (is.null(dp.cluster.sigma.T0) || 
      !is.matrix(dp.cluster.sigma.T0) ||
      nrow(dp.cluster.sigma.T0) != 3 || ncol(dp.cluster.sigma.T0) != 3 ||
      !is.positive.definite(dp.cluster.sigma.T0)) {
    stop("Prior options error :: invalid scale matrix (hyperparameter) for the cluster-specific variance of the Dirichlet process")
  }
  
  opts = list(
    # mcmc iterations
    iterations = iterations,
    # burn in period
    burnin = burnin,
    # thinning interval
    thin = thin,
    # starting value for mixing weights
    start.weights.mixing = start.weights.mixing,
    # number of clusters in the Dirchlet Process
    n.clusters = n.clusters,
    # times at which the dropout time dependent slopes will be estimated
    dropout.estimationTimes = dropout.estimationTimes,
    # offset to avoid dropout times of 0
    dropout.offset = dropout.offset,
    # prior concentration of the DP
    dp.concentration = dp.concentration,
    # hyperprior for Gamma distribution of the DP concentration
    dp.concentration.alpha = dp.concentration.alpha,
    dp.concentration.beta = dp.concentration.beta,
    
    # prior value for the common cluster covariance
    dp.cluster.sigma = dp.cluster.sigma,
    # hyperprior for the common cluster covariance
    dp.cluster.sigma.nu0 = dp.cluster.sigma.nu0,
    dp.cluster.sigma.T0 = dp.cluster.sigma.T0,
    
    # prior value for the mean of the baseline distribution of the DP
    dp.dist.mu0 = dp.dist.mu0,
    # hyperprior for the mean of the baseline distribution
    dp.dist.mu0.mb = dp.dist.mu0.mb,
    dp.dist.mu0.Sb = dp.dist.mu0.Sb,
    
    # prior value for the covariance of the baseline distribution of the DP
    dp.dist.sigma0 = dp.dist.sigma0,
    # hyperprior for the covariance of the baseline distribution
    dp.dist.sigma0.nub = dp.dist.sigma0.nub,
    dp.dist.sigma0.Tb = dp.dist.sigma0.Tb,
    
    # prior for the regression coefficients related to covariates
    betas.covariates = betas.covariates,
    # hyperprior for the covariate regression coefficients
    betas.covariates.mu = betas.covariates.mu,
    betas.covariates.sigma = betas.covariates.sigma,
    
    # residual error (Gaussian outcomes only)
    sigma.error = sigma.error,
    sigma.error.tau = sigma.error.tau
    
  )
  
  class(opts) <- append(class(opts), "dirichlet.model.options")
  return(opts)
  
}







#'
#' infortmativeDropout.bayes.dirichlet
#'
#' Fit a bayesian model which accounts for dropout by modeling the 
#' releationship between dropout time and slope with natural cubic B-splines
#'
#' @param data the data set 
#' @param ids.var column of the data set containing participant identifiers
#' @param outcomes.var column of the data set containing the outcome variable
#' @param groups.var column of the data set indicating the treatment groups
#' @param covariates.var list of columns in the data set containing covariates
#' @param times.dropout.var column of the data set containing dropout times
#' @param times.observation.var column of the data set containing 
#' @param dist the distribution of the outcome, valid values are "gaussian" or "binary"
#' @param model.options model options (see dirichlet.model.options for details)
#'
#'
informativeDropout.bayes.dirichlet <- function(data, ids.var, outcomes.var, groups.var, 
                                               covariates.var, 
                                               times.dropout.var, times.observation.var, 
                                               dist, model.options) {
  
  # make sure we have the correct type of mcmc opts
  if (!is(model.options,"dirichlet.model.options")) {
    stop("Model options error :: options must be of type dirichlet.model.options")
  }

  # prior for covariate effects
  if (!is.null(covariates.var)) {
    if (is.na(model.options$betas.covariates.mu) || 
        length(model.options$betas.covariates.mu) != length(covariates.var)) {
      stop("Prior options error :: invalid prior mean for fixed effects related to covariates")
    }
    if (is.na(model.options$betas.covariates.R0) || 
        !is.matrix(model.options$betas.covariates.R0) ||
        nrow(model.options$betas.covariates.R0) != length(covariates.var) || 
        ncol(model.options$betas.covariates.R0) != length(covariates.var) ||
        !is.positive.definite(model.options$betas.covariates.R0)) {
      stop("Prior options error :: invalid prior variance for fixed effects related to covariates")
    }
  }
  
  if (dist == 'gaussian') {
    if (is.null(model.options$sigma.error.tau)) {
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
    return(model.options$dp.concentration)
  })
  
  # initialize cluster specific means
  cluster.mu = lapply(groupList, function(group) {
    means = matrix(0,n.clusters,3)
    return(means)
  })
  # initialize the cluster specific variance
  cluster.sigma = lapply(groupList, function(group) {
    covar = model.options$dp.cluster.sigma
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
    if (!is.null(model.options$betas.covariates)) {
      start.betas.covariates = model.options$betas.covariates
    } else {
      start.betas.covariates = getInitialEstimatesCovariates(dist, 
                                                             as.data.frame(data[,covariates.var]), 
                                                             data[,outcomes.var])
    }
  }
  
  # starting values for the mean and variance of the Dirichlet baseline distribution
  dp.dist.mu0 = lapply(groupList, function(group) {
    return(model.options$dp.dist.mu0)
  })
  dp.dist.sigma0 = lapply(groupList, function(group) {
    return(model.options$dp.dist.sigma0)
  })
  # starting values for the common cluster variance
  dp.cluster.sigma = lapply(groupList, function(group) {
    return(model.options$dp.cluster.sigma)
  })
  
  # for Gaussian outcomes, the residual error
  if (dist == "gaussian") {
    sigma.error = model.options$sigma.error
  }
  
  
  # initialize the first model iteration
  modelIterationList <- vector(mode = "list", length = model.options$iterations)
  modelIterationList[[1]] = dirichlet.iteration(
    weights.mixing=weights.mixing, 
    weights.conditional=weights.conditional,
    betas = betas, 
    betas.deviations = betas.deviations,
    betas.covariates = start.betas.covariates, 
    betas.covariates.mu = model.options$betas.covariates.mu,
    betas.covariates.sigma = model.options$betas.covariates.sigma,
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
                       model.options$dp.concentration + 
                         sum(group.cluster.N[(h+1):n.clusters])))
        } else {
          return(1)
        }
      })
      
      # number of clusters that have at least one subject (common to have empty clusters)
      numNonEmptyClusters = length(group.cluster.N[group.cluster.N != 0])
      
      # Alpha is concentration parameter of the dirichlet process (step 7)
      # sample latent beta distributed variable (Escobar and West 1995)
      eta <- rbeta(1, model.options$dp.concentration + 1, group.n)
      # sample the concentration parameter from a mixture of Gamma distributions
      pi.eta <- ((model.options$dp.concentration.alpha + numNonEmptyClusters - 1) /
                   (group.n * (model.options$dp.concentration.beta - log(eta)) +
                      model.options$dp.concentration.alpha + numNonEmptyClusters - 1))
      # now draw the new concentration parameter for the Dirichlet process
      model.current$dp.concentration[[group.index]] <- 
        ((pi.eta * rgamma(1, model.options$dp.concentration.alpha + numNonEmptyClusters, 
                          model.options$dp.concentration.beta - log(eta))) + 
           ((1 - pi.eta) * rgamma(1, model.options$dp.concentration.alpha + numNonEmptyClusters - 1, 
                                  model.options$dp.concentration.beta - log(eta))))
      
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
        
        # get the data for subjects in this cluster
        tempid = rep(group.cluster.assignments, group.nobs)
        data.currentCluster = group.data[tempid==h,]
        
        if (dist == 'gaussian') {
          ## calculate the posterior mean and variance for the random intercept
          posterior.var = 1 / ((1 / prior.var) + (group.nobs[group.cluster.assignments==h] / model.current$sigma.error))
          
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
        } else {
          # binary distribution
          # Current Observation-Level Log-Likelihood
          eta<-rep(betas[c==h,1], nobs[c==h])+t[tempid==h]*rep(betas[c==h,2],nobs[c==h])
          Lt<-Lstar<-rep(0, length(tempid[tempid==h]))
          Lt[y[tempid==h]==0]<-log(1-inv.logit(eta[y[tempid==h]==0]))  
          Lt[y[tempid==h]==1]<-log(inv.logit(eta[y[tempid==h]==1]))
          # Draw Candidate from Symmetric Rand Walk Proposal 
          bs<-betas[c==h,1]+rnorm(ns[h])      # or b+rmvnorm(n,sigma=Sigma)
          etastar<-rep(bs,nobs[c==h])+t[tempid==h]*rep(betas[c==h,2],nobs[c==h])
          pstar<-inv.logit(etastar)
          Lstar[y[tempid==h]==0]<-log(1-pstar[y[tempid==h]==0])
          Lstar[y[tempid==h]==1]<-log(pstar[y[tempid==h]==1])
          # Subject-Level Log-Likelihood
          l<-tapply(Lt,patid[tempid==h],sum)
          ls<-tapply(Lstar, patid[tempid==h], sum)
          
          ratio<-ls+dnorm(bs,pm,sqrt(pv),log=TRUE)-(l+dnorm(betas[c==h,1],pm,sqrt(pv),log=TRUE))
          
          un<-1*(log(runif(ns[h]))<ratio)
          betas[c==h,1][un==1]<-bs[un==1]
          
        }
        
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
        model.options$dp.cluster.sigma.nu0 + group.n, # <- check typo
        model.options$dp.cluster.sigma.T0 + crossprod(group.betas.deviations)
      )
      
      ### Update the hyprparameters of the baseline distribution of 
      ### the Dirichlet process (Step 6)
      # update the mean
      Sb.inv = solve(model.options$dp.dist.mu0.Sb)
      var = solve(Sb.inv + numNonEmptyClusters * dp.dist.sigma0.inv)
      if (numNonEmptyClusters == 1) {
        m <- var %*% (Sb.inv %*% model.options$dp.dist.mu0.mb + 
                        dp.dist.sigma0.inv %*% ((model.current$cluster.mu[[group.index]][group.cluster.N > 0,])) )
      } else {
        m <- var %*% (model.options$dp.dist.mu0.Sb %*% model.options$dp.dist.mu0.mb + 
                        dp.dist.sigma0.inv %*% (colSums(model.current$cluster.mu[[group.index]][group.cluster.N > 0,])) )
      }
      model.current$dp.dist.mu0[[group.index]] = as.vector(rmvnorm(1,m,var))
      
      # update the variance of the baseline distribution
      model.current$dp.dist.sigma0[[group.index]] = riwish(
        model.options$dp.dist.sigma0.nub + numNonEmptyClusters,
        model.options$dp.dist.sigma0.Tb + 
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
      g <- model.options$sigma.error.tau + crossprod(residual)/2
      tau <- rgamma(1, model.options$sigma.error.tau + n.total / 2, g)
      model.current$sigma.error = 1 / tau
    } 
    
    # update fixed effects associated with covariates
    if (!is.null(covariates.var)) {
      if (dist == "gaussian") {
        sigma.error.inv = 1/model.current$sigma.error
        prior.sigma.inv = solve(model.options$betas.covariates.sigma)
        
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





