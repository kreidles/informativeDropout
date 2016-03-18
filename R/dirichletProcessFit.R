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
#' @export
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

#' Priors for all components of the Dirichlet process model
#'
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
#' @export
#' 
dirichlet.prior.options <- function(dp.concentration=NULL,
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
  
  class(opts) <- append(class(opts), "dirichlet.prior.options")
  return(opts)
  
}

#'
#' Simulation and model options for the Dirichlet Process model
#' 
#' @param iterations number of iterations for the MCMC simulation
#' @param burnin burn in period for the simulation, i.e. the number of 
#' iterations to throw away at the beginning of the simulation
#' @param thin thinning interval, i.e. if thin=n, only keep every nth iteration
#' @param numClusters number of clusters for the Dirichlet Process stick breaking model
#' 
#' @export
#' 
dirichlet.model.options = function(iterations=10000, burnin=500, thin=NULL,
                             start.weights.mixing=NULL, n.clusters=15,
                             dropout.estimationTimes=NULL, dropout.offset=1/15) {
  # TODO: add na.rm to ignore subjects with missing outcomes or covariates
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
    dropout.offset = dropout.offset
  )
  
  class(opts) <- append(class(opts), "dirichlet.model.options")
  return(opts)
  
}





test.sim <- function() {
  data <- read.table("../Rnsv/code/sim_sml_1.dat")
  #data$day = data$years * 365
  
  names(data) <- c("patid", "alpha", "drptm", "b1", "b2",
                   "b2ui", "b2uii", "b2uiii", "t", "e", "yi", "yii", "yiii")
  data$group <- rep(1,nrow(data))
  # for debugging
  ids.var = "patid"
  outcomes.var = "yi"
  groups.var = "group"
  covariates.var = NULL
  times.dropout.var = "drptm"
  times.observation.var = "t"
  method="bayes.dirichlet"
  dist = "gaussian"

  model.options=dirichlet.model.options(iterations=5000, n.clusters=15, burnin=0,
                                       dropout.estimationTimes = seq(1/15,1,1/15))
  
  prior.options = dirichlet.prior.options(dp.concentration=1,
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
                                          sigma.error.tau = 0.01)
    
  
  set.seed(1066)
  result = informativeDropout.bayes.dirichlet(data, ids.var, 
                                              outcomes.var, groups.var,
                                              covariates.var, 
                                              times.dropout.var, times.observation.var,
                                              dist, model.options, prior.options)
#   result = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
#                               times.dropout.var, times.observation.var, 
#                               method, dist,
#                               knots.options = knots.options, 
#                               mcmc.options = mcmc.options,
#                               prior.options = prior.options,
#                               dropoutEstimationTimes = dropoutEstimationTimes)
#   
#   
#   acceptanceProbability(result, "knot.add")
#   acceptanceProbability(result, "knot.remove")
#   acceptanceProbability(result, "knot.move")
#   acceptanceProbability(result, "fixedEffects")
#   #acceptanceProbability(result, "fixedEffectsCovariates")
#   
#   
#   nknots = unlist(lapply(result, function(x) { return(length(x$knots[[1]])) } ))
#   ts.plot(nknots)
#   summary(nknots)
#   
#   slopes = unlist(lapply(result, function(x) { return(x$slope.marginal[[1]]) }))
#   summary(slopes)
#   ts.plot(slopes)
#   
#   sum.sigma.error = unlist(lapply(result, function(x) { return(x$sigma.error) }))
#   summary(sum.sigma.error)
#   ts.plot(sum.sigma.error)
#   
#   sum.sigma.randomIntercept = unlist(lapply(result, function(x) { return(x$sigma.randomIntercept) }))
#   summary(sum.sigma.randomIntercept)
#   ts.plot(sum.sigma.randomIntercept)
#   
#   sum.sigma.randomSlope = unlist(lapply(result, function(x) { return(x$sigma.randomSlope) }))
#   summary(sum.sigma.randomSlope)
#   ts.plot(sum.sigma.randomSlope)
#   
#   sum.sigma.randomInterceptSlope = unlist(lapply(result, function(x) { return(x$sigma.randomInterceptSlope) }))
#   summary(sum.sigma.randomInterceptSlope)
#   ts.plot(sum.sigma.randomInterceptSlope)
#   
#   
#   dropout.slopes1 = unlist(lapply(result, function(x) { return(x$slope.dropoutSpecific[[1]][1])}))
#   summary(dropout.slopes1)
#   ts.plot(dropout.slopes1)
#   
#   dropout.slopes2 = unlist(lapply(result, function(x) { return(x$slope.dropoutSpecific[[1]][2])}))
#   summary(dropout.slopes2)
#   ts.plot(dropout.slopes2)
#   
#   dropout.slopes3 = unlist(lapply(result, function(x) { return(x$slope.dropoutSpecific[[1]][3])}))
#   summary(dropout.slopes3)
#   ts.plot(dropout.slopes3)
}






