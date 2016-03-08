#
# Functions for fitting the Dirichlet process approach to controlling
# for non-informative dropout
#
#

#' Data stored with each iteration of the Dirichlet process model
#'
#' @param weights.mixing vector containing the probability of belonging to cluster k, 
#'      for k = 1 to the number of clusters
#' @param weights.conditional vector containing the probability of belonging to cluster k, given
#' that the subject was not in clusters 1 to k - 1
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
                                betas = NULL, 
                                betas.covariates = NULL, 
                                betas.covariates.mu = NULL,
                                betas.covariates.sigma = NULL,
                                betas.deviations = NULL,
                                dp.dist.mu0 = NULL, dp.dist.sigma0 = NULL,
                                dp.cluster.sigma = NULL,
                                dp.concentration=NULL, 
                                cluster.N = NULL,
                                cluster.mu=NULL, 
                                expected.intercept=NULL, expected.slope=NULL,
                                sigma.error = NULL) {
  mi = list(
    # mixing weights
    weights.mixing = weights.mixing,
    # posterior probabilities that subject i belongs to cluster h
    weights.conditional = weights.conditional,
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
    # expected values of the random intercept and slope
    expected.intercept = expected.intercept,
    expected.slope = expected.slope,
    # residual model error
    sigma.error = sigma.error
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
                                    sigma.error = NULL) {
  mi = list(
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
    sigma.error = sigma.error
  )
  
  class(mi) <- append(class(mi), "dirichlet.prior.options")
  return(mi)
  
}

#' Starting values for the MCMC loop in the Dirichlet Process model
#' 
#' @param weights.mixing vector containing starting values for the probability of 
#' belonging to cluster k, for k = 1 to the number of clusters
#' 
#' @export
#' 
dirichlet.start.options = function(weights.mixing=NULL) {
  mi = list(
    # mixing weights
    weights.mixing = weights.mixing
  )
  
  class(mi) <- append(class(mi), "dirichlet.start.options")
  return(mi)
}



test.sim <- function() {
  data <- read.table("../Rnsv/code/sim_sml_1.dat")
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
  method="bayes.dirichlet"
  dist = "gaussian"
  knots.options=list(birthProbability=0.2, min=1, max=10, stepSize=0.1,
                     startPositions=c(0,7/30,0.5, 23/30,1), candidatePositions=seq(0,1,0.1/3)) 
  mcmc.options=list(iterations=40000, burnIn=10, numClusters=15)
  prior.options=list(dp.concentration=1,
                     dp.concentration.alpha=1,
                     dp.concentration.beta=1,
                     dp.cluster.sigma.nu0 = 5,
                     dp.cluster.sigma.T0 = diag(3),
                     betas.covariates.mu = NULL,
                     betas.covariates.sigma = NULL)
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






