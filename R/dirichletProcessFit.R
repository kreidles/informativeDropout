#
# Functions for fitting the Dirichlet process approach to controlling
# for non-informative dropout
#
#

dirichlet.iteration <- function(mixingWeights=NULL, conditionalWeights=NULL,
                                betas = NULL, betas.covariates = NULL, 
                                betas.deviations = NULL,
                                perClusterN = NULL,
                                dp.dist.mu0 = NULL, dp.dist.sigma0 = NULL,
                                dp.cluster.sigma = NULL,
                                clusterMeans=NULL, clusterCovariance=NULL, 
                                expected.intercept=NULL, expected.slope=NULL,
                                concentration=NULL, 
                                sigma.error = 1) {
  mi = list(
    # mixing weights
    mixingWeights = mixingWeights,
    # posterior probabilities that subject i belongs to cluster h
    conditionalWeights = conditionalWeights,
    # regression coefficients for random intercept and slope,
    # concatenated with the log of the dropout time for convenience
    betas = betas,
    # deviations from the cluster mean
    betas.deviations = betas.deviations,
    # regression coefficients associated with the shared covariates 
    # (includes group offsets)
    betas.covariates = betas.covariates,
    # number of subjects per cluster
    perClusterN = perClusterN,
    # mean and variance of the baseline distribution of the Dirichlet process
    dp.dist.mu0 = dp.dist.mu0,
    dp.dist.sigma0 = dp.dist.sigma0,
    # common covariance for cluster random effects
    dp.cluster.sigma = dp.cluster.sigma,
    # cluster specific means
    clusterMeans = clusterMeans,
    # cluster specific covariance
    clusterCovariance = clusterCovariance,
    # expected values of the random intercept and slope
    expected.intercept = expected.intercept,
    expected.slope = expected.slope,
    # concentration parameter
    concentration = concentration,
    # residual model error
    sigma.error = sigma.error
  )
  
  class(mi) <- append(class(mi), "dirichlet.iteration")
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
    shape.tau = 0.001, rate.tau = 0.001, lambda.numKnots = 5,
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






