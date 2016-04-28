#
# Example simulations to accompany the manuscript
#
#

#
# Dirichlet process model examples corresponding to the vignettes
#

gendata.gaussian_1group_nocovar <- function() {
  
}






example.dp_gaussian_1group_nocovar <- function() {
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
                                        sigma.error.tau = 0.01)
  
  
  set.seed(1066)
  result = informativeDropout.bayes.dirichlet(data, ids.var, 
                                              outcomes.var, groups.var,
                                              covariates.var, 
                                              times.dropout.var, times.observation.var,
                                              dist, model.options)
  #   result = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
  #                               times.dropout.var, times.observation.var, 
  #                               method, dist,
  #                               knots.options = knots.options, 
  #                               mcmc.options = mcmc.options,
  #                               model.options = model.options,
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



example.dp_binary_1group_nocovar <- function() {
  data <- read.table("../Rnsv/code/sim_sml_1.dat")
  #data$day = data$years * 365
  
  names(data) <- c("patid", "alpha", "drptm", "b1", "b2",
                   "b2ui", "b2uii", "b2uiii", "t", "e", "yi", "yii", "yiii")
  data$group <- rep(1,nrow(data))
  data$yi_bin = (data$yi> 0)
  # for debugging
  ids.var = "patid"
  outcomes.var = "yi_bin"
  groups.var = "group"
  covariates.var = NULL
  times.dropout.var = "drptm"
  times.observation.var = "t"
  method="bayes.dirichlet"
  dist = "binary"
  
  model.options=dirichlet.model.options(iterations=10, n.clusters=15, burnin=0,
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
                                        betas.covariates.sigma = NULL)
  
  
  set.seed(1066)
  result = informativeDropout.bayes.dirichlet(data, ids.var, 
                                              outcomes.var, groups.var,
                                              covariates.var, 
                                              times.dropout.var, times.observation.var,
                                              dist, model.options)
  #   result = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
  #                               times.dropout.var, times.observation.var, 
  #                               method, dist,
  #                               knots.options = knots.options, 
  #                               mcmc.options = mcmc.options,
  #                               model.options = model.options,
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




runExamples <- function() {
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
                                        sigma.error.tau = 0.01)
  
  
  set.seed(1066)
  result = informativeDropout.bayes.dirichlet(data, ids.var, 
                                              outcomes.var, groups.var,
                                              covariates.var, 
                                              times.dropout.var, times.observation.var,
                                              dist, model.options)
  #   result = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
  #                               times.dropout.var, times.observation.var, 
  #                               method, dist,
  #                               knots.options = knots.options, 
  #                               mcmc.options = mcmc.options,
  #                               model.options = model.options,
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











