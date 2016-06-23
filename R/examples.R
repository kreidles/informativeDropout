#
# Example simulations to accompany the manuscript
#
#' @include informativeDropout.R





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
  
  model.options <- bayes.splines.model.options(iterations=100, burnin=10, thin=NA,
                                               knots.prob.birth=0.5, knots.min=3, knots.max=10, knots.setpSize=3,
                                               knots.positions.start=c(330,550,1060), 
                                               knots.positions.candidate=seq(10,max(data$day),10),
                                               dropout.estimationTimes=c(1,2,3),
                                               shape.tau=0.001, rate.tau=0.001,
                                               sigma.beta=1, lambda.numKnots=1,
                                               sigma.error=1,
                                               sigma.error.df = 3,
                                               sigma.error.scale = diag(2),
                                               eta.null=NULL)
  
  
  set.seed(1066)
  result = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
                              times.dropout.var, times.observation.var, 
                              method, dist,
                              knots.options = knots.options, 
                              mcmc.options = mcmc.options,
                              prior.options = prior.options)
}

gaussian_1group_nocovar <- function() {
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
  
  model.options <- bayes.splines.model.options(iterations=100, burnin=10, thin=1,
                                               knots.prob.birth=0.2, knots.min=1, knots.max=10, 
                                               knots.stepSize=0.1,
                                               knots.positions.start=c(0,7/30,0.5, 23/30,1), 
                                               knots.positions.candidate=seq(0,1,0.1/3),
                                               dropout.estimationTimes=seq(0,1,1/15),
                                               sigma.error.shape.tau=0.001, sigma.error.rate.tau=0.001,
                                               sigma.beta=25, lambda.numKnots=5,
                                               sigma.residual = 1,
                                               sigma.error=1,
                                               sigma.randomIntercept = 1,
                                               sigma.randomSlope = 1,
                                               sigma.randomInterceptSlope = 0.001,
                                               sigma.randomEffects.df = 3,
                                               sigma.randomEffects.scale = diag(2),
                                               eta.null=NULL)
  model.options$sigma.randomIntercept
  
  
  set.seed(1066) 
  fit = informativeDropout(data, ids.var, outcomes.var, groups.var, covariates.var, 
                           times.dropout.var, times.observation.var, 
                           method, dist, model.options)
  
  summary(fit)
  
  
  
  nknots = unlist(lapply(fit$iterations, function(x) { return(length(x$knots[[1]])) } ))
  ts.plot(nknots)
  summary(nknots)
  
  slopes = unlist(lapply(fit$iterations, function(x) { return(x$slope.marginal[[1]]) }))
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

binary_1group_nocovar <- function() {
  data <- read.table("../../Rnsv/code/sim_sml_1.dat")
  #data$day = data$years * 365
  
  names(data) <- c("patid", "alpha", "drptm", "b1", "b2",
                   "b2ui", "b2uii", "b2uiii", "t", "e", "yi", "yii", "yiii")
  data$group <- rep(1,nrow(data))
  data$y_bin <- as.numeric(data$yiii > 0)
  data$drop <- data$drptm + 1/15
  # for debugging
  ids.var = "patid"
  outcomes.var = "y_bin"
  groups.var = "group"
  covariates.var = NULL
  times.dropout.var = "drop"
  times.observation.var = "t"
  method="bayes.splines"
  dist = "binary"
  
  model.options <- bayes.splines.model.options(iterations=50, burnin=0, thin=1, print=1,
                                               knots.prob.birth=0.2, knots.min=1, knots.max=10, 
                                               knots.stepSize=0.1,
                                               knots.positions.start=c(0,7/30,0.5, 23/30,1), 
                                               knots.positions.candidate=seq(0,1,0.1/3),
                                               dropout.estimationTimes=seq(0,1,1/15),
                                               sigma.error.shape.tau=0.001, sigma.error.rate.tau=0.001,
                                               sigma.beta=25, lambda.numKnots=5,
                                               sigma.residual = 1,
                                               sigma.error=1,
                                               sigma.randomIntercept = 1,
                                               sigma.randomSlope = 1,
                                               sigma.randomInterceptSlope = 0.001,
                                               sigma.randomEffects.df = 3,
                                               sigma.randomEffects.scale = diag(2),
                                               eta.null=0.2)
  model.options$sigma.randomIntercept
  
  
  set.seed(1066) 
  fit = informativeDropout.bayes.splines(data, ids.var, outcomes.var, groups.var, covariates.var, 
                           times.dropout.var, times.observation.var, 
                           dist, model.options)
  summary(fit)
  
  
  prob.acceptance(result, "knot.add")
  prob.acceptance(result, "knot.remove")
  prob.acceptance(result, "knot.move")
  prob.acceptance(result, "fixedEffects")
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











