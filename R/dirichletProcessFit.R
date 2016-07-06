#
# Functions for fitting the Dirichlet process approach to controlling
# for non-informative dropout
#
#

#' @include generics.R


#' 
#' Data stored with each iteration of the Dirichlet process model
#'
#' @title dirichlet.iteration
#' @param weights.mixing vector containing the probability of belonging to cluster k, 
#'      for k = 1 to the number of clusters
#' @param weights.conditional vector containing the probability of belonging to cluster k, given
#' that the subject was not in clusters 1 to k - 1
#' @param cluster.assignments current cluster assignments
#' @param betas A (k x 3) matrix of regression coefficients for the random intercept, slope,
#' and log of the dropout time for each cluster, with k = number of clusters
#' @param betas.deviations An (N x k) matrix of subject specific deviations from the cluster means
#' @param betas.covariates A (c x 1) vector of regression coefficients for covariates, 
#' with c = number of covariates
#' @param betas.covariates.mu A (c x 1) vector representing the mean of the distribution of 
#' regression coefficients related to covariates
#' @param betas.covariates.sigma A (c x c) vector representing the covariance of the distribution of 
#' regression coefficients related to covariates
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
#' @param expected.intercept expected value of the random intercepts
#' @param expected.slope expected value of the random slopes
#' @param slope.dropoutTimes estimated slopes by dropout times
#' @param density.intercept estimated density of the random intercepts
#' @param density.slope estimated density of the random slopes
#' 
#' @export dirichlet.iteration
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

#' dirichlet.model.options
#' 
#' Model options for the Dirichlet Process model.  Includes starting values, priors and
#' simulation/MCMC parameters.
#'
#' @param iterations number of iterations for the MCMC simulation
#' @param burnin burn in period for the simulation, i.e. the number of 
#' iterations to throw away at the beginning of the simulation
#' @param thin thinning interval, i.e. if thin=n, only keep every nth iteration
#' @param print printing interval, i.e. if print=n, print a counter on every nth iteration
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
#' @param sigma.error.tau Hyperprior for the residual error (Gaussian outcomes only)
#' @param density.intercept.domain vector of values at which to calculate the density of the random intercepts. If
#' NULL, no density will be calculated
#' @param density.slope.domain vector of values at which to calculate the density of the random slopes.  If NULL,
#' no density will be calculated
#' 
#' @importFrom matrixcalc is.positive.definite
#' 
#' @export dirichlet.model.options
#' 
dirichlet.model.options <- function(iterations=10000, burnin=500, thin=1, print=1,
                                    start.weights.mixing=NULL, n.clusters=15,
                                    dropout.estimationTimes=NULL, dropout.offset=0,
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
                                    sigma.error.tau = NULL,
                                    density.intercept.domain = NULL,
                                    density.slope.domain = NULL) {
  
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
    # print the ith iteration
    print = print,
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
    sigma.error.tau = sigma.error.tau,
    
    # density domains
    density.intercept.domain = density.intercept.domain,
    density.slope.domain = density.slope.domain
    
  )
  
  class(opts) <- append(class(opts), "dirichlet.model.options")
  return(opts)
  
}

#' 
#' Model fit for a Dirichlet Process model run
#' 
#' @title dirichlet.fit
#' @param model.options the original model options
#' @param dist the distribution of the outcome ("gaussian" or "binary")
#' @param groups the list of groups
#' @param covariates.var list of covariate column names
#' @param iterations the model run iterations after removing burn-in iterations and thinning
#'
#' @export dirichlet.fit
#' 
dirichlet.fit <- function(model.options, dist, groups, covariates.var, iterations) {
  fit = list(
    
    # store the model options for reference
    model.options = model.options,
    # distribution of the outcome    
    dist = dist,
    # list of group names
    groups = groups,
    # list of covariate names
    covariates.var = covariates.var,
    # list of dirichlet.iteration objects
    iterations = iterations
  )
  
  class(fit) <- append(class(fit), "dirichlet.fit")
  return(fit)
}


#' summary.dirichlet.fit
#' 
#' Summarize a Dirichlet model run
#'
#' @param fit dirichlet.fit object from an informativeDropout run.
#'
#' @export summary.dirichlet.fit
#' 
summary.dirichlet.fit <- function(fit, upper_tail=0.975, lower_tail=0.025) {
  if (length(fit$iterations) <= 0) {
    stop("No results found in Dirichlet model fit")
  }
  dist = fit$dist
  model.options = fit$model.options
  iterations = fit$iterations
  groups = fit$groups
  covariates.var = fit$covariates.var
  
  result.summary = list()
  
  for (group.index in 1:length(groups)) {
    group = groups[group.index]
    
    # total clusters for this group
    total_clusters = length(iterations[[1]]$cluster.N[[group.index]])
    
    ## subject specific effects
    # calculate mean, median, and 95% credible interval for the expected intercept
    expected_slope_sample = unlist(lapply(iterations, function(x) { 
      expected.slope = x$expected.slope[[group.index]]
      return (ifelse(is.null(expected.slope), 0, expected.slope))
    }))
    # calculate mean, median, and 95% credible interval for the expected intercept
    expected_intercept_sample = unlist(lapply(iterations, function(x) { 
      expected.intercept = x$expected.intercept[[group.index]]
      return (ifelse(is.null(expected.intercept), 0, expected.intercept))
    }))
    
    subject_specific_effects = data.frame(
      mean=c(mean(expected_slope_sample),mean(expected_intercept_sample)),
      median=c(median(expected_slope_sample), median(expected_intercept_sample)),
      ci_lower=c(quantile(expected_slope_sample, probs=lower_tail),
                 quantile(expected_intercept_sample, probs=lower_tail)),
      ci_upper=c(quantile(expected_slope_sample, probs=upper_tail),
                 quantile(expected_intercept_sample, probs=upper_tail))
    )
    row.names(subject_specific_effects) = c("slope", "intercept")
    
    ## covariance parameters
    total_covar = 16
    covariance_parameters = data.frame(
      mean=numeric(total_covar),
      median=numeric(total_covar),
      ci_lower=numeric(total_covar),
      ci_upper=numeric(total_covar)
    )
    row = 1
    covar_names <- vector()
    # summarize each component of the mean of the DP distribution
    for(i in 1:3) {
      covar_names = c(covar_names, paste("dp.dist.mu0[", i, "]", sep=""))
      param_sample <- unlist(lapply(iterations, function(x) { 
        return(x$dp.dist.mu0[[group.index]][i])
      }))
      covariance_parameters$mean[row] = mean(param_sample)
      covariance_parameters$median[row] = median(param_sample)
      covariance_parameters$ci_lower[row] = quantile(param_sample, probs=lower_tail)
      covariance_parameters$ci_upper[row] = quantile(param_sample, probs=upper_tail)
      row = row + 1
    }
    # summarize the covariance matrices (just the upper triangle)
    for(param in c("dp.dist.sigma0", "dp.cluster.sigma")) {
      for(i in 1:3) {
        for(j in i:3) {
          covar_names = c(covar_names, paste(param, "[", i, ",", j, "]", sep=""))
          param_sample <- unlist(lapply(iterations, function(x) { 
            return(x[[param]][[group.index]][i,j])
          }))
          covariance_parameters$mean[row] = mean(param_sample)
          covariance_parameters$median[row] = median(param_sample)
          covariance_parameters$ci_lower[row] = quantile(param_sample, probs=lower_tail)
          covariance_parameters$ci_upper[row] = quantile(param_sample, probs=upper_tail)
          row = row + 1
        }
      }
    }
    # add the DP concentration parameter
    covar_names <- c(covar_names, "dp.concentration")
    concentration_sample <- unlist(lapply(iterations, function(x) { 
      return(x$dp.concentration[[group.index]])
    }))
    covariance_parameters$mean[row] = mean(param_sample)
    covariance_parameters$median[row] = median(param_sample)
    covariance_parameters$ci_lower[row] = quantile(param_sample, probs=lower_tail)
    covariance_parameters$ci_upper[row] = quantile(param_sample, probs=upper_tail)
    row = row + 1
    # set the row names for the covariance param data frame
    row.names(covariance_parameters) = covar_names
    
    ## dropout time specific slopes
    slopes_by_dropout_time = data.frame(
      time=model.options$dropout.estimationTimes,
      mean=numeric(length(model.options$dropout.estimationTimes)),
      median=numeric(length(model.options$dropout.estimationTimes)),
      ci_lower=numeric(length(model.options$dropout.estimationTimes)),
      ci_upper=numeric(length(model.options$dropout.estimationTimes))
    )
    for (time.index in 1:(length(model.options$dropout.estimationTimes))) {
      slope_sample <- unlist(lapply(iterations, function(x) { 
        return(x$slope.dropoutTimes[[group.index]][time.index])
      }))
      slopes_by_dropout_time$mean[time.index] = mean(slope_sample)
      slopes_by_dropout_time$median[time.index] = median(slope_sample)
      slopes_by_dropout_time$ci_lower[time.index] = quantile(slope_sample, na.rm=TRUE, probs=lower_tail)
      slopes_by_dropout_time$ci_upper[time.index] = quantile(slope_sample, na.rm=TRUE, probs=upper_tail)
    }  
    row.names(slopes_by_dropout_time) = NULL
    
    # build the summary list
    result.summary[[paste("group", group, sep="_")]] = list(
      total_clusters = total_clusters,
      subject_specific_effects = subject_specific_effects,
      covariance_parameters = covariance_parameters,
      slopes_by_dropout_time = slopes_by_dropout_time
    )
  }
  
  ## summarize covariate effects
  if (!is.null(covariates.var)) {
    num_covariates = length(covariates.var)
    covariate_effects = data.frame(mean=numeric(num_covariates), median=numeric(num_covariates),
                                   ci_lower=numeric(num_covariates), ci_upper=numeric(num_covariates))
    for(i in 1:length(covariates.var)) {
      beta_covariate_sample = unlist(lapply(iterations, function(x) { 
        return(x$betas.covariates[i])
      }))
      covariate_effects$mean[i] = mean(beta_covariate_sample)
      covariate_effects$median[i] = median(beta_covariate_sample)
      covariate_effects$ci_lower[i] = quantile(beta_covariate_sample, probs=lower_tail)
      covariate_effects$ci_upper[i] = quantile(beta_covariate_sample, probs=upper_tail)
    }
    row.names(covariate_effects) <- covariates.var
    result.summary$covariate_effects = covariate_effects
  }
  
  if (dist == "gaussian") {
    # for gaussian outcomes, summarize the error variance
    sigma_error_sample = unlist(lapply(iterations, function(x) { 
      return(x$sigma.error)
    }))
    sigma_error_result = data.frame(
      mean = mean(sigma_error_sample),
      median = median(sigma_error_sample),
      ci_lower = quantile(sigma_error_sample, probs=lower_tail),
      ci_upper = quantile(sigma_error_sample, probs=upper_tail)
    )
    row.names(sigma_error_result) = "sigma.error"
    result.summary$sigma.error = sigma_error_result
  }
  
  return(result.summary)
  
}

#' plot.trace.dirichlet.fit
#' 
#' Produce a trace plot of the specified variable in Dirichlet model fit
#' 
#' @param fit the Dirichlet model fit 
#' @param type the type of information to plot.  
#' @param group (optional) the group to plot.  If not specified, separate plots will
#'    be produced for each group
#'
#' @export plot.trace.dirichlet.fit
#' 
plot.trace.dirichlet.fit <- function (fit, type="expectation", groups=NULL, params=NULL) {
  # determine the groups to plot
  if (is.null(groups)) {
    group_list <- fit$groups
  } else {
    group_list = groups
  }
  
  if (type == "expectation") {
    ### trace plot of expected slopes and intercepts
    # determine the parameters to plot
    if (is.null(params)) {
      params_list <- c("slope", "intercept")
    } else {
      params_list <- params
    }
    
    for(group_idx in 1:length(group_list)) {
      for(param_idx in 1:length(params_list)) {
        varname = paste("expected", params_list[param_idx], sep=".")
        sample = unlist(lapply(fit$iterations, function(x) { 
          return(x[[varname]][[group_idx]])
        }))
        ts.plot(sample, ylab=paste("Expected ", params_list[param_idx],
                                   ", group ", group_list[group_idx], sep=""))
      }
    }
    
  } else if (type == "betas.covariates") {
    ### trace plot of covariate effects
    if (is.null(fit$covariates.var)) {
      stop("No covariates were included in the specified model fit")
    }
    # determine which covariates to plot
    if (is.null(params)) {
      params_list <- fit$covariates.var
    } else {
      params_list <- params
    }
    
    for(i in 1:length(params_list)) {
      sample = unlist(lapply(fit$iterations, function(x) { 
        index = which(fit$covariates == params_list[i])
        return(x$betas.covariates[[index]])
      }))
      ts.plot(sample, ylab=paste(params_list[i], " effect", sep=""))
    }
    
  } else if (type == "clusters") {
    ### trace plot of number of non-empty clusters
    for(i in 1:length(group_list)) {
      sample = unlist(lapply(fit$iterations, function(x) { 
        return(sum(as.numeric(x$cluster.N[[i]] > 0)))
      }))
      ts.plot(sample, ylab=paste("Total Non-Empty Clusters, Group", group_list[i], sep=" "))
    }
    
  } else if (type == "covariance") {
    ### trace plot of covariance parameters
    if (is.null(params)) {
      params <- c("dp.dist.mu0", "dp.dist.sigma0", "dp.cluster.sigma", "dp.concentration")
      if (dist == "gaussian") {
        params <- c(params, "sigma.error")
      }
    }
    # plot group specific params
    for(group_idx in 1:length(group_list)) {
      group = group_list[group_idx]
      # plot the complonents of the mean of the DP distribution
      if ("dp.dist.mu0" %in% params) {
        for(i in 1:3) {
          group = group_list[group_idx]
          sample = unlist(lapply(fit$iterations, function(x) { 
            return(x$dp.dist.mu0[[group_idx]][i])
          }))
          ts.plot(sample, ylab=paste("dp.dist.mu0[", i,"], group ", group, sep=""))
        }
      }
      # plot the components of the covariance of the DP distribution and cluster covariance
      for(param in c("dp.dist.sigma0", "dp.cluster.sigma")) {
        if (param %in% params) {
          for(i in 1:3) {
            for(j in 1:3) {
              sample = unlist(lapply(fit$iterations, function(x) { 
                return(x[[param]][[group_idx]][i,j])
              }))
              ts.plot(sample, ylab=paste(param, "[", i, ",", j, "], group ", group, sep=""))
            }
          }
        }
      }
      # plot the concentration param
      if ("dp.concentration" %in% params) {
        sample = unlist(lapply(fit$iterations, function(x) { 
          return(x$dp.concentration[[group_idx]])
        }))
        ts.plot(sample, ylab=paste("dp.concentration, group ", group, sep=""))
      }
    }
    if (dist == "gaussian" && "sigma.error" %in% params) {
      sample = unlist(lapply(fit$iterations, function(x) { 
        return(x$sigma.error)
      }))
      ts.plot(sample, ylab="sigma.error")
    }
    
  } else {
    stop("Invalid type")
  }
  
}

#' plot.density.dirichlet.fit
#' 
#' density plot for a parameter
plot.density.dirichlet.fit <- function (fit, name="slope") {
  # take mean at each point across the iterations
  # histogram the means
}


#' plot.slopeByDropout.dirichlet.fit
#' 
#' plot the slope by dropout time
plot.slopeByDropout.dirichlet.fit <- function (fit, group=1, xlim=NULL, ylim=NULL,
                                               lower_tail=0.025, upper_tail=0.975) {
  
  slopes_by_dropout_time = data.frame(
    time=model.options$dropout.estimationTimes,
    mean=numeric(length(model.options$dropout.estimationTimes)),
    median=numeric(length(model.options$dropout.estimationTimes)),
    ci_lower=numeric(length(model.options$dropout.estimationTimes)),
    ci_upper=numeric(length(model.options$dropout.estimationTimes))
  )
  for (time.index in 1:(length(model.options$dropout.estimationTimes))) {
    slope_sample <- unlist(lapply(fit$iterations, function(x) { 
      return(x$slope.dropoutTimes[[group]][time.index])
    }))
    slopes_by_dropout_time$mean[time.index] = mean(slope_sample)
    slopes_by_dropout_time$median[time.index] = median(slope_sample)
    slopes_by_dropout_time$ci_lower[time.index] = quantile(slope_sample, probs=lower_tail)
    slopes_by_dropout_time$ci_upper[time.index] = quantile(slope_sample, probs=upper_tail)
  }  
  row.names(slopes_by_dropout_time) = NULL
  
  plot(slopes_by_dropout_time$time, slopes_by_dropout_time$mean, "l", xlim=xlim, ylim=ylim,
       xlab="Dropout time", ylab="Expected Slope")
  lines(slopes_by_dropout_time$time, slopes_by_dropout_time$ci_lower, "l", lty=3)
  lines(slopes_by_dropout_time$time, slopes_by_dropout_time$ci_upper, "l", lty=3)
}

# perform sensitivity analysis on the slope results
sensitivity.slope <- function(x, ...) {
  UseMethod("sensitivity.slope")
}


#'
#' Fit a bayesian model which accounts for dropout by modeling the 
#' releationship between dropout time and slope with natural cubic B-splines
#'
#' @title infortmativeDropout.bayes.dirichlet
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
#' @importFrom MCMCpack riwish
#' @importFrom gtools inv.logit
#' @importFrom Matrix Diagonal nearPD
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom Hmisc rMultinom
#'
#' @export
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
    if (is.na(model.options$betas.covariates.sigma) || 
        !is.matrix(model.options$betas.covariates.sigma) ||
        nrow(model.options$betas.covariates.sigma) != length(covariates.var) || 
        ncol(model.options$betas.covariates.sigma) != length(covariates.var) ||
        !is.positive.definite(model.options$betas.covariates.sigma)) {
      stop("Prior options error :: invalid prior variance for fixed effects related to covariates")
    }
  }
  
  if (dist == 'gaussian') {
    if (is.null(model.options$sigma.error.tau)) {
      stop("Prior options error :: invalid tau (hyperparameter) for the error variance")
    }
    if (is.null(model.options$sigma.error) || model.options$sigma.error <= 0) {
      stop("Prior options error :: invalid sigma.error")
    }
  }
  
  
  ## reorder the data and cache information on the groups and subjects
  # number of clusters
  n.clusters <- model.options$n.clusters
  
  # if no group is specified, create a bogus column with a single group value, 
  # making sure it doesn't already appear in the data set
  if (is.null(groups.var)) {
    groups.var <- "generated_group_0"
    idx <- 1
    while(groups.var %in% names(data)) {
      groups.var = paste("generated_group_", idx, collapse="")
      idx <- idx + 1
    }
    data[,groups.var] = rep(1, nrow(data))
  }
  
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
      return (rep(0,length(model.options$dropout.estimationTimes)))
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
                                                             data[,covariates.var, drop=FALSE], 
                                                             data[,outcomes.var, drop=FALSE])
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
  } else {
    sigma.error = NULL
  }
  
  # initialize the first model iteration
  iterations.saved = ceiling((model.options$iterations - model.options$burnin) / model.options$thin)
  modelIterationList <- vector(mode = "list", length = iterations.saved)
  
  model.previous = dirichlet.iteration(
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
    expected.intercept=rep(0, length(groupList)), expected.slope=rep(0, length(groupList)),
    slope.dropoutTimes = start.slope.dropoutTimes,
    sigma.error = sigma.error
  )
  
  iterations.saved.next = 1
  for (i in 1:model.options$iterations) {
    
    if (i %% model.options$print == 0) {
      print(paste("Dirichlet process model iteration = ", i, sep=""))
    }
    
    # make a copy of the previous iteration which will be modified as we move through the iteration
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
      eta <- rbeta(1, model.current$dp.concentration[[group.index]] + 1, group.n)
      # sample the concentration parameter from a mixture of Gamma distributions
      pi.eta <- ((model.options$dp.concentration.alpha + numNonEmptyClusters - 1) /
                   (group.n * (model.options$dp.concentration.beta - log(eta)) +
                      model.options$dp.concentration.alpha + numNonEmptyClusters - 1))
      # now draw the new concentration parameter for the Dirichlet process
      model.current$dp.concentration[[group.index]] <- 
        ((pi.eta * rgamma(1, model.options$dp.concentration.alpha + numNonEmptyClusters, 
                          model.options$sigma.error.tau - log(eta))) + 
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
                        dp.cluster.sigma.inv %*% ((group.betas[group.cluster.assignments==h,])))
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
              as.matrix(data.currentCluster[, covariates.var]) %*% matrix(model.current$betas.covariates)
          }
          posterior.mean = (posterior.var * 
                              ((prior.mean/prior.var) + 
                                 ((1/model.current$sigma.error) * 
                                    tapply(randomInts, data.currentCluster[, ids.var],sum)))) 
          # draw the random intercepts
          group.betas[group.cluster.assignments == h, 1] = 
            rnorm(group.cluster.N[h], posterior.mean, sqrt(posterior.var))
        } else {
          # binary distribution
          # Current Observation-Level Log-Likelihood
          eta <- (rep(group.betas[group.cluster.assignments == h, 1], group.nobs[group.cluster.assignments == h]) + 
                    data.currentCluster[, times.observation.var] * 
                    rep(group.betas[group.cluster.assignments == h, 2], group.nobs[group.cluster.assignments == h]))
          loglikelihood.previous <- rep(0, nrow(data.currentCluster))
          loglikelihood.previous[data.currentCluster[, outcomes.var] == 0] <- 
            log(1 - inv.logit(eta[data.currentCluster[, outcomes.var] == 0]))  
          loglikelihood.previous[data.currentCluster[, outcomes.var] == 1] <- 
            log(inv.logit(eta[data.currentCluster[, outcomes.var] == 1]))  
          # draw a proposal candidate beta from a symmetric random walk 
          beta.star <- group.betas[group.cluster.assignments == h, 1] + rnorm(group.cluster.N[h])
          eta.star <- (rep(beta.star, group.nobs[group.cluster.assignments == h]) + 
                         data.currentCluster[, times.observation.var] * 
                         rep(group.betas[group.cluster.assignments == h, 2], group.nobs[group.cluster.assignments == h]))
          prob.star <- inv.logit(eta.star)
          loglikelihood.star <- rep(0, nrow(data.currentCluster))
          loglikelihood.star[data.currentCluster[, outcomes.var] == 0] <- 
            log(1 - prob.star[data.currentCluster[, outcomes.var] == 0])  
          loglikelihood.star[data.currentCluster[, outcomes.var] == 1] <- 
            log(prob.star[data.currentCluster[, outcomes.var] == 1])  
          
          # subject level log likelihood
          subject.loglikelihood.previous <- tapply(loglikelihood.previous, data.currentCluster[, ids.var], sum)
          subject.loglikelihood.star <- tapply(loglikelihood.star, data.currentCluster[, ids.var], sum)
          
          ratio <- (subject.loglikelihood.star + 
                      dnorm(beta.star, prior.mean, sqrt(prior.var), log=TRUE) - 
                      (subject.loglikelihood.previous + 
                         dnorm(group.betas[group.cluster.assignments == h, 1], prior.mean, 
                               sqrt(prior.var), log=TRUE)))
          
          update <- 1 * (log(runif(group.cluster.N[h])) < ratio)
          group.betas[group.cluster.assignments == h, 1][update == 1] <- beta.star[update == 1]
        }
        
        ######  B1|B0 and U (Step 5 continued, draw random slope given intercept and u)
        inv.intU = solve(group.dp.cluster.sigma[-2,-2])
        covar = matrix(c(group.dp.cluster.sigma[1,2], group.dp.cluster.sigma[2,3]),1)
        # conditional prior mean/var
        prior.mean = as.numeric((means.currentCluster[2] + covar %*% inv.intU %*% 
                                   t(group.betas[group.cluster.assignments==h, c(1,3), drop=FALSE] - 
                                       subjectMeans.currentCluster[,c(1,3),drop=FALSE])))
        prior.var = as.numeric((group.dp.cluster.sigma[2,2] - covar %*% inv.intU %*% t(covar)))
        
        if (dist == "gaussian") {
          ## calculate the posterior mean and variance for the random intercept
          posterior.var = 1 / ((1 / prior.var) + 
                                 (tapply(data.currentCluster[, times.observation.var]^2, 
                                         data.currentCluster[, ids.var], sum) / 
                                    model.current$sigma.error))
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
                 as.matrix(data.currentCluster[, covariates.var]) %*% matrix(model.current$betas.covariates))
          }
          posterior.mean = (posterior.var * 
                              ((prior.mean/prior.var) + 
                                 ((1/model.current$sigma.error) * 
                                    tapply(randomSlopes, data.currentCluster[, ids.var],sum)))) 
          # draw the random intercepts
          group.betas[group.cluster.assignments == h, 2] = 
            rnorm(group.cluster.N[h], posterior.mean, sqrt(posterior.var))
        } else {
          # binary update
          # Current Observation-Level Log-Likelihood
          eta <- (rep(group.betas[group.cluster.assignments == h, 1], group.nobs[group.cluster.assignments == h]) + 
                    data.currentCluster[, times.observation.var] * 
                    rep(group.betas[group.cluster.assignments == h, 2], group.nobs[group.cluster.assignments == h]))
          loglikelihood.previous <- rep(0, nrow(data.currentCluster))
          loglikelihood.previous[data.currentCluster[, outcomes.var] == 0] <- 
            log(1 - inv.logit(eta[data.currentCluster[, outcomes.var] == 0]))  
          loglikelihood.previous[data.currentCluster[, outcomes.var] == 1] <- 
            log(inv.logit(eta[data.currentCluster[, outcomes.var] == 1]))  
          # draw a proposal candidate beta from a symmetric random walk 
          beta.star <- group.betas[group.cluster.assignments == h, 2] + rnorm(group.cluster.N[h])
          eta.star <- (rep(group.betas[group.cluster.assignments == h, 1], group.nobs[group.cluster.assignments == h]) + 
                         data.currentCluster[, times.observation.var] * 
                         rep(beta.star, group.nobs[group.cluster.assignments == h]))
          prob.star <- inv.logit(eta.star)
          loglikelihood.star <- rep(0, nrow(data.currentCluster))
          loglikelihood.star[data.currentCluster[, outcomes.var] == 0] <- 
            log(1 - prob.star[data.currentCluster[, outcomes.var] == 0])  
          loglikelihood.star[data.currentCluster[, outcomes.var] == 1] <- 
            log(prob.star[data.currentCluster[, outcomes.var] == 1])  
          
          # subject level log likelihood
          subject.loglikelihood.previous <- tapply(loglikelihood.previous, data.currentCluster[, ids.var], sum)
          subject.loglikelihood.star <- tapply(loglikelihood.star, data.currentCluster[, ids.var], sum)
          
          ratio <- (subject.loglikelihood.star + 
                      dnorm(beta.star, prior.mean, sqrt(prior.var), log=TRUE) - 
                      (subject.loglikelihood.previous + 
                         dnorm(group.betas[group.cluster.assignments == h, 2], prior.mean, sqrt(prior.var), log=TRUE)))
          
          update <- 1 * (log(runif(group.cluster.N[h])) < ratio)
          group.betas[group.cluster.assignments == h, 2][update == 1] <- beta.star[update == 1]
        }
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
      if (!is.null(model.options$density.intercept.domain)) {
        sim.int = vector()
        for(h in 1:n.clusters) {
          sim.int <-rbind(sim.int,
                          (clusterProbabilities[h] * 
                             dnorm(model.options$density.intercept.domain, 
                                   model.current$cluster.mu[[group.index]][h,1], 
                                   sqrt(model.current$dp.dist.sigma0[[group.index]][1,1]))))
        }
        model.current$density.intercept[[group.index]] <- colSums(sim.int)
      }
      if (!is.null(model.options$density.slope.domain)) {
        sim.slope = vector()
        for (h in 1:n.clusters) {
          sim.slope <-rbind(sim.slope,
                            (clusterProbabilities[h] * 
                               dnorm(model.options$density.slope.domain, 
                                     model.current$cluster.mu[[group.index]][h,2], 
                                     sqrt(model.current$dp.dist.sigma0[[group.index]][2,2]))))
        }
        model.current$density.slope[[group.index]] <- colSums(sim.slope)
      }
      
      # Equation 19 
      #Estimate slope at each dropout time
      if (!is.null(model.options$dropout.estimationTimes)) {
        tmp = matrix(0,length(model.options$dropout.estimationTimes), n.clusters)
        for (h in 1:n.clusters) {
          tmp[,h] = (
            group.weights.mixing[h] * 
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
                     (log(dropoutTime) - model.current$cluster.mu[[group.index]][,3]) /
                     covar[3,3]))
          )
        }
      } else {
        model.current$slope.dropoutTimes = NULL
      }
      
      # Estimate E(B1), E(B0) - section 2.4 expectation of random effects (int, slope)
      model.current$expected.intercept[[group.index]] = 
        sum(group.weights.mixing * model.current$cluster.mu[[group.index]][,1])
      model.current$expected.slope[[group.index]] = 
        sum(group.weights.mixing * model.current$cluster.mu[[group.index]][,2])
      
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
    
    # update fixed effects associated with covariates
    if (!is.null(covariates.var)) {
      if (dist == "gaussian") {
        sigma.error.inv = 1/model.current$sigma.error
        prior.sigma.inv = solve(model.options$betas.covariates.sigma)
        
        residuals = data[, outcomes.var] - intercepts - slopes * data[, times.observation.var]
        var = solve(prior.sigma.inv + crossprod(as.matrix(data[,covariates.var])) * sigma.error.inv)
        
        m <- var %*% (crossprod(as.matrix(data[, covariates.var]), residuals) * sigma.error.inv)
        model.current$betas.covariates = rmvnorm(1,m,var)
        
      } else {
        # binary case, need to do metropolis hastings
        # build components of eta
        y <- as.matrix(data[, outcomes.var])
        C = as.matrix(data[,covariates.var])
        cBeta = as.vector(C %*% model.current$betas.covariates)
        XTheta.previous = intercepts + slopes * data[, times.observation.var]
        # calculate the previous eta and associated probability
        eta.previous = as.vector(XTheta.previous + cBeta)
        prob.previous = inv.logit(eta.previous)
        # adjust probabilities within tolerance levels
        prob.previous[prob.previous < model.options$prob.min] <-  model.options$prob.min
        prob.previous[prob.previous > model.options$prob.max] <-  model.options$prob.max
        loglikelihood.previous <- sum(log((1 - prob.previous[y==0]))) + sum(log(prob.previous[y==1]))    
        
        # build y-tilde
        r0.inv = solve(model.options$betas.covariates.sigma)
        yt = cBeta + (y - prob.previous) * (1 / (prob.previous * (1 - prob.previous)))
        # build the weight matrix
        weight <- Diagonal(x = prob.previous * (1-prob.previous))
        # build the covariance of the beta coefficients and make sure it is positive definite
        covar <- as.matrix(nearPD(solve(r0.inv + crossprod(C, as.matrix(weight %*% C))))$mat)
        # build the mean
        mu <- covar  %*% (crossprod(C, as.matrix(weight %*% yt)))
        # draw the proposed coefficients for the fixed effects
        beta.covariates.star <- rmvnorm(1, mu, covar)
        cBeta.star <- as.vector(C %*% beta.covariates.star)
        # get proposal probabilities
        eta.star <- XTheta.previous + C %*% beta.covariates.star
        prob.star = inv.logit(eta.star)
        # adjust probabilities within tolerance levels
        prob.star[prob.star < model.options$prob.min] <-  model.options$prob.min
        prob.star[prob.star > model.options$prob.max] <-  model.options$prob.max
        loglikelihood.star <- sum(log((1 - prob.star[y == 0]))) + sum(log(prob.star[y == 1]))    
        
        y.star <- cBeta.star + (y - prob.star) * (1 / (prob.star * (1 - prob.star)))
        weight.star <- Diagonal(x=prob.star * (1 - prob.star))
        covar.star <- as.matrix(nearPD(solve(r0.inv + crossprod(C, as.matrix(weight.star %*% C))))$mat)
        mu.star <- covar.star %*% (crossprod(C,as.matrix(weight.star%*%y.star)))
        
        # Metropolis hastings step
        rho <- (loglikelihood.star - loglikelihood.previous + 
                  log(dmvnorm(as.vector(model.current$betas.covariates), as.vector(mu.star), covar.star)) - 
                  log(dmvnorm(as.vector(beta.covariates.star), as.vector(mu), covar)) +
                  log(dmvnorm(as.vector(beta.covariates.star), as.vector(model.options$betas.covariates.mu), model.options$betas.covariates.sigma)) -
                  log(dmvnorm(as.vector(model.current$betas.covariates), as.vector(model.options$betas.covariates.mu), model.options$betas.covariates.sigma)))
        
        
        
        if (rho > log(runif(1))) {
          model.current$betas.covariates = betaCovariates.star
        }
      }
    }
    
    if (dist == "gaussian") {
      ## TODO: remove group specific variables, use betas (random effects)
      # Step 9 Update sigma.error with inverse gamma
      residual = data[, outcomes.var] - (intercepts + data[, times.observation.var] * slopes)
      if (!is.null(covariates.var)) {
        residual = residual - as.vector(as.matrix(data[,covariates.var]) %*% 
                                          matrix(model.current$betas.covariates))
      }
      g <- model.options$sigma.error.tau + crossprod(residual)/2
      tau <- rgamma(1, model.options$sigma.error.tau + n.total / 2, g)
      model.current$sigma.error = 1 / tau
    } 
    
    # save the current iteration
    model.previous = model.current
    if ((i %% model.options$thin == 0 && i > model.options$burnin &&
         iterations.saved.next <= length(modelIterationList)) || 
        iterations.saved.next == length(modelIterationList)) {
      modelIterationList[[iterations.saved.next]] = model.current
      iterations.saved.next = iterations.saved.next + 1
    }
    
    
  } # END ITERATION LOOP
  
  return(dirichlet.fit(model.options, dist, groupList, covariates.var, modelIterationList))
  
} # END FUNCTION 





