#
# Demo reproducing the results of the manuscript
#
# Moore CM, Carlson NE, Kreidler SM, MaWhinney S.
# "A Bayesian natural cubic B-spline varying coefficient method for non-ignorable dropout.",
# 2018, under review.
#
#
require(pbapply)

##################################################
# Analysis of untreated subjects
##################################################
# load the data set of untreated hiv participants
data("untreated_hiv")

# Calculate dropout year
untreated_hiv$drop_years = untreated_hiv$drop_day / 365

# set model options
# Note the number of iterations has been reduced so the example code will
# run quickly.  The models may not converge with this number of iterations.  
# Due to the stochastic nature of MCMC, results may not be reproduced exactly. 
model.options <- bayes.splines.model.options(
  iterations=10000, burnin=5000, thin=1, print=1000,
  knots.prob.birth=0.2, knots.min=1, knots.max=10, knots.stepSize=15*3,
  knots.positions.start=list(c(142, 562, 982, 1402, 1822)/365, c(139, 559, 979, 1399, 1819)/365),
  knots.positions.candidate=list(seq(142,1826,15)/365,seq(139,1826,15)/365),
  dropout.estimationTimes=c(1,2,3),
  sigma.error=1,
  sigma.error.shape.tau=0.001, sigma.error.rate.tau=0.001,
  sigma.beta=25, lambda.numKnots=5,
  sigma.residual=1.25^2,
  sigma.randomIntercept = 1,
  sigma.randomSlope = 1,
  sigma.randomInterceptSlope = 0,
  sigma.randomEffects.df = 3,
  sigma.randomEffects.scale = diag(2))

# select columns to include in the model
ids.var = "patid"
outcomes.var = "logcd4"
groups.var = "hard_drug"
covariates.var = c("baselogcd4", "baselogcd4_years")
times.dropout.var = "drop_years"
times.observation.var = "years"

# set the model fitting method
method="bayes.splines"

# select the distribution of the outcome
dist = "gaussian"

# set a random seed
set.seed(1066)

# fit the model
fit = informativeDropout(untreated_hiv, model.options, ids.var, outcomes.var, groups.var, covariates.var, 
                         times.dropout.var, times.observation.var, 
                         method=method, dist=dist) 

# summarize the result
summary(fit)

# To calculate the slope for subjects in each group assuming baseline CD4 of 478.5:
others.slope<-do.call(rbind, pblapply(fit$iterations, function(x){
  x$slope.marginal[[1]]+x$betas.covariates[2]*log(478.5)
}))

mean(others.slope)
quantile(others.slope, c(0.025, 0.975))

druguse.slope<-do.call(rbind, pblapply(fit$iterations, function(x){
  x$slope.marginal[[2]]+x$betas.covariates[2]*log(478.5)
}))

mean(druguse.slope)
quantile(druguse.slope, c(0.025, 0.975))

# Calculate and summarize difference in slopes between groups
diff <- druguse.slope - others.slope

mean(diff)
quantile(diff, c(0.025, 0.975))

# Calculate Posterior P Value for difference in slopes < 0
length(diff[diff<0])/length(diff)

# sensitivity analysis
# Estimate slopes at 1 to 4 years, if we adjust by 
# deltas of 0, 0.25, 0.5, and 0.75
#
# To perform sensitivity analysis, we first create a dataframe with one row
# per subject. 

# number of observations per subject
numObservations = sapply(unique(untreated_hiv[,ids.var]), function(id) {
  return(nrow(untreated_hiv[untreated_hiv[,ids.var] == id,]))
})
# index of the first observation per subject
firstObsPerSubject = c(1,1+cumsum(numObservations)[-length(numObservations)])
data.onePerSubject = untreated_hiv[firstObsPerSubject,]

# we perform the sensitivity analysis assuming a baseline CD4 count of 478.5
data.onePerSubject$baselogcd4 = log(478.5)

# for time interacted covariates, set the value of the covariate*time when time = 1
# so they can be included as part of the slope
data.onePerSubject$baselogcd4_years = log(478.5)

sens.fit = sensitivity(fit, 
                       times.estimation=(1:4), deltas=c(0, 0.25, 0.5, 0.75),
                       data.onePerSubject, 
                       times.dropout.var="drop_years", group.var=groups.var, 
                       covariates.time.var=c("baselogcd4_years"), 
                       covariates.nontime.var=c("baselogcd4"))
summary(sens.fit)

#######################################################
# Analysis of Treated Subjects: Viral Load Suppression
#######################################################
# load the data set
data("treated_hiv")

# Get set of candidate knot locations and specify starting positions
candidates <- seq(min(treated_hiv$drop_year), max(treated_hiv$drop_year),1/12)
starting <- candidates[c(seq(1, length(candidates), length(candidates)/4),length(candidates))]

# set the model options
# Note the number of iterations has been reduced so the example code will
# run quickly.  The models may not converge with this number of iterations.  
# Due to the stochastic nature of MCMC, results may not be reproduced exactly. 
model.options <- bayes.splines.model.options(iterations=10000, burnin=5000, thin=1,print=1000,
                                             knots.prob.birth=0.2, knots.min=1, knots.max=10, 
                                             knots.stepSize=0.25,
                                             knots.positions.start=list(starting, starting), 
                                             knots.positions.candidate=list(candidates, candidates),
                                             dropout.estimationTimes=seq(1,11,2),
                                             sigma.error.shape.tau=0.001, sigma.error.rate.tau=0.001,
                                             sigma.beta=100, lambda.numKnots=5,
                                             sigma.residual = 2.5^2,
                                             sigma.error=1,
                                             sigma.randomIntercept = 1,
                                             sigma.randomSlope = 1,
                                             sigma.randomInterceptSlope = 0.001,
                                             sigma.randomEffects.df = 3,
                                             sigma.randomEffects.scale = diag(2),
                                             eta.null=NULL)

# set the columns to use in the model
ids.var = "patid"
outcomes.var = "suppressed"
groups.var = "rec_drug"
covariates.var = c("baselogcd4", "baselog10vl", "baselogcd4_years", "baselog10vl_years")
times.dropout.var = "drop_year"
times.observation.var = "years"

# set the model fitting method
method="bayes.splines"
# set the distribution of the outcome
dist = "binary"

# set a random seed
set.seed(1066)
# fit the model
fit = informativeDropout(treated_hiv, model.options, ids.var, 
                         outcomes.var, groups.var,
                         covariates.var, 
                         times.dropout.var, times.observation.var,
                         method=method, dist=dist)

# summarize the result
summary(fit)

# To calculate slope for subject with baseline CD4 of 267, log10VL of 4.2:
others.slope<-do.call(rbind, pblapply(fit$iterations, function(x){
  x$slope.marginal[[1]]+x$betas.covariates[3]*log(267) + x$betas.covariates[4]*4.2
}))

quantile(others.slope, c(0.025, 0.975))
mean(others.slope[20000:49999])

druguse.slope<-do.call(rbind, pblapply(fit$iterations, function(x){
  x$slope.marginal[[2]]+x$betas.covariates[3]*log(267) + x$betas.covariates[4]*4.2
}))

quantile(druguse.slope, c(0.025, 0.975))
mean(druguse.slope[20000:49999])

covar <- do.call(rbind, pblapply(fit$iterations, function(x){
  x$betas.covariates
}))

# Calculate and summarize difference in slopes between groups
diff <- druguse.slope - others.slope

summary(diff)
quantile(diff, c(0.025, 0.975))

# Calculate Posterior P Value for difference in slopes < 0
length(diff[diff<0])/length(diff)

# sensitivity analysis
# Estimate slopes at 1 to 4 years, if we adjust by 
# deltas of 0, 0.25, 0.5, and 0.75
#
# To perform sensitivity analysis, we first create a dataframe with one row
# per subject. 

# number of observations per subject
numObservations = sapply(unique(treated_hiv[,ids.var]), function(id) {
  return(nrow(treated_hiv[treated_hiv[,ids.var] == id,]))
})
# index of the first observation per subject
firstObsPerSubject = c(1,1+cumsum(numObservations)[-length(numObservations)])
data.onePerSubject = treated_hiv[firstObsPerSubject,]

# we perform the sensitivity analysis assuming a baseline CD4 count of 267
# and baseline log10vl of 4.2
data.onePerSubject$baselogcd4 = log(267)
data.onePerSubject$baselog10vl = 4.2

# for time interacted covariates, set the value of the covariate*time when time = 1
# so they can be included as part of the slope
data.onePerSubject$baselogcd4_years = log(267)
data.onePerSubject$baselog10vl_years = 4.2

sens.fit = sensitivity(fit, 
                       times.estimation=(1:4), deltas=c(0, 0.25, 0.5, 0.75),
                       data.onePerSubject, 
                       times.dropout.var=times.dropout.var, group.var=groups.var, 
                       covariates.time.var=c("baselogcd4_years", "baselog10vl_years"), 
                       covariates.nontime.var=c("baselogcd4", "baselog10vl"))
summary(sens.fit) 
