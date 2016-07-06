#
# Functions / classes to support the mixed model spline fit
#
#

#' Model options for the mixed model fit.  A convenience wrapper
#' around lme
#'
#' @param knots.positions knot positions
#' 
#' @export mixed.model.options
#' 
mixed.model.options <- function(knots.positions=NULL) {
  
  opts = list(
    # knot positions
    knots.positions = knots.positions
  )
  
  class(opts) <- append(class(opts), "mixed.model.options")
  return(opts)
  
}

#'
#' infortmativeDropout.mixed
#'
#' Fit a mixed model which accounts for dropout by modeling the 
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
#' @export
#'
informativeDropout.mixed <- function(data, ids.var, outcomes.var, groups.var, 
                                     covariates.var, times.dropout.var, times.observation.var, 
                                     dist, model.options) {
  
  # make sure we have the correct type of mcmc opts
  if (!is(model.options,"mixed.model.options")) {
    stop("Model options error :: options must be of type mixed.model.options")
  }
  
  if (!is.null(groups.var)) {
    
  }
  knots = model.options$knots
  # get the coefficients for the knots
#   
#   formula.string <- paste(
#     paste(outcomes.var, "~"),
#     ns
#     paste("ns(", times.dropout.var, ", knots=c(", paste(knots.interior, collapse=",", Boundary.knots=knots.boundary,
#        "intercept=T)" 
#     
#     paste(
#       paste(covariates.var, collapse=" + "), 
#       paste("(1+", times.observation.var, "|", ids.var, ")", sep=""),
#       sep = " + "),
#     sep=" ")
#   
  
  if (dist == "gaussian") {
    fit = lmer(as.formula(formula.string), data=data, REML=FALSE)
    
  } else {
    
  }
  
  return(fit)
  
} # END FUNCTION 
