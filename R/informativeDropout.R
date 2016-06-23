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

#' @include bayesianSplineFit.R
#' @include dirichletProcessFit.R
#' @include mixedModelFit.R

#' 
#' Fit a varying coefficient model for longitudinal studies with
#' informative dropout. 
#' 
#' @title informativeDropout
#' @param data the data set 
#' @param ids.var column of the data set containing participant identifiers
#' @param outcomes.var column of the data set containing the outcome variable
#' @param groups.var column of the data set indicating the treatment groups
#' @param covariates.var list of columns in the data set containing covariates
#' @param times.dropout.var column of the data set containing dropout times
#' @param times.observation.var column of the data set containing 
#' @param method the modeling method
#' @param dist the distribution of the outcome, valid values are "gaussian" or "binary"
#' @param model.options model options (see dirichlet.model.options for details)
#' 
#' @export informativeDropout
#' 
informativeDropout <- function(data, ids.var, outcomes.var, groups.var, covariates.var, 
                               times.dropout.var, times.observation.var,
                               method="bayes.splines", dist="normal",
                               model.options) {
  
  if (method == 'bayes.splines') {
    # model the relationship between dropout time and slope using natural splines
    return (informativeDropout.bayes.splines(data, ids.var, outcomes.var, groups.var, covariates.var, 
                                             times.dropout.var, times.observation.var, dist, model.options))
    
  } else if (method == 'dirichlet') {
    # account for informative dropout using a dirichlet process 
    return (informativeDropout.bayes.dirichlet(data, ids.var, outcomes.var, groups.var, covariates.var, 
                                               times.dropout.var, times.observation.var, dist, model.options))
  } else if (method == 'mixed.splines') {
    # fit a mixed model which models the relationship between dropout time and slope using natural splines
    return (informativeDropout.mixed(data, ids.var, outcomes.var, groups.var, covariates.var, 
                                     times.dropout.var, times.observation.var, dist, 
                                     model.options))
  } else {
    stop(paste("Invalid method: ", method, sep=""))
  }
  
}



