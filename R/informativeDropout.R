#
# 
#
#
#
#


#
# Reverse Jump MCMC 
#
reverseJumpMCMC <- function(iterations, pAddKnot, knots.max, knots.min,
                            knots.startPositions, knots.candidatePositions) {
  
  knotPositions <- knots.startPositions
  
  for (i in 1:iterations) {
    # randomly decide to add/remove a knot
    if (runif(1) < pAddKnot) {
      # add a knot by randomly selecting a candidate knot
      index = sample(1:length(knots.candidatePositions), 1)
      knotPositions <- sort(c(knotPositions, knots.candidatePositions[index]))
    } else {
      # randomly remove an existing knot
      index = sample(1:length(knotPositions), 1)
      knotPositions <- knotPositions[-index]
    }
    
    # estimate fixed effects with new set of knots
    
    
  }
}

#
#
#
#
#
informativeDropout <- function() {
  
 
  # return an nsv class
}


# functions
# - summary: show estimates
# - print
# - predict: get distribution for 


