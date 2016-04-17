#' Global specs for rcsurplus
#' 
#' \code{rcsurplus1d.defaults} is a list that contains default 
#' parameters for this package
#' @param priorK prior for carrying capacity
#' @param priorr prior for intrinsic growth rate
#' @param priorq prior for catchability
#' @param priors prior for observational variance
#' @param mcmc_n number of Markov Chain Monte Carlo iterations
#' @param mcmc_t MCMC thinning
#' @param mcmc_b MCMC burn-in fraction
#' @param mcmc_c number of MCMC chains
#' @export rcsurplus1d.default
rcsurplus1d.default <- list(
      priorK   = c(100, 15000), 
	  priorr   = c(0,1), 
      priorPHI = c(0,1), 
	  priorq   = c(-20,-1), 
	  priors   = c(-20,20),
      mcmc_n   = c(1000, 10000),
      mcmc_t   = c(1,100),
      mcmc_b   = c(0.05, 0.5),
      mcmc_c   = c(1,3)   
  )

#' Function name
#'		
#' \code{function_name} is a function that does
#' @param par1 Explanation of the role of parameter 1
#' @return theExplanation of what is returned (z)
function_name <- function( par1 ){
	print('working on it')
    z <- 0
	return(z)
}

