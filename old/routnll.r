# routnll: neg log lik for routing | v0.2
# * Change log:
#    - v0.2: added option in datalist for routingpred: "add", "replace", or
#            "convcomb" for a convex combination with an extra param
#    - v0.1: initial version

# * Convention for vector stacking: time = outer loop, space = inner loop

# # # scaled logistic transfo R -> (lb,ub)
# # logit_transfo <- function(x,lb=0,ub=1){lb+(ub-lb)/(1+exp(-x))}
# # # scaled inverse logistic transfo (lb,ub) -> R
# # logit_unstranfo <- function(x,lb=0,ub=1){-log((ub-lb)/(x-lb)-1)}
# transfo01 <- function(x){1/(1+exp(-x))} # R -> (0,1)
# untransfo01 <- function(p){log(p/(1-p))} # (0,1) -> R

# alr <- function(pi){ # maps (n-1)-Simplex to R^(n-1), length(pi) > 1
# 	return(log(pi[-length(pi)]/pi[length(pi)])) # last pi as reference
# }
# alrinv <- function(alrpi){ # maps R^(n-1) to (n-1)-Simplex
# 	alrpi <- ifelse(alrpi>400,400,alrpi) # otherwise explodes, 1 anyway
# 	pi <- exp(c(alrpi,0)) # last pi as reference
# 	return(pi/sum(pi)) # closure operation
# }
alrinv <- function(alrpi){ # maps R^(n-1) to (n-1)-Simplex
	# alrpi <- ifelse(alrpi>400,400,alrpi) # otherwise explodes, 1 anyway
	pi <- exp(c(alrpi,0)) # last pi as reference
	return(pi/sum(pi)) # closure operation
}


rout_nll <- function(par){
	# * data in R global envir assumed to be stored in datalist
	# * par vector order:
	#   - log_sigma (1)
	#   - log_wscale (1)
	#   - wshapebeta (p)
	#   - alr_prob (1)
	
	
	
	#----------------------------------------------------------------------------#
	# Inputs
	#----------------------------------------------------------------------------#
	
	getAll(par, datalist)
	# ^ datalist must include:
	#   * obsmat: numeric matrix of discharge (m3/s), can incl NAs, nS x nT
	#   * obsindmat: Boolean matrix of same dim as obsmat
	#   * predmat: numeric matrix of pred discharge, no NAs, same dim as obsmat
	#   * maxlag: integer >=1, max lag in routing from spatial neighbors
	#   * routingorder: integer vector routing ustr -> dstr, within 1:nS but of
	#                   length <=nS since excl loc most ustr where fitted=predmat
	#   * neighlist: list of length nS, each element being an integer vector of
	#                values within 1:nS of direct/Markov neighbors ustr, not incl
	#                itself, so each entry is of length within {0, ..., nS-1}
	#   * wshapecovlist: list of length nS, same ordering as neighlist, each
	#                    element being a list of numeric vectors (number of
	#                    vectors is length of the corresponding element in
	#                    neighlist), each vector being of same length p =
	#                    length(wshapebeta) where the values are static cov
	#   * lag0: lag 0 = same day for gamma kernel, but >0, e.g. 1e-3
	
	
	#----------------------------------------------------------------------------#
	# Setup and init
	#----------------------------------------------------------------------------#
	
	nS <- nrow(predmat) # nb loc overall (polyg and stations)
	nT <- ncol(predmat) # nb total time steps
	
	sigma <- exp(log_sigma)
	wscale <- exp(log_wscale) # cst, gamma shape has all the cov
	probvec <- alrinv(alr_prob) # dim 2
	
	
	#----------------------------------------------------------------------------#
	# routing
	#----------------------------------------------------------------------------#
	
	fitted <- predmat # matrix(data=qst,nrow=nS,ncol=nT)
	
	if (routingpred=='add'){
		for (t in (1+maxlag):nT){
			# excl first maxlag*nS obs for routing
			for (s in routingorder){
				# loop over all loc, excl the ones most ustr where fitted[s,]=predmat[s,]
				# if (length(neighlist[[s]])>0){
				for (ss in 1:length(neighlist[[s]])){ # loop over direct ustr neighbors
					# gammadens <- dgamma(x=c(lag0,1:maxlag), # lag0 = same day = e.g. 1e-3
					# 										shape=wshape*distlist[[s]][ss],scale=wscale)
					whshape_ss <- exp(wshapebeta%*%wshapecovlist[[s]][[ss]])
					# ^ lin comb (p->1) with log link for shape>0
					gammadens <- dgamma(x=c(lag0,1:maxlag), # lag0 = same day = e.g. 1e-3
															shape=whshape_ss,
															scale=wscale)
					gammadens <- gammadens/sum(gammadens) # rescale, so sum(weights) = 1
					
					fitted[s,t] <- fitted[s,t] +
						+ sum(gammadens*fitted[neighlist[[s]][[ss]], t-(0:maxlag)]) # add
					# fitted[s,t] <- sum(gammadens*fitted[neighlist[[s]][[ss]], t-(0:maxlag)])
					# ^ add to predmat or replace it completely
					
					# ^ weighted comb of pred from ustr neighbors at lag 0:maxlag
				}
				# } # else no loc ustr from s and thus fitted[s,] = predmat[s,]
			}
		}
	} else if (routingpred=='replace'){
		for (t in (1+maxlag):nT){
			# excl first maxlag*nS obs for routing
			for (s in routingorder){
				# loop over all loc, excl the ones most ustr where fitted[s,]=predmat[s,]
				# if (length(neighlist[[s]])>0){
				for (ss in 1:length(neighlist[[s]])){ # loop over direct ustr neighbors
					# gammadens <- dgamma(x=c(lag0,1:maxlag), # lag0 = same day = e.g. 1e-3
					# 										shape=wshape*distlist[[s]][ss],scale=wscale)
					whshape_ss <- exp(wshapebeta%*%wshapecovlist[[s]][[ss]])
					# ^ lin comb (p->1) with log link for shape>0
					gammadens <- dgamma(x=c(lag0,1:maxlag), # lag0 = same day = e.g. 1e-3
															shape=whshape_ss,
															scale=wscale)
					gammadens <- gammadens/sum(gammadens) # rescale, so sum(weights) = 1
					
					# fitted[s,t] <- fitted[s,t] +
					# 	+ sum(gammadens*fitted[neighlist[[s]][[ss]], t-(0:maxlag)]) # add
					fitted[s,t] <- sum(gammadens*fitted[neighlist[[s]][[ss]], t-(0:maxlag)])
					# ^ add to predmat or replace it completely
					
					# ^ weighted comb of pred from ustr neighbors at lag 0:maxlag
				}
				# } # else no loc ustr from s and thus fitted[s,] = predmat[s,]
			}
		}
	} else if (routingpred=='convcomb'){
		for (t in (1+maxlag):nT){
			# excl first maxlag*nS obs for routing
			for (s in routingorder){
				# loop over all loc, excl the ones most ustr where fitted[s,]=predmat[s,]
				# if (length(neighlist[[s]])>0){
				for (ss in 1:length(neighlist[[s]])){ # loop over direct ustr neighbors
					# gammadens <- dgamma(x=c(lag0,1:maxlag), # lag0 = same day = e.g. 1e-3
					# 										shape=wshape*distlist[[s]][ss],scale=wscale)
					whshape_ss <- exp(wshapebeta%*%wshapecovlist[[s]][[ss]])
					# ^ lin comb (p->1) with log link for shape>0
					gammadens <- dgamma(x=c(lag0,1:maxlag), # lag0 = same day = e.g. 1e-3
															shape=whshape_ss,
															scale=wscale)
					gammadens <- gammadens/sum(gammadens) # rescale, so sum(weights) = 1
					
					# fitted[s,t] <- fitted[s,t] +
					# 	+ sum(gammadens*fitted[neighlist[[s]][[ss]], t-(0:maxlag)]) # add
					# fitted[s,t] <- sum(gammadens*fitted[neighlist[[s]][[ss]], t-(0:maxlag)])
					# ^ add to predmat or replace it completely
					fitted[s,t] <- probvec[1]*fitted[s,t] +
						+ probvec[2]*sum(gammadens*fitted[neighlist[[s]][[ss]], t-(0:maxlag)])
					# ^ convex combination with extra param in probvec
					
					# ^ weighted comb of pred from ustr neighbors at lag 0:maxlag
				}
				# } # else no loc ustr from s and thus fitted[s,] = predmat[s,]
			}
		}
	} else {stop('routingpred must be either "add", "replace" or "convcomb".')}
	
	
	
	
	#----------------------------------------------------------------------------#
	# pnll eval at fitted
	#----------------------------------------------------------------------------#
	
	pnll <- 0 # ini pen neg loglik
	
	for (s in 1:nS){ # loop over all loc
		for (t in (1+maxlag):nT){ # loop over all time points after burn-in
			# time = outer loop, space = inner loop
			if (obsindmat[s,t]){ # lkhd contrib only when obs available (tr set)
				pnll <- pnll - dnorm(x=obsmat[s,t], mean=fitted[s,t], sd=sigma, log=T)
			}
		}
	}
	
	
	#----------------------------------------------------------------------------#
	# Outputs
	#----------------------------------------------------------------------------#
	
	REPORT(sigma)
	REPORT(wscale)
	REPORT(wshapebeta)
	REPORT(probvec) # report even if not used, routingpred = "add" or "replace"
	
	REPORT(fitted) # fitted values on modeling scale
	
	return(pnll)
}
# end routnll
