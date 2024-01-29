# routnll blockwise: neg log lik for routing | v0.2
# * Change log:
#    - v0.2: initial version, forked from routmod/rout_nll v0.2

# * Convention for vector stacking: time = outer loop, space = inner loop

rout_nll_block_1 <- function(par){
	# * data in R global envir assumed to be stored in datalist
	# * par vector order:
	#   - log_sigma (1)
	#   - log_wscale (1)
	#   - wshapebeta (p)
	
	#----------------------------------------------------------------------------#
	# Inputs
	#----------------------------------------------------------------------------#
	
	# getAll(par, data)
	getAll(par)
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
	
	nS <- nrow(datalist$predmat) # nb loc overall (polyg and stations)
	nT <- ncol(datalist$predmat) # nb total time steps
	
	# sigma <- exp(parvec[1])
	# wscale <- exp(parvec[2]) # cst, gamma shape has all the cov
	# wshapebeta <- parvec[3:length(parvec)]
	sigma <- exp(log_sigma)
	wscale <- exp(log_wscale) # cst, gamma shape has all the cov
	
	
	#----------------------------------------------------------------------------#
	# routing
	#----------------------------------------------------------------------------#
	
	fitted <- datalist$predmat # matrix(data=qst,nrow=nS,ncol=nT)
	
	for (t in (1+datalist$maxlag):nT){
		# for (t in 1:nT){ # blockwise after block 1: nT=maxlag, keep first maxlag*nS obs
		for (s in datalist$routingorder){
			# loop over all loc, excl the ones most ustr where fitted[s,]=predmat[s,]
			for (ss in 1:length(datalist$neighlist[[s]])){ # direct ustr neighbors
				whshape_ss <- exp(wshapebeta%*%datalist$wshapecovlist[[s]][[ss]])
				# ^ lin comb (p->1) with log link for shape>0
				gammadens <- dgamma(x=c(datalist$lag0,# lag0 = same day = e.g. 1e-3
																1:datalist$maxlag),
														shape=whshape_ss,
														scale=wscale)
				gammadens <- gammadens/sum(gammadens) # rescale, so sum(weights) = 1
				
				fitted[s,t] <- fitted[s,t] +
					+ sum(gammadens*fitted[datalist$neighlist[[s]][[ss]],
																 t-(0:datalist$maxlag)]) # add
				# ^ add to predmat or replace it completely
				# ^ weighted comb of pred from ustr neighbors at lag 0:maxlag
			}
		}
	}
	
	
	
	#----------------------------------------------------------------------------#
	# pnll eval at fitted
	#----------------------------------------------------------------------------#
	
	pnll <- 0 # ini pen neg loglik
	
	for (s in 1:nS){ # loop over all loc
		for (t in (1+datalist$maxlag):nT){ # time points after burn-in
			# for (t in 1:nT){ # blockwise after block 1: nT=maxlag, no burn-in here
			# time = outer loop, space = inner loop
			if (datalist$obsindmat[s,t]){ # lkhd contrib only when obs available (tr set)
				pnll <- pnll - dnorm(x=datalist$obsmat[s,t], mean=fitted[s,t],
														 sd=sigma, log=T)
			}
		}
	}
	
	
	#----------------------------------------------------------------------------#
	# Outputs
	#----------------------------------------------------------------------------#
	
	REPORT(sigma)
	REPORT(wscale)
	REPORT(wshapebeta)
	REPORT(fitted) # fitted values on modeling scale
	
	return(pnll)
}

rout_nll_block <- function(par){
	# * data in R global envir assumed to be stored in datalist
	# * par vector order:
	#   - log_sigma (1)
	#   - log_wscale (1)
	#   - wshapebeta (p)
	#   - b (1), index for blocks in DataEval
	
	#----------------------------------------------------------------------------#
	# Inputs
	#----------------------------------------------------------------------------#
	
	# parms <- relist(par)
	
	# getAll(par, datalist)
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
	
	# sigma <- exp(parvec[1])
	# wscale <- exp(parvec[2]) # cst, gamma shape has all the cov
	# wshapebeta <- parvec[3:length(parvec)]
	sigma <- exp(par$log_sigma)
	wscale <- exp(par$log_wscale) # cst, gamma shape has all the cov
	wshapebeta <- par$wshapebeta
	
	nS <- nrow(datalist$obsmat) # nb loc overall (polyg and stations)
	ncolpredmatprev <- ncol(datalist$predmatprev)
	
	# predmat1 <- DataEval(f=function(b){
	# 	bvec <- 1:datalist$maxlag + (b-2)*datalist$maxlag + 2*datalist$maxlag # subset time points (cols)
	# 	return(as.numeric(cbind(
	# 		datalist$predmatprev[,(ncol(datalist$predmatprev)-datalist$maxlag+1):ncol(datalist$predmatprev)],
	# 		datalist$predmat[,bvec]
	# 	)))
	# }, x=par$b)
	# obsmat1 <- DataEval(f=function(b){
	# 	bvec <- 1:datalist$maxlag + (b-2)*datalist$maxlag + 2*datalist$maxlag # subset time points (cols)
	# 	return(as.numeric(datalist$obsmat[,bvec]))
	# }, x=par$b)
	# obsindmat1 <- DataEval(f=function(b){
	# 	bvec <- 1:datalist$maxlag + (b-2)*datalist$maxlag + 2*datalist$maxlag # subset time points (cols)
	# 	return(as.numeric(datalist$obsindmat[,bvec]))
	# }, x=par$b)
	# # ^ TODO: need exception for last block if nT is not a multiple of maxlag
	# 
	# predmat1 <- matrix(predmat1, nS)
	# obsmat1 <- matrix(obsmat1, nS)
	# obsindmat1 <- matrix(obsindmat1, nS)
	
	
	# bvec <- DataEval(function(i){1:datalist$maxlag +
	# 		+ (i-2)*datalist$maxlag + # as.integer(par$b)
	# 		+ 2*datalist$maxlag
	# 	}, x=par$b)
	
	# begin good
	predmat1 <- DataEval(f=function(i){
		bvec <- 1:datalist$maxlag +
			+ (i-2)*datalist$maxlag +
			+ 2*datalist$maxlag
		return(cbind(
			datalist$predmatprev[,(ncolpredmatprev-datalist$maxlag+1):ncolpredmatprev],
			datalist$predmat[,bvec]))
	},x=par$b)
	obsmat1 <- DataEval(f=function(i){
		bvec <- 1:datalist$maxlag +
			+ (i-2)*datalist$maxlag +
			+ 2*datalist$maxlag
		return(datalist$obsmat[,bvec])
	},x=par$b)
	obsindmat1 <- DataEval(f=function(i){
		bvec <- 1:datalist$maxlag +
			+ (i-2)*datalist$maxlag +
			+ 2*datalist$maxlag
		return(datalist$obsindmat[,bvec])
	},x=par$b)
	# ^ TODO: need exception for last block if nT is not a multiple of maxlag
	# end good
	
	# # b <- 3 # for tests, hardcoding to check correct obj$fn() value
	# b <- datalist$b
	# # b <- as.integer(par$b)
	# bvec <- 1:datalist$maxlag + (b-2)*datalist$maxlag + # as.integer(par$b)
	# 	+ 2*datalist$maxlag # subset time points (cols)
	# #
	# predmat1 <- cbind(
	# 	datalist$predmatprev[,(ncolpredmatprev-datalist$maxlag+1):ncolpredmatprev],
	# 	datalist$predmat[,bvec]
	# )
	# obsmat1 <- datalist$obsmat[,bvec]
	# obsindmat1 <- datalist$obsindmat[,bvec]
	# # ^ TODO: need exception for last block if nT is not a multiple of maxlag
	
	
	nT <- ncol(predmat1) # nb total time steps
	# ^ nb col = 2*maxlag
	
	
	#----------------------------------------------------------------------------#
	# routing
	#----------------------------------------------------------------------------#
	
	fitted <- predmat1
	
	# for (t in 1:nT){ # blockwise after block 1: nT=maxlag, keep first maxlag*nS obs
	for (t in (1+datalist$maxlag):nT){
		for (s in datalist$routingorder){
			# loop over all loc, excl the ones most ustr where fitted[s,]=predmat[s,]
			for (ss in 1:length(datalist$neighlist[[s]])){ # direct ustr neighbors
				whshape_ss <- exp(wshapebeta%*%wshapecovlist[[s]][[ss]])
				# ^ lin comb (p->1) with log link for shape>0
				gammadens <- dgamma(x=c(datalist$lag0,# lag0 = same day = e.g. 1e-3
																1:datalist$maxlag),
														shape=whshape_ss,
														scale=wscale)
				gammadens <- gammadens/sum(gammadens) # rescale, so sum(weights) = 1
				
				fitted[s,t] <- fitted[s,t] +
					+ sum(gammadens*fitted[datalist$neighlist[[s]][[ss]],
																 t-(0:datalist$maxlag)]) # add
				# ^ add to predmat or replace it completely
				# ^ weighted comb of pred from ustr neighbors at lag 0:maxlag
			}
		}
	}
	fitted <- fitted[,(1+datalist$maxlag):nT] # remove previous block
	
	
	
	#----------------------------------------------------------------------------#
	# pnll eval at fitted
	#----------------------------------------------------------------------------#
	
	# pnll <- 0 # ini pen neg loglik
	# 
	# for (s in 1:nS){ # loop over all loc
	# 	# for (t in 1:nT){ # blockwise after block 1: nT=maxlag, no burn-in here
	# 	# for (t in (1+data$maxlag):nT){ # time points after burn-in
	# 	for (t in 1:datalist$maxlag){ # time points after burn-in
	# 		# time = outer loop, space = inner loop
	# 		# if (obsindmat1[s,t]){ # lkhd contrib only when obs available (tr set)
	# 		if (obsindmat1[s,t]==1){ # lkhd contrib only when obs available (tr set)
	# 			# if (!is.na(obsmat1[s,t])){
	# 			pnll <- pnll - dnorm(x=obsmat1[s,t], mean=fitted[s,t],
	# 													 sd=sigma, log=T)
	# 		}
	# 	}
	# }
	
	# # obsind <- which(!is.na(obsmat1))
	# # pnll <- sum((obsmat1[obsind]-fitted[obsind])^2) # cannot have na.rm=T
	# obsind <- which(!is.na(obsmat1))
	# pnll <- sum((obsmat1[obsind]-fitted[obsind])^2) # cannot have na.rm=T
	
	pnll <- sum((obsmat1-fitted)^2*obsindmat1)
	# ^ NAs in obsmat replaced by arbitrary numeric (e.g. zero) but then
	#  multiplied by 0 in sum so ok
	
	
	#----------------------------------------------------------------------------#
	# Outputs
	#----------------------------------------------------------------------------#
	
	REPORT(sigma)
	REPORT(wscale)
	REPORT(wshapebeta)
	REPORT(par$b)
	REPORT(fitted) # fitted values on modeling scale
	
	# REPORT(obsind)
	# REPORT(obsindmat1)
	# REPORT(obsmat1)
	
	return(pnll)
}
# end routnll_blockwise
