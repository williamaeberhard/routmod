# routnll blockwise: neg log lik for routing | v0.3
# * Change log:
#    - v0.3: clean up notation
#    - v0.2: initial version, forked from routmod/rout_nll v0.2

# * Convention for vector stacking: time = outer loop, space = inner loop

rout_nll_block_ini <- function(par){
	# * par vector order:
	#   - log_wscale (1)
	#   - wshapebeta (p)
	# * datalist_ini must include:
	#   - obsmat: numeric matrix of discharge (m3/s) starting after first maxlag
	#     time points, no NAs, nS x maxlag
	#   - obsindmat: 0-1 matrix, same dim as obsmat (to multiply elementwise)
	#   - predmat: numeric matrix of pred discharge starting at first time point,
	#     no NAs, nS x (2*maxlag)
	#   - maxlag: integer >=1, max lag in routing from spatial neighbors
	#   - routingorder: integer vector routing ustr -> dstr, within 1:nS but of
	#     length <=nS since excl loc most ustr where fitted=predmat
	#   - neighlist: list of length nS, each element being an integer vector of
	#     values within 1:nS of direct/Markov neighbors ustr, not incl itself, so
	#     each entry is of length within {0, ..., nS-1}
	#   - wshapecovlist: list of length nS, same ordering as neighlist, each
	#     element being a list of numeric vectors (number of vectors is length of
	#     the corresponding element in neighlist), each vector being of same
	#     length p = length(wshapebeta) where the values are static cov
	#   - lag0: lag 0 = same day for gamma kernel, but >0, e.g. 1e-3
	
	
	#----------------------------------------------------------------------------#
	# Inputs
	#----------------------------------------------------------------------------#
	
	# getAll(par, data)
	# getAll(par)
	
	
	#----------------------------------------------------------------------------#
	# Setup and init
	#----------------------------------------------------------------------------#
	
	wscale <- exp(par$log_wscale) # cst, gamma shape has all the cov
	wshapebeta <- par$wshapebeta
	
	nS <- nrow(datalist_ini$obsmat) # nb loc overall (polyg and stations)
	nT <- ncol(datalist_ini$predmat) # not the total nb of time steps
	# ^ ini: ncols(predmat)=2*maxlag, but ncols(obsmat) = latter maxlag
	
	
	#----------------------------------------------------------------------------#
	# routing
	#----------------------------------------------------------------------------#
	
	fitted <- datalist_ini$predmat
	for (t in (1+datalist_ini$maxlag):nT){ # ini: cond on 1st maxlag obs
		for (s in datalist_ini$routingorder){
			# loop over all loc excl the ones most ustr where fitted[s,]=predmat[s,]
			for (ss in 1:length(datalist_ini$neighlist[[s]])){ # direct ustr neighbors
				whshape_ss <- exp(wshapebeta%*%datalist_ini$wshapecovlist[[s]][[ss]])
				# ^ lin comb (p->1) with log link for shape>0
				gammadens <- dgamma(
					x=c(datalist_ini$lag0, 1:datalist_ini$maxlag), # lag0 = same day
					shape=whshape_ss, scale=wscale
				)
				gammadens <- gammadens/sum(gammadens) # rescale so sum(weights) = 1
				fitted[s,t] <- fitted[s,t] +
					+ sum(gammadens*fitted[datalist_ini$neighlist[[s]][[ss]],
																 t-(0:datalist_ini$maxlag)]) # add to pred
				# ^ weighted comb of pred from ustr neighbors at lag 0:maxlag
			}
		} # for s in routingorder
	} # for t in (maxlag+1):nT
	fitted <- fitted[, (1+datalist_ini$maxlag):nT] # remove previous block
	# ^ output fitted is nS x maxlag, same dim as obsmat
	
	
	#----------------------------------------------------------------------------#
	# pnll eval at fitted
	#----------------------------------------------------------------------------#
	
	# pnll <- 0 # ini pen neg loglik
	# for (s in 1:nS){ # loop over all loc
	# 	for (t in (1+datalist$maxlag):nT){ # time points after burn-in
	# 		# time = outer loop, space = inner loop
	# 		if (datalist$obsindmat[s,t]){ # lkhd contrib only when obs available (tr set)
	# 			pnll <- pnll - dnorm(x=datalist$obsmat[s,t], mean=fitted[s,t],
	# 													 sd=sigma, log=T)
	# 		}
	# 	}
	# }
	
	# pnll <- sum((datalist_ini$obsmat-fitted)^2*datalist_ini$obsindmat)
	# # ^ sum of squared residuals
	pnll <- sum((datalist_ini$obsmat-fitted)^2*datalist_ini$obsindmat)/
		sum(datalist_ini$obsindmat) # mean squared loss
	# ^ NAs in obsmat replaced by arbitrary numeric (e.g. zero) but then
	#   multiplied by 0 in sum so no contribution to loss function
	
	
	#----------------------------------------------------------------------------#
	# Outputs
	#----------------------------------------------------------------------------#
	
	# REPORT(wscale)
	# REPORT(wshapebeta)
	REPORT(fitted) # fitted values on modeling scale
	
	return(pnll)
}

rout_nll_block <- function(par){
	# * par vector order:
	#   - log_wscale (1)
	#   - wshapebeta (p)
	#   - b (1), index for blocks in DataEval
	# * datalist must include:
	#   - obsmat: numeric matrix of discharge (m3/s), no NAs, nS x nT
	#   - obsindmat: 0-1 matrix, same dim as obsmat (to multiply elementwise)
	#   - predmat: numeric matrix of pred discharge, no NAs, same dim as obsmat
	#   - maxlag: integer >=1, max lag in routing from spatial neighbors
	#   - routingorder: integer vector routing ustr -> dstr, within 1:nS but of
	#     length <=nS since excl loc most ustr where fitted=predmat
	#   - neighlist: list of length nS, each element being an integer vector of
	#     values within 1:nS of direct/Markov neighbors ustr, not incl itself, so
	#     each entry is of length within {0, ..., nS-1}
	#   - wshapecovlist: list of length nS, same ordering as neighlist, each
	#     element being a list of numeric vectors (number of vectors is length of
	#     the corresponding element in neighlist), each vector being of same
	#     length p = length(wshapebeta) where the values are static cov
	#   - lag0: lag 0 = same day for gamma kernel, but >0, e.g. 1e-3
	
	#----------------------------------------------------------------------------#
	# Inputs
	#----------------------------------------------------------------------------#
	
	# getAll(par, datalist)
	
	
	#----------------------------------------------------------------------------#
	# Setup and init
	#----------------------------------------------------------------------------#
	
	wscale <- exp(par$log_wscale) # cst, gamma shape has all the cov
	wshapebeta <- par$wshapebeta
	
	nS <- nrow(datalist$obsmat) # nb loc overall (polyg and stations)
	nT <- ncol(datalist$predmat) # total nb time points
	
	predmat1 <- DataEval(f=function(i){
		bvec <- (1-datalist$maxlag):datalist$maxlag +
			+ (i-2)*datalist$maxlag +
			+ 2*datalist$maxlag
			# ^ replaces former predmatprev, now previous block's fitted overwrite
			#   values in rpredmat at each iteration
		if (bvec[datalist$maxlag]>nT){
			bvec <- bvec[-which(bvec>nT)] # last block: cap at nT
		}
		return(datalist$predmat[,bvec])
	},x=par$b)
	obsmat1 <- DataEval(f=function(i){
		bvec <- 1:datalist$maxlag +
			+ (i-2)*datalist$maxlag +
			+ 2*datalist$maxlag
		if (bvec[datalist$maxlag]>nT){
			bvec <- bvec[-which(bvec>nT)] # last block: cap at nT
		}
		return(datalist$obsmat[,bvec])
	},x=par$b)
	obsindmat1 <- DataEval(f=function(i){
		bvec <- 1:datalist$maxlag +
			+ (i-2)*datalist$maxlag +
			+ 2*datalist$maxlag
		if (bvec[datalist$maxlag]>nT){
			bvec <- bvec[-which(bvec>nT)] # last block: cap at nT
		}
		return(datalist$obsindmat[,bvec])
	},x=par$b)
	
	# ### testing with hardcoding b, to check correct obj$fn() value
	# # b <- 2
	# bvec <- 1:datalist$maxlag + (b-2)*datalist$maxlag + # as.integer(par$b)
	# 	+ 2*datalist$maxlag # subset time points (cols)
	# #
	# predmat1 <- cbind(
	# 	datalist$predmatprev[,(ncolpredmatprev-datalist$maxlag+1):ncolpredmatprev],
	# 	datalist$predmat[,bvec]
	# )
	# obsmat1 <- datalist$obsmat[,bvec]
	# obsindmat1 <- datalist$obsindmat[,bvec]
	
	nT1 <- ncol(predmat1) # nb time points in this block
	# ^ nb col = 2*maxlag
	
	
	#----------------------------------------------------------------------------#
	# routing
	#----------------------------------------------------------------------------#
	
	fitted <- predmat1
	for (t in (1+datalist$maxlag):nT1){
		for (s in datalist$routingorder){
			# loop over all loc, excl the ones most ustr where fitted[s,]=predmat[s,]
			for (ss in 1:length(datalist$neighlist[[s]])){ # direct ustr neighbors
				whshape_ss <- exp(wshapebeta%*%wshapecovlist[[s]][[ss]])
				# ^ lin comb (p->1) with log link for shape>0
				gammadens <- dgamma(
					x=c(datalist$lag0, 1:datalist$maxlag), # lag0 = same day
					shape=whshape_ss, scale=wscale
				)
				gammadens <- gammadens/sum(gammadens) # rescale, so sum(weights) = 1
				fitted[s,t] <- fitted[s,t] +
					+ sum(gammadens*fitted[datalist$neighlist[[s]][[ss]],
																 t-(0:datalist$maxlag)]) # add to pred
				# ^ weighted comb of pred from ustr neighbors at lag 0:maxlag
			}
		} # for s in routingorder
	} # for t in (maxlag+1):nT1
	fitted <- fitted[, (1+datalist$maxlag):nT1]
	# ^ remove previous block, now maxlag cols, matching obsmat1
	
	
	#----------------------------------------------------------------------------#
	# pnll eval at fitted
	#----------------------------------------------------------------------------#
	
	# pnll <- sum((obsmat1-fitted)^2*obsindmat1) # sum of squared residuals
	pnll <- sum((obsmat1-fitted)^2*obsindmat1)/sum(obsindmat1) # mean squared loss
	# ^ NAs in obsmat replaced by arbitrary numeric (e.g. zero) but then
	#   multiplied by 0 in sum so no contribution to loss function
	
	
	#----------------------------------------------------------------------------#
	# Outputs
	#----------------------------------------------------------------------------#
	
	# REPORT(wscale)
	# REPORT(wshapebeta)
	# REPORT(par$b)
	REPORT(fitted) # fitted values on modeling scale
	
	return(pnll)
}
# end routnll_blockwise
