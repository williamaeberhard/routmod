# routnll blockwise: neg log lik for routing | v0.7
# * Change log:
#    - v0.7: added datalist_ini$losscode==2 for runoff scale (mm/day), assuming
#      datalist_ini$areamat has same dim as datalist_ini$obsmat
#    - v0.6: added option for different loss function, MSE on sqrt discharge
#      scale if datalist_ini$losscode==1. But fitted remains on original (raw)
#      discharge scale
#    - v0.5:
#      * removed as.numeric coercion for par$predmatprev, was creating weird
#            differences between ob$fn() and REPORTed objects.
#      * default loss function now is sum of squared residuals, so that
#        routnll_blockwise_firsteval.r and routnll_blockwise.r now process
#        obj_ini$fn and obj_b$fn as a sum and then divide by sum(obsindmat01)
#        so that overall fn is the correct mean squared error.
#    - v0.4: added intercept to wshapebeta
#    - v0.3: clean up notation
#    - v0.2: initial version, forked from routmod/rout_nll v0.2

# * Convention for vector stacking: time = outer loop, space = inner loop

is.whole <- function(x,tol=.Machine$double.eps^0.5){abs(x-round(x))<tol}

rout_nll_block_ini <- function(par){
	# * par vector order:
	#   - log_wscale (1)
	#   - wshapebeta (p, now incl intercept)
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
	#     length p-1 = length(wshapebeta)-1 where the values are static cov
	#   - lag0: lag 0 = same day for gamma kernel, but >0, e.g. 1e-3
	
	
	#----------------------------------------------------------------------------#
	# Inputs
	#----------------------------------------------------------------------------#
	
	# getAll(par, data)
	# getAll(par)
	
	# datalist_ini <- parent.frame(n=2)[['datalist_ini']] # works
	# assign(x='datalist_ini',value=parent.frame(n=2)[['datalist_ini']],pos=1) #
	# datalist_ini <- get(x='datalist_ini',envir=parent.frame(n=2)) # pos=2
	# assign(x='vv',value=2,pos=1) # envir=parent.frame()
	
	
	#----------------------------------------------------------------------------#
	# Setup and init
	#----------------------------------------------------------------------------#

	wscale <- exp(par$log_wscale) # cst, gamma shape has all the cov
	wshapebeta <- par$wshapebeta # v0.4: now includes intercept as [1]

	nS <- nrow(datalist_ini$obsmat) # nb loc overall (polyg and stations)
	nT <- ncol(datalist_ini$predmat) # not the total nb of time steps
	# ^ ini: ncol(predmat)=2*maxlag, but ncols(obsmat) = latter maxlag
	
	
	#----------------------------------------------------------------------------#
	# routing
	#----------------------------------------------------------------------------#

	fitted <- datalist_ini$predmat # v=0.6: remains on raw discharge scale
	# for (t in (1+datalist_ini$maxlag):nT){ # ini: cond on 1st maxlag obs
	for (t in datalist_ini$tvec){ # ini: cond on 1st maxlag obs
		for (s in datalist_ini$routingorder){
			# loop over all loc excl the ones most ustr where fitted[s,]=predmat[s,]
			for (ss in 1:length(datalist_ini$neighlist[[s]])){ # direct ustr neighbors
				# whshape_ss <- exp(wshapebeta%*%datalist_ini$wshapecovlist[[s]][[ss]])
				# whshape_ss <- exp(sum(wshapebeta*datalist_ini$wshapecovlist[[s]][[ss]]))
				whshape_ss <- exp(sum(wshapebeta[1] +
																wshapebeta[-1]*datalist_ini$wshapecovlist[[s]][[ss]]))
				# ^ v0.4: added intercept
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
	fitted <- fitted[, datalist_ini$tvec] # remove previous block
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

	if (datalist_ini$losscode==0){
		# 0 = SSE on raw discharge scale
		pnll <- sum((datalist_ini$obsmat-fitted)^2*datalist_ini$obsindmat)
	} else if (datalist_ini$losscode==1){
		# 1 = SSE on sqrt discharge scale, assuming obsmat is sqrt discharge
		pnll <- sum((datalist_ini$obsmat-sqrt(fitted))^2*datalist_ini$obsindmat)
	} else if (datalist_ini$losscode==2){
		# 2 = SSE on runoff=discharge/area*8.64e7 scale (mm/day), assuming obsmat is
		#     runoff and areamat is same dim as obsmat (rows are same area
		#     values, by rgs)
		pnll <- sum((datalist_ini$obsmat-fitted/datalist_ini$areamat)^2*
									datalist_ini$obsindmat)
		# ^ areamat already has 8.64e7 coef in it for mm/day
	}
	
	# pnll <- sum((datalist_ini$obsmat-fitted)^2*datalist_ini$obsindmat)/
	# 	sum(datalist_ini$obsindmat) # mean squared loss

	# ^ NAs in obsmat replaced by arbitrary numeric (e.g. zero) but then
	#   multiplied by 0 in sum so no contribution to loss function


	#----------------------------------------------------------------------------#
	# Outputs
	#----------------------------------------------------------------------------#
	
	
	# REPORT(wscale)
	# REPORT(wshapebeta)
	# REPORT(nT)
	# REPORT(datalist_ini$tvec)
	
	# # REPORT(ls(parent.frame(n=1)))
	# REPORT(parent.frame()[['p']])
	# REPORT(parent.frame()[['par']])
	
	# testout <- get(x='nT',pos=-1)
	# testout <- 1:get(x='nT',pos=1)
	# testout <- seq(1,get(x='nT'),1)
	# testout <- vv:1
	# testout <- seq(from=1, to=10, by=1)
	# REPORT(testout)
	
	# REPORT((1 + datalist_ini$maxlag):(nT))
	
	# REPORT(str(datalist_ini,1))
	# REPORT(maxl+1)
	
	# REPORT(sum(datalist_ini$obsmat))
	# REPORT(pnll)
	
	# pnll <- 12.3
	
	REPORT(fitted) # fitted values on modeling scale
	return(pnll)
}

rout_nll_block <- function(par){
	# * par vector order:
	#   - log_wscale (1)
	#   - wshapebeta (p, now incl intercept)
	#   - b (1), index for blocks in DataEval
	#   - predmatprev (nS x maxlag), fitted from previous block
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
	#     length p-1 = length(wshapebeta)-1 where the values are static cov
	#   - lag0: lag 0 = same day for gamma kernel, but >0, e.g. 1e-3
	
	#----------------------------------------------------------------------------#
	# Inputs
	#----------------------------------------------------------------------------#
	
	# getAll(par, datalist)
	# datalist <- parent.frame(n=2)[['datalist']] # 
	
	
	#----------------------------------------------------------------------------#
	# Setup and init
	#----------------------------------------------------------------------------#
	
	wscale <- exp(par$log_wscale) # cst, gamma shape has all the cov
	wshapebeta <- par$wshapebeta # v0.4: now includes intercept as [1]
	
	nS <- nrow(datalist$obsmat) # nb loc overall (polyg and stations)
	nT <- ncol(datalist$predmat) # total nb time points
	
	# predmatprev1 <- DataEval(f=function(i){
	# 	matrix(as.double(i),nS)
	# },x=par$predmatprev)
	# # ^ just to declare it as data

	# predmatprev1 <- matrix(as.numeric(par$predmatprev),nS) # <= bad!
	predmatprev1 <- par$predmatprev
	
	predmat11 <- DataEval(f=function(i){
		# bvec <- (1-datalist$maxlag):datalist$maxlag +
		# 	+ (i-2)*datalist$maxlag +
		# 	+ 2*datalist$maxlag
		bvec <- 1:datalist$maxlag +
			+ (i-2)*datalist$maxlag +
			+ 2*datalist$maxlag
		# if (bvec[datalist$maxlag]>nT){
		# 	bvec <- bvec[-which(bvec>nT)] # last block: cap at nT
		# }
		# return(cbind(predmatprev1, datalist$predmat[,bvec])) # par$predmatprev
		return(datalist$predmat[,bvec]) # par$predmatprev
	},x=par$b)
	predmat1 <- cbind(predmatprev1, predmat11)

	obsmat1 <- DataEval(f=function(i){
		bvec <- 1:datalist$maxlag +
			+ (i-2)*datalist$maxlag +
			+ 2*datalist$maxlag
		# if (bvec[datalist$maxlag]>nT){
		# 	bvec <- bvec[-which(bvec>nT)] # last block: cap at nT
		# }
		return(datalist$obsmat[,bvec])
	},x=par$b)

	obsindmat1 <- DataEval(f=function(i){
		bvec <- 1:datalist$maxlag +
			+ (i-2)*datalist$maxlag +
			+ 2*datalist$maxlag
		# if (bvec[datalist$maxlag]>nT){
		# 	bvec <- bvec[-which(bvec>nT)] # last block: cap at nT
		# }
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
				# whshape_ss <- exp(wshapebeta%*%wshapecovlist[[s]][[ss]]) # wrong! missing datalist$
				# whshape_ss <- exp(sum(wshapebeta*datalist$wshapecovlist[[s]][[ss]]))
				whshape_ss <- exp(sum(wshapebeta[1] +
																wshapebeta[-1]*datalist$wshapecovlist[[s]][[ss]]))
				# ^ v0.4: added intercept
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
	
	fitted <- fitted[, (1+datalist$maxlag):nT1] # breaks if truncated bvec
	# fitted <- fitted[, (nT1-ncol(obsmat1)+1):nT1]
	# ^ remove previous block, now nb cols matches obsmat1 (= maxlag unless
	#   truncated bvec)
	
	
	#----------------------------------------------------------------------------#
	# pnll eval at fitted
	#----------------------------------------------------------------------------#
	
	if (datalist$losscode==0){
		# 0 = SSE on raw discharge scale
		pnll <- sum((obsmat1-fitted)^2*obsindmat1)
	} else if (datalist$losscode==1){
		# 1 = SSE on sqrt discharge scale, assuming obsmat is sqrt discharge
		pnll <- sum((obsmat1-sqrt(fitted))^2*obsindmat1)
	} else if (datalist_ini$losscode==2){
		# 2 = SSE on runoff=dischareg/area scale, assuming obsmat is runoff and
		#     areamat is same dim as obsmat1 (rows are same area values, by rgs)
		pnll <- sum((obsmat1-fitted/datalist$areamat)^2*obsindmat1)
		# ^ areamat already has 8.64e7 coef in it for mm/day
	}
	
	# pnll <- sum((obsmat1-fitted)^2*obsindmat1)/sum(obsindmat1) # mean squared loss

	# pnll <- sum(fitted) # test
	# pnll <- fitted[4393,5] # test

	# ^ NAs in obsmat replaced by arbitrary numeric (e.g. zero) but then
	#   multiplied by 0 in sum so no contribution to loss function
	
	
	#----------------------------------------------------------------------------#
	# Outputs
	#----------------------------------------------------------------------------#
	
	# REPORT(wscale)
	# REPORT(wshapebeta)
	# REPORT(par$b)
	REPORT(fitted) # fitted values on modeling scale
	
	# REPORT(predmatprev1) # to check if prematprev passed correctly as param
	# REPORT(obsmat1)
	# # REPORT(obsindmat1)
	# REPORT(predmat1)

	return(pnll)
}
