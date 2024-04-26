# routnll_blockwise_eval_nocov: eval fn and gr of rounll_blockwise | v0.9
# * Change log:
#    - v0.9: initial version, forked from routnll_blockwise_eval_nolake v0.9

routnll_blockwise_eval_nocov <- function(parvec, predmat, maxlag, outputfitted=FALSE){
	
	### // setup ----
	nT <- ncol(predmat)
	if (!is.whole(nT/maxlag)){
		stop('nT must be a multiple of maxlag (for now)')
	}
	rpredmat <- predmat # routed predictions
	

	# wallclock <- proc.time()[3]
	
	
	### // ini ----
	# ini: first maxlag time points on which we condition, then first block of
	# maxlag time points for first contrib to loss function and fitted values
	
	parlist_ini <- list(
		'log_wscale'=parvec[1],
		# 'wshapebeta'=parvec[2:length(parvec)] #
		'wshapebeta'=parvec[2] #
	)
	lenbeta <- length(parlist_ini$wshapebeta)
	
	# tvec <- (maxlag+1):(2*maxlag)
	# 
	# datalist_ini <- list(
	# 	'obsmat'=obsmat01[, tvec],
	# 	'obsindmat'=obsindmat01[, tvec],
	# 	'predmat'=rpredmat[, 1:(2*maxlag)],
	# 	'routingorder'=routingorder,
	# 	'neighlist'=neighlist,
	# 	'wshapecovlist'=wshapecovlist,
	# 	'maxlag'=maxlag,
	# 	'lag0'=lag0,
	# 	'tvec'=tvec
	# )
	
	# assign(x='datalist_ini',value=datalist_ini,pos=1) # assign to calling envir
	
	# obj_ini <- MakeADFun(
	# 	func = rout_nll_block_ini,
	# 	parameters = parlist_ini,
	# 	silent = T
	# )
	
	
	objfn <- obj_ini$fn(unlist(parlist_ini))             # ini, sum squared resid
	objgr <- as.numeric(obj_ini$gr(unlist(parlist_ini))) # ini

	
	# system.time(
	rep_ini <- obj_ini$rep(unlist(parlist_ini))
	# )
	# str(rep_ini,1)
	
	# print(str(sdreport(obj_ini)))
	
	# update routed pred
	predmatprev <- rep_ini$fitted
	rpredmat[,(maxlag+1):(2*maxlag)] <- predmatprev
	# ^ not update first maxlag time points on which we condition, they remain equal
	#   to the non-routed predmat (supplied initial predictions)
	
	
	### // start iterations over blocks with b=2 ----
	b <- 2L
	
	parlist <- list(
		'log_wscale'=parlist_ini$log_wscale, # parvec[1],
		'wshapebeta'=parlist_ini$wshapebeta, # parvec[2],
		'b'=b,
		'predmatprev'=predmatprev # update fitted from previous block
	)
	
	# datalist <- list(
	# 	'obsmat'=obsmat01,
	# 	'obsindmat'=obsindmat01,
	# 	'predmat'=rpredmat,
	# 	'routingorder'=routingorder,
	# 	'neighlist'=neighlist,
	# 	'wshapecovlist'=wshapecovlist,
	# 	'maxlag'=maxlag,
	# 	'lag0'=lag0
	# )
	
	# system.time(
	# suppressWarnings(
	# 	obj_b <- MakeADFun(
	# 		func = rout_nll_block,
	# 		parameters = parlist,
	# 		silent = T
	# 	)
	# ) # ^ imaginary parts discarded in as.numeric(par$predmatprev)
	# )
	
	objfn <- objfn + obj_b$fn(unlist(parlist)) # sum of squared resid
	# objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)])
	objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,2)])

	
	# system.time(
	rep_b <- obj_b$rep(unlist(parlist))
	# )
	
	# update routed pred
	predmatprev <- rep_b$fitted
	rpredmat[,1:datalist$maxlag+(b-2)*datalist$maxlag+2*datalist$maxlag] <- predmatprev
	# ^ not update first maxlag time points on which we condition, they remain equal
	#   to the non-routed predmat (supplied initial predictions)
	
	# predmatprev1_prevb <- rep_b$predmatprev1
	
	
	
	### // iterations over blocks, b>=3 ----
	
	# wallclock <- proc.time()[3]

	maxb <- nT/maxlag - 1L
	
	# b <- 3L
	for (b in 3:maxb){
		parlist <- list(
			'log_wscale'=parlist_ini$log_wscale, # parvec[1],
			'wshapebeta'=parlist_ini$wshapebeta, # parvec[2],
			'b'=b,
			'predmatprev'=predmatprev # update fitted from previous block
		)
		
		objfn <- objfn + obj_b$fn(unlist(parlist)) # sum of squared resid
		# objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)])
		objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,2)])
		
		# system.time(
		rep_b <- obj_b$rep(unlist(parlist))
		# )
		
		# update routed pred
		predmatprev <- rep_b$fitted
		rpredmat[,1:datalist$maxlag+(b-2)*datalist$maxlag+2*datalist$maxlag] <- predmatprev
		# ^ not update first maxlag time points on which we condition, they remain equal
		#   to the non-routed predmat (supplied initial predictions)
		
		# predmatprev1_prevb <- rep_b$predmatprev1
	} # for b in 3:maxb
	
	# print(proc.time()[3]-wallclock)
	
	
	### // output ----
	n.active <- sum(datalist$obsindmat[,(datalist$maxlag+1):nT])
	objfn <- objfn/n.active
	objgr <- objgr/n.active
	# ^ sum of squared resid => MSE

	if (outputfitted){
		out <- list('fn'=objfn,'gr'=objgr,'fitted'=rpredmat)
	}	else {
		out <- list('fn'=objfn,'gr'=objgr)
		# better within optim
	}
	
	return(out)	
}
# end routnll_blockwise_eval_nocov
