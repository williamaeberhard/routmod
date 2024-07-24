# routnll_blockwise_eval: eval fn and gr of rounll_blockwise | v1.0
# * Change log:
#    - v1.0: different wshapecovlist for lake=0 and lake=1, including different
#      intercepts. Now supply distinct wshapebeta0 and wshapebeta1. Only applies
#      to dischargeinshape=0 (will see later if needed for dischargeinshape=1).
#    - v0.9.3: nothing changed here, just to match routnll_blockwise.r
#    - v0.9.2: nothing changed here, just to match routnll_blockwise.r
#    - v0.9:
#      - matching routnll_blockwise.r v0.9 with conditioning on lake status
#      - adapted to work with pred_routmod by calling only objects in datalist
#        in addition to parvec1
#    - v0.5: routnll_blockwise_ini and routnll_blockwise now output sum of
#      squared residuals as defualt loss, so that here we sum them and then
#      output $fn as an overall mean squared error
#    - v0.4: added intercept to wshapebeta (doesn't change anything below)
#    - v0.3: initial version

routnll_blockwise_eval <- function(parvec, predmat, maxlag, outputfitted=FALSE){
	
	### // setup ----
	nT <- ncol(predmat)
	if (!is.whole(nT/maxlag)){
		stop('nT must be a multiple of maxlag (for now)')
	}
	rpredmat <- predmat # routed predictions
	
	polyg_ref <- which(lapply(datalist$neighlist,'length')>1)[1]
	# ^ [1] arbitrary, just a polyg with >=1 neighb ustr, to extract cov dim
	p0 <- length(datalist$wshapecovlist0[[polyg_ref]][[1]]) + 1
	p1 <- length(datalist$wshapecovlist1[[polyg_ref]][[1]]) + 1
	# ^ cov dim for lake=0 and lake=1, incl intercept as in routnll_blockwise.r
	
	if (length(parvec) != (p0+p1+1)){ # +1 scale
		stop('length(parvec) != p0+p1+1, where p0 and p1 are from wshapecovlist0 ',
				 'and wshapecovlist1 (+1 for each intercept) in datalist.')
	}

	# wallclock <- proc.time()[3]
	
	
	### // ini ----
	# ini: first maxlag time points on which we condition, then first block of
	# maxlag time points for first contrib to loss function and fitted values
	
	parlist_ini <- list(
		'log_wscale'=parvec[1],
		# 'wshapebeta'=parvec[2:length(parvec)] # v0.4: now incl intercept as [1]
		'wshapebeta0'=parvec[2:(p0+1)],        # v1.0
		'wshapebeta1'=parvec[(p0+2):(p0+p1+1)] # v1.0
	)
	# lenbeta <- length(parlist_ini$wshapebeta)
	lenbeta <- length(parvec) # v1.0
	
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
	
	# print(str(obj_ini,1))
	
	objfn <- obj_ini$fn(unlist(parlist_ini))             # ini
	objgr <- as.numeric(obj_ini$gr(unlist(parlist_ini))) # ini

	# print(objfn)
	# print(str(objgr,1))
	
	# print(eval(expr=obj_ini$fn(unlist(parlist_ini)),envir=obj_ini$env))
	
	# system.time(
	rep_ini <- obj_ini$rep(unlist(parlist_ini))
	# )
	# ^ 0.3s | MBP13 Toy4 nS=9445 maxlag=2
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
		# 'wshapebeta'=parlist_ini$wshapebeta, # parvec[2],
		'wshapebeta0'=parlist_ini$wshapebeta0, # v1.0
		'wshapebeta1'=parlist_ini$wshapebeta1, # v1.0
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
	# ^ 11-13 s | MBP13 Toy4 nS=9445 maxlag=2
	# ^  | MBP13 Toy4 nS=9445 maxlag=5
	
	objfn <- objfn + obj_b$fn(unlist(parlist))
	# objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)])
	objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[1:lenbeta]) # v1.0
	# ^ v0.5: sum of squared resid
	
	# system.time(
	rep_b <- obj_b$rep(unlist(parlist))
	# )
	# ^ 0.3s | MBP13 Toy4 nS=9445 maxlag=2
	
	# head(rpredmat[,1:datalist$maxlag+(b-2)*datalist$maxlag+2*datalist$maxlag])
	# head(rep_b$fitted)
	
	# summary(predmatprev)
	# summary(rep_b$predmatprev1)
	# summary(rep_ini$fitted)
	# # ^ all equal, ok
	
	# update routed pred
	predmatprev <- rep_b$fitted
	rpredmat[, 1:maxlag+(b-2)*maxlag+2*maxlag] <- predmatprev
	# ^ not update first maxlag time points on which we condition, they remain equal
	#   to the non-routed predmat (supplied initial predictions)
	
	# predmatprev1_prevb <- rep_b$predmatprev1
	
	
	
	### // iterations over blocks, b>=3 ----
	
	# wallclock <- proc.time()[3]
	# if (is.whole(nT/maxlag)){
	# nT is a multiple of maxlag, last block is complete
	maxb <- nT/maxlag - 1L
	
	# b <- 3L
	for (b in 3:maxb){
		parlist <- list(
			'log_wscale'=parlist_ini$log_wscale, # parvec[1],
			# 'wshapebeta'=parlist_ini$wshapebeta, # parvec[2],
			'wshapebeta0'=parlist_ini$wshapebeta0, # v1.0
			'wshapebeta1'=parlist_ini$wshapebeta1, # v1.0
			'b'=b,
			'predmatprev'=predmatprev # update fitted from previous block
		)
		
		objfn <- objfn + obj_b$fn(unlist(parlist))
		# objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)])
		objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[1:lenbeta]) # v1.0
		# ^ v0.5: sum of squared resid
		
		# system.time(
		rep_b <- obj_b$rep(unlist(parlist))
		# )
		# ^ 0.3s | MBP13 Toy4 nS=9445 maxlag=2
		
		# summary(predmatprev)
		# summary(rep_b$predmatprev1)
		# summary(predmatprev1_prevb)
		# # ^ ok, predmatprev1 got updated
		
		# update routed pred
		predmatprev <- rep_b$fitted
		rpredmat[, 1:maxlag+(b-2)*maxlag+2*maxlag] <- predmatprev
		# ^ not update first maxlag time points on which we condition, they remain equal
		#   to the non-routed predmat (supplied initial predictions)
		
		# predmatprev1_prevb <- rep_b$predmatprev1
	} # for b in 3:maxb
	# } else {
	# nT is not a multiple of maxlag, last block is not complete, need to call
	# rout_nll_block_last
	# 
	# stop('nT must be a multiple of maxlag (for now)')
	# 
	# maxb <- floor(nT/maxlag) #
	# 
	# for (b in 3:(maxb-1)){
	# 	parlist <- list(
	# 		'log_wscale'=parlist_ini$log_wscale, # parvec[1],
	# 		'wshapebeta'=parlist_ini$wshapebeta, # parvec[2],
	# 		'b'=b,
	# 		'predmatprev'=predmatprev # update fitted from previous block
	# 	)
	# 
	# 	objfn <- objfn + obj_b$fn(unlist(parlist))
	# 	objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)])
	# 
	# 	# system.time(
	# 	rep_b <- obj_b$rep(unlist(parlist))
	# 	# )
	# 	# ^ 0.3s | MBP13 Toy4 nS=9445 maxlag=2
	# 
	# 	# update routed pred
	# 	predmatprev <- rep_b$fitted
	# 	rpredmat[,1:datalist$maxlag+(b-2)*datalist$maxlag+2*datalist$maxlag] <- predmatprev
	# 	# ^ not update first maxlag time points on which we condition, they remain equal
	# 	#   to the non-routed predmat (supplied initial predictions)
	# 
	# } # for b in 3:(maxb-1)
	# 
	# 
	# #### last (incomplete block with b=maxb)
	# bvec <- 1:datalist$maxlag+(maxb-2)*datalist$maxlag+2*datalist$maxlag
	# bvec <- bvec[-which(bvec>nT)]
	# 
	# parlist_last <- parlist_ini
	# 
	# # datalist_last <- list(
	# # 	'obsmat'=as.matrix(obsmat01[, bvec]),
	# # 	'obsindmat'=as.matrix(obsindmat01[, bvec]),
	# # 	'predmat'=cbind(predmatprev, rpredmat[, bvec]),
	# # 	'routingorder'=routingorder,
	# # 	'neighlist'=neighlist,
	# # 	'wshapecovlist'=wshapecovlist, # alternative with additional sqrt
	# # 	'maxlag'=maxlag,
	# # 	'lag0'=lag0
	# # )
	# 
	# # system.time(
	# obj_last <- MakeADFun(
	# 	func = rout_nll_block_last,
	# 	parameters = parlist_last,
	# 	silent = T
	# )
	# # )
	# 
	# objfn <- objfn + obj_last$fn(unlist(parlist_last))
	# objgr <- objgr + as.numeric(obj_last$gr(unlist(parlist_last)))
	# 
	# # system.time(
	# rep_last <- obj_last$rep(unlist(parlist_last))
	# # )
	# 
	# # update routed pred
	# rpredmat[,bvec] <- rep_last$fitted
	# }
	# print(proc.time()[3]-wallclock)
	
	
	### // output ----
	n.active <- sum(datalist$obsindmat[,(maxlag+1):nT])
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
# end routnll_blockwise_eval
