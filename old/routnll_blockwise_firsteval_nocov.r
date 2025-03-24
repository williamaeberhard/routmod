# routnll_blockwise_firsteval_nocov: eval fn and gr of rounll_blockwise | v0.9.2
# * Change log:
#    - v0.9.2:
#      - adapted to work with pred_routmod by calling only objects in datalist
#        in addition to parvec1
#      - should work with routnll_blockwise_nocov.r v0.9.2 with routed discharge
#        as cov in wshape.
#    - v0.9: initial version, forked from routnll_blockwise_firsteval_nolake v0.9


### // setup ----
# nT <- ncol(predmat) # bad, should call datalist
nT.eval <- ncol(datalist$predmat)
maxlag.eval <- as.integer(datalist$maxlag)

if (!is.whole(nT.eval/maxlag.eval)){
	stop('ncol(datalist$predmat) must be a multiple of datalist$maxlag.')
}

# rpredmat <- predmat # bad, should call datalist
rpredmat.eval <- datalist$predmat

parvec <- parvec1
# ^ v0.9 nocov: c(log_wscale, intercept=log_wshape)


### // ini ----
# ini: first maxlag time points on which we condition, then first block of
# maxlag time points for first contrib to loss function and fitted values

parlist_ini <- list(
	'log_wscale'=parvec[1],
	# 'wshapebeta'=parvec[2] # v0.9 nocov: only intercept
	'wshapebeta'=parvec[2:length(parvec)] # v0.9.2: better
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

obj_ini <- MakeADFun(
	func = rout_nll_block_nocov_ini,
	parameters = parlist_ini,
	silent = T
)

# print(str(obj_ini,1))

objfn <- obj_ini$fn(unlist(parlist_ini))             # ini, sum of squared resid
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
rpredmat.eval[,(maxlag.eval+1):(2*maxlag.eval)] <- predmatprev
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
suppressWarnings(
	obj_b <- MakeADFun(
		func = rout_nll_block_nocov,
		parameters = parlist,
		silent = T
	)
) # 
# )

objfn <- objfn + obj_b$fn(unlist(parlist)) # sum of squared resid
# objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,2)])
objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)]) # v0.9.2

# system.time(
rep_b <- obj_b$rep(unlist(parlist))
# )

# head(rpredmat[,1:datalist$maxlag+(b-2)*datalist$maxlag+2*datalist$maxlag])
# head(rep_b$fitted)

# summary(predmatprev)
# summary(rep_b$predmatprev1)
# summary(rep_ini$fitted)
# # ^ all equal, ok

# update routed pred
predmatprev <- rep_b$fitted
rpredmat.eval[,1:maxlag.eval+(b-2)*maxlag.eval+2*maxlag.eval] <- predmatprev
# ^ not update first maxlag time points on which we condition, they remain equal
#   to the non-routed predmat (supplied initial predictions)

# predmatprev1_prevb <- rep_b$predmatprev1



### // iterations over blocks, b>=3 ----

# wallclock <- proc.time()[3]
# if (is.whole(nT/maxlag)){
# nT is a multiple of maxlag, last block is complete
maxb <- nT.eval/maxlag.eval - 1L

# b <- 3L
for (b in 3:maxb){
	parlist <- list(
		'log_wscale'=parlist_ini$log_wscale, # parvec[1],
		'wshapebeta'=parlist_ini$wshapebeta, # parvec[2],
		'b'=b,
		'predmatprev'=predmatprev # update fitted from previous block
	)
	
	objfn <- objfn + obj_b$fn(unlist(parlist)) # sum of squared resid
	# objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,2)])
	objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)]) # v0.9.2
	
	# system.time(
	rep_b <- obj_b$rep(unlist(parlist))
	# )
	
	# summary(predmatprev)
	# summary(rep_b$predmatprev1)
	# summary(predmatprev1_prevb)
	# # ^ ok, predmatprev1 got updated
	
	# update routed pred
	predmatprev <- rep_b$fitted
	rpredmat.eval[,1:maxlag.eval+(b-2)*maxlag.eval+2*maxlag.eval] <- predmatprev
	# ^ not update first maxlag time points on which we condition, they remain equal
	#   to the non-routed predmat (supplied initial predictions)
	
	# predmatprev1_prevb <- rep_b$predmatprev1
} # for b in 3:maxb

# print(proc.time()[3]-wallclock)


### // output ----
# if (outputfitted){
# 	out <- list('fn'=objfn,'gr'=objgr,'fitted'=rpredmat)
# }	else {
# 	out <- list('fn'=objfn,'gr'=objgr)
# 	# better within optim
# }

# return(out)

n.active <- sum(datalist$obsindmat[,(maxlag.eval+1):nT.eval])
objfn <- objfn/n.active
objgr <- objgr/n.active
# ^ sum of squared resid => MSE

fit1 <- list('fn'=objfn,'gr'=objgr)
# }

# end routnll_blockwise_firsteval_nocov
