# routnll_blockwise_firsteval_nolake: eval fn and gr of rounll_blockwise | v0.9
# * Change log:
#    - v0.9: initial version, forked from routnll_blockwise_firsteval v0.9


### // setup ----
nT <- ncol(predmat)
if (!is.whole(nT/maxlag)){
	stop('nT must be a multiple of maxlag (for now)')
}
rpredmat <- predmat # routed predictions

parvec <- parvec1
# ^ v0.9: c(log_wscale, intercept, coef for cov)


### // ini ----
# ini: first maxlag time points on which we condition, then first block of
# maxlag time points for first contrib to loss function and fitted values

parlist_ini <- list(
	'log_wscale'=parvec[1],
	'wshapebeta'=parvec[2:length(parvec)] # v0.4: now incl intercept as [1]
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
	func = rout_nll_block_nolake_ini,
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
suppressWarnings(
	obj_b <- MakeADFun(
		func = rout_nll_block_nolake,
		parameters = parlist,
		silent = T
	)
) # 
# )

objfn <- objfn + obj_b$fn(unlist(parlist)) # sum of squared resid
objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)])

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
rpredmat[,1:datalist$maxlag+(b-2)*datalist$maxlag+2*datalist$maxlag] <- predmatprev
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
		'wshapebeta'=parlist_ini$wshapebeta, # parvec[2],
		'b'=b,
		'predmatprev'=predmatprev # update fitted from previous block
	)
	
	objfn <- objfn + obj_b$fn(unlist(parlist)) # sum of squared resid
	objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)])
	
	# system.time(
	rep_b <- obj_b$rep(unlist(parlist))
	# )
	
	# summary(predmatprev)
	# summary(rep_b$predmatprev1)
	# summary(predmatprev1_prevb)
	# # ^ ok, predmatprev1 got updated
	
	# update routed pred
	predmatprev <- rep_b$fitted
	rpredmat[,1:datalist$maxlag+(b-2)*datalist$maxlag+2*datalist$maxlag] <- predmatprev
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

n.active <- sum(datalist$obsindmat[,(datalist$maxlag+1):nT])
objfn <- objfn/n.active
objgr <- objgr/n.active
# ^ sum of squared resid => MSE

fit1 <- list('fn'=objfn,'gr'=objgr)
# }

# end routnll_blockwise_firsteval_nolake
