# routnll_blockwise_eval: eval fn and gr of rounll_blockwise | v0.4
# * Change log:
#    - v0.4: added intercept to wshapebeta (doesn't change anything below)
#    - v0.3: initial version

# routnll_blockwise_eval <- function(parvec,
# 																	 # obsmat01, obsindmat01,
# 																	 predmat,
# 																	 # routingorder, neighlist, wshapecovlist,
# 																	 # lag0=1e-3,
# 																	 maxlag, outputfitted=FALSE){

### // setup ----
nT <- ncol(predmat)
if (!is.whole(nT/maxlag)){
	stop('nT must be a multiple of maxlag (for now)')
}
rpredmat <- predmat # routed predictions

parvec <- parvec1
# ^ v0.4: now incl intercept as [2], log_wscale remains [1]


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

# assign(x='datalist_ini',value=datalist_ini,pos=1) # assign to calling envir

obj_ini <- MakeADFun(
	func = rout_nll_block_ini,
	parameters = parlist_ini,
	silent = T
)

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
		func = rout_nll_block,
		parameters = parlist,
		silent = T
	)
) # ^ imaginary parts discarded in as.numeric(par$predmatprev)
# )
# ^ 11-13 s | MBP13 Toy4 nS=9445 maxlag=2
# ^  | MBP13 Toy4 nS=9445 maxlag=5

objfn <- objfn + obj_b$fn(unlist(parlist))
objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)])

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
	
	objfn <- objfn + obj_b$fn(unlist(parlist))
	objgr <- objgr + as.numeric(obj_b$gr(unlist(parlist))[c(1,1:lenbeta+1)])
	
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
	rpredmat[,1:datalist$maxlag+(b-2)*datalist$maxlag+2*datalist$maxlag] <- predmatprev
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
# if (outputfitted){
# 	out <- list('fn'=objfn,'gr'=objgr,'fitted'=rpredmat)
# }	else {
# 	out <- list('fn'=objfn,'gr'=objgr)
# 	# better within optim
# }

# return(out)

fit1 <- list('fn'=objfn,'gr'=objgr)
# }
# end routnll_blockwise_eval
