routmod: R code for routing module to improve river discharge predictions along a river network
-----------------------------------------------------------------------------------------------

### Todo

* [ ] adapt *_nolake vertsion of routnll_blockwise_firsteval.r and routnll_blockwise_eval.r to work with pred_routmod by calling only objects in datalist and parvec1
* [ ] adapt *_nocov vertsion of routnll_blockwise_firsteval.r and routnll_blockwise_eval.r to work with pred_routmod by calling only objects in datalist and parvec1
* [ ] adapt *_onlylake vertsion of routnll_blockwise_firsteval.r and routnll_blockwise_eval.r to work with pred_routmod by calling only objects in datalist and parvec1



### Version History

This is routmod version 0.9.

* v0.9:
  - routnll_blockwise.r: now conditioning on lake dummy variable, distinct sets of param in wshapebeta for lake=0 and lake=1.
  - created *_nolake versions of routnll_blockwise.r, routnll_blockwise_firsteval.r, and routnll_blockwise_eval.r, for routing without lake conditioning, i.e. equivalent to v0.8.
  - created *_nocov versions of routnll_blockwise.r, routnll_blockwise_firsteval.r, and routnll_blockwise_eval.r, for routing with same param everywhere, i.e., constant scale and shape parameters for the gamma density kernel.
 - created *_onlylake versions of routnll_blockwise.r, routnll_blockwise_firsteval.r, and routnll_blockwise_eval.r, for routing with lake conditioning but nother covariates, i.e., no wshapecovlist is used and there are only two distinct intercepts for lake=0 and lake=1 in wshapebeta.
 - adapted routnll_blockwise_firsteval.r and routnll_blockwise_eval.r to work with pred_routmod by calling only objects in datalist in addition to parvec1
* v0.8:
  - routnll_blockwise.r: added option datalist$losscode=3 for MSE on sqrt runoff scale. $fitted remains on raw discharge scale.
* v0.7:
  - routnll_blockwise.r: added option datalist$losscode=2 for MSE on runoff=discharge/area scale. as in v0.6, $fitted remains on raw discharge scale.
* v0.6:
  - added option for different loss functions in routnll_blockwise, through datalist_ini$losscode and datalist$losscode. For now only 0 (MSE on raw discharge scale) and 1 (MSE on sqrt discharge scale) are implemented. $fitted remains on raw discharge scale regardless of losscode, so that all post-processing is unaffected.
* v0.5:
  - fixed as.numeric coercion in rout_nll_block_ini and rout_nll_block functions.
  - changed default loss in routnll_blockwise_firsteval.r and routnll_blockwise_eval.r to be overall MSE (before was a sum of blockwise MSE, not great).
* v0.4:
  - added intercept to wshapebeta in routnll_blockwise.r. Made notes in routnll_blockwise_firsteval.r and routnll_blockwise_eval.r although nothing changed there.
* v0.3:
  - created routnll_blockwise.r which defines rout_nll_block_ini and rout_nll_block
  - created routnll_blockwise_eval.r which defines routnll_blockwise_eval which evaluates the fn and gr by block, assuming datalist_ini and datalist are defined in the global envir
  - created routnll_blockwise_eval1.r which is to be sourced once after defining datalist_ini and datalist are defined in the global envir to define obj_ini and obj_b, to then be called by routnll_blockwise_eval <= replaced by routnll_blockwise_firsteval.r
* v0.2: routnll.r added option in datalist for routingpred: "add", "replace", or "convcomb"
* v0.1: initial release.


