routmod: R code for routing module to improve river discharge predictions along a river network
-----------------------------------------------------------------------------------------------

### Todo

* [ ] 


### Version History

This is routmod version 0.5.

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


