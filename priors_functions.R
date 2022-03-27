#################################################################
#################################################################
#################################################################
## Functions to Borrow using different Priors: 
#################################################################
#################################################################

## for all input a dataset which is a dataframe with the following columns:
## Trt: 0 indicates control, 1, .... A indicates the treatment arms 1 to A
## Trial: 0 indicates RWD membership, 1 indicates RCT membership
## 

get_draws <- function(dat, prior.text)

horseshoe_prior <- function(dat){}

spike_slab_prior <- function(dat){}

mem_prior <- function(dat){}

mem_classic <- function(dat){}

