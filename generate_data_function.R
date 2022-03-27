#################################################################
#################################################################
#################################################################
## Functions to be used for sims for platform RCT
#################################################################
#################################################################

## generate data as follows: 
## generate all X values from MVN with 
###         N = n_{RCT}+n_{RWD} (i.e. all observations from the whole trial)
###         lambda = vector of treatment effects
###         gamma = intercept
###         psi = (to start) bias of control RWD 
###         drift = treatment drift across time (if relevant)



generate_data <- function(NRCT, NRWD, 
                          num_X = 5, A,
                          lambda, gamma, psi = 0, drift = 0){
  
  N <- NRCT+NRWD ## total N
  
  trt_inds <- c(paste0("TRT_", seq(0, A) ) ) ## treatment indicators
  
  library(MASS); library(tidyr); library(dplyr); library(boot)
  X <-  data.frame( mvrnorm(n = N, 
                mu = rep(0, num_X), 
                Sigma = diag(1, nrow = num_X, ncol = num_X) ) ) ## generate all X values
  
  
  X_P_Trial <- mutate(X, P_trial = inv.logit( 1*X1+0.2*X2-0.5*X3+X4+X5) ) ## get prob of being in the trial P(RCT|X)
  
  RCT <- slice_sample(X_P_Trial, n = NRCT,
                      weight_by = X_P_Trial$P_trial, replace = F) ## sample from RCT with P(RCT|X) probability
  
  RWD <- setdiff(X_P_Trial, RCT) ## the rest are the RWD
  
  RCT$Trt <- as.vector(sapply(1:(NRCT/(A+1) ), 
                              function(p) sample(trt_inds ) ) ) ## set the Trt indicator for the RCT
  
  RCT$S <- 1 ##Study is 1
  
  RWD$Trt <- c(paste0("TRT_", seq(0) ) ) ## only borrowing from controls so all RWD is untreated
  
  RWD$S <- 0 ##Study is 0 for all RWD
  
  all_dat <- rbind(RCT, RWD) ## all data together
  
  X_dat <- model.matrix(as.formula( paste0("~-1+Trt+S+", 
                                           paste0(paste0("X", seq(1,num_X) ), collapse = "+")) ), 
                                    data = all_dat) 
  
  beta <-  c(lambda, gamma, psi, rep(0, num_X)) ## 0 beta for the X values, lambda: ctrl; gamma: trt; psi: bias
  
  # print( dim(X_dat) ); 
  # print( length(beta) )
  
  all_dat$prob <- inv.logit( X_dat%*%beta )
  
  all_dat$Y <- rbinom(N, 1, prob = all_dat$prob )
  
  ret.list <- list(all_dat = all_dat, 
                   beta = beta, 
                   X_dat = X_dat)
        
  
  return(ret.list)
  
}

test <- generate_data(5*100, 5*1000, 5, 4, 
                      lambda = log( 0.3/0.7) , 
                      gamma = c(-1, 1, 0.50, 
                                0.70), 
                      psi = 0.8, drift = 0)

## Check on things: 
#################################################################
#################################################################
all_dat_check = test$all_dat
print( paste0( "TRT == 0, both: ", mean(filter(all_dat_check, Trt == "TRT_0")$Y) ))
print( paste0( "TRT == 0, S == 1: ", mean(filter(all_dat_check, Trt == "TRT_0", S == 1)$Y) ) )
print( paste0( "TRT == 0, S == 0: ", mean(filter(all_dat_check, Trt == "TRT_0", S == 0)$Y) ) )
print( paste0( "TRT == 1: ", mean(filter(all_dat_check, Trt == "TRT_1")$Y) ) )
print( paste0( "TRT == 2: ", mean(filter(all_dat_check, Trt == "TRT_2")$Y) ) )
print( paste0( "TRT == 3: ", mean(filter(all_dat_check, Trt == "TRT_3")$Y) ) )
print( paste0( "TRT == 4: ", mean(filter(all_dat_check, Trt == "TRT_4")$Y) )  )
print( paste0( "TRT == 0, prob: ", mean(filter(all_dat_check, Trt == "TRT_0")$prob) ))
print( paste0( "TRT == 0, S == 1, prob: ", mean(filter(all_dat_check, Trt == "TRT_0", S == 1)$prob) ) )
print( paste0( "TRT == 0, S == 0, prob: ", mean(filter(all_dat_check, Trt == "TRT_0", S == 0)$prob) ) )
print( paste0( "TRT == 1, prob: ", mean(filter(all_dat_check, Trt == "TRT_1")$prob) ) )
print( paste0( "TRT == 2, prob: ", mean(filter(all_dat_check, Trt == "TRT_2")$prob) ) )
print( paste0( "TRT == 3, prob: ", mean(filter(all_dat_check, Trt == "TRT_3")$prob) ) )
print( paste0( "TRT == 4, prob: ", mean(filter(all_dat_check, Trt == "TRT_4")$prob) )  )

# print(head(X_dat) )
print(test$beta)