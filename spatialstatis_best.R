#################
#Packages
#################
library(mnormt)
library(msm)
library(MCMCpack)

##################################
#read in required data
##################################
y <- readRDS("D:/re-emergent RSV timing/sample_requiredata/y.rds")
x <- readRDS("D:/re-emergent RSV timing/sample_requiredata/x.rds")
neighbors_mat <- readRDS("D:/re-emergent RSV timing/sample_requiredata/neighbors_mat.rds")
censored <- readRDS("D:/re-emergent RSV timing/sample_requiredata/censored.rds")
censoring_time <- readRDS("D:/re-emergent RSV timing/sample_requiredata/censoring_time.rds") 

#############################################################
#Function
#############################################################
rho_fun<-function(rho_trans_val,
                  CAR_corr_mat){
  
  log_h<-determinant(CAR_corr_mat,
                     logarithm = TRUE)$modulus/2.00 +
    -0.50*crossprod(b, (CAR_corr_mat%*%b))/tau2 +
    rho_trans_val +
    -2.00*log(1.00 + exp(rho_trans_val))
  return(log_h)
  
}

####################
#Global Settings
####################
n<-length(y)
p_x<-ncol(x)

mcmc_samples<-200000

##############################################################
#Parameters + Initial Values
##############################################################
z<-rep(0.00,
       times = n)
eta<-rep(0.00,
         times = p_x)
b<-rep(0.00,
       times = n)
sigma2_epsilon<-1.00
tau2<-1.00
rho<-0.50
CAR_corr<-rho*(diag(rowSums(neighbors_mat)) - neighbors_mat) +
  (1.00 - rho)*diag(n)

eta_keep<-matrix(0.00,
                 nrow = mcmc_samples,
                 ncol = p_x)
b_keep<-matrix(0.00,
               nrow = mcmc_samples,
               ncol = n)
sigma2_epsilon_keep<-rep(0.00,
                         times = mcmc_samples)
sigma2_epsilon_keep[1]<-sigma2_epsilon
tau2_keep<-rep(0.00,
               times = mcmc_samples)
tau2_keep[1]<-tau2
rho_keep<-rep(0.00,
              times = mcmc_samples)
rho_keep[1]<-rho

####################
#Metropolis Settings
####################
acctot_rho<-1
rho_metrop_sd<-1.25

#########################
#Main Sampling Loop
#########################
for(i in 2:mcmc_samples){
  
  ###########################################################
  #z
  ###########################################################
  z[censored == 0]<-y[censored == 0]
  z[censored == 1]<-rtnorm(n = sum(censored == 1),
                           mean = ((x%*%eta)[censored == 1] +
                                     b[censored == 1]),
                           sd = sqrt(sigma2_epsilon),
                           lower = censoring_time,
                           upper = Inf)
  
  ########################################################################
  #\eta
  ########################################################################
  cov_eta<-chol2inv(chol(crossprod(x)/sigma2_epsilon + diag(p_x)/(100^2)))
  mean_eta<-cov_eta%*%crossprod(x, (z - b))/sigma2_epsilon
  eta<-rmnorm(n = 1,
              mean = mean_eta,
              varcov = cov_eta)
  
  #############################################################
  #b
  #############################################################
  cov_b<-chol2inv(chol(diag(n)/sigma2_epsilon + CAR_corr/tau2))
  mean_b<-cov_b%*%(z - x%*%eta)/sigma2_epsilon
  b<-rmnorm(n = 1,
            mean = mean_b,
            varcov = cov_b)
  b<-b -
    mean(b)
  
  ####################################################
  #\sigma^2_{\epsilon}
  ####################################################
  rate_new<-crossprod(z - x%*%eta - b)/2.00 +
    0.01
  sigma2_epsilon<-1.00/rgamma(n = 1,
                              shape = (n/2.00 + 0.01),
                              rate = rate_new)
  
  ###########################################
  #\tau^2
  ###########################################
  rate_new<-crossprod(b, CAR_corr%*%b)/2.00 +
    0.01
  tau2<-1.00/rgamma(n = 1,
                    shape = (n/2.00 + 0.01),
                    rate = rate_new)
  
  ##############################################################
  #\rho
  ##############################################################
  rho_old<-rho   
  CAR_corr_old<-CAR_corr
  
  second<-rho_fun(log(rho_old/(1.00 - rho_old)),
                  CAR_corr_old)
  
  rho_trans<-rnorm(n = 1, 
                   mean = log(rho_old/(1.00 - rho_old)), 
                   sd = rho_metrop_sd)
  rho<-exp(rho_trans)/(1.00 + exp(rho_trans))
  CAR_corr<-rho*(diag(rowSums(neighbors_mat)) - neighbors_mat) +
    (1.00 - rho)*diag(n)
  first<-rho_fun(rho_trans,
                 CAR_corr)
  
  ratio<-exp(first - second)   
  acc<-1
  if(ratio < runif(n = 1, min = 0.00, max = 1.00)){
    
    rho<-rho_old
    CAR_corr<-CAR_corr_old
    acc<-0
    
  }
  
  acctot_rho<-acctot_rho + 
    acc
  
  ######################################
  #Information to Keep
  ######################################
  eta_keep[i,]<-eta
  b_keep[i,]<-b
  sigma2_epsilon_keep[i]<-sigma2_epsilon
  tau2_keep[i]<-tau2
  rho_keep[i]<-rho
  
  ########################################
  #Printing to Screen
  ########################################
  if((i %% 1000) == 0){
    
    print(c("Completion %: ",
            round((100*i/mcmc_samples), 2)))
    
    print(c("\rho Acceptance %: ",
            round((100*acctot_rho/i), 2)))
    
  }
  
}






