library(MASS)
library(gam)
library(nleqslv)
library(MESS)
library(glmnet)
library(numDeriv)
source("TwoPhaseHighD_Functions.R")
set.seed(333)

########################################################
##### Set Parameters to Generate Example Data Set ######
########################################################

# Total number of individuals in Phase I
N <- 10000

# Prevalence of y
prob_y <- 0.2

# Number of Phase I predictors, X
dim_X <- 300

# Number of non-zero predictors, X
nonzero_X <- 10

# Prevalence of binary X1
prev_X_efficient <- 0.2

# Dimension of Phase II predictors, Z
dim_Z <- 3

# probabilty of being sampled into phase II
prob_n_II <- 0.05

# rho in covariance matrix
rho=0.3

# set log odds ratios for the nonzero X predictors
nonzeroX_logOR <- c(log(2),-log(1.8),log(1.6),-log(1.4),log(1.2),log(2),-log(1.8),log(1.6),-log(1.4),log(1.2))

# set log odds ratios for phase II predictors, Z
Z_logOR <- c(log(1.5),log(2),-log(1.8))

# total number of truly associated variables
numbeta <- length(nonzeroX_logOR)+length(Z_logOR)+1

########################################################
#############  Generate Example Data Set ###############
########################################################

datafull_list <- data_gen(N,prob_y,dim_X,dim_Z,nonzero_X,Z_logOR,rho)
R <- Dat_gensrs(prob_n_II,N)
data.full <- datafull_list$data_P1
data.full <- cbind(data.full,R)
gamma_logOR <- datafull_list$gamma_logOR

##########################################################################
#############  Fit Data Generating Model and Calculate AUC ###############
##########################################################################

fullFit_list <- Dat_fullFit(data.full,gamma_logOR)
full_coef <- fullFit_list$beta.full
full_coef
fullauc_list <- Dat_fullauc(fullFit_list$yhat.full)
auc_full <- fullauc_list$auc
auc_full

##########################################################################
##############  Fit Pseudo-Score Method and Calculate AUC ################
##########################################################################

# fit the stage I adaptive lasso model
P1OnlyFit_List <- Dat_adaLasso(as.matrix(data.full[,c(1:dim_X)]),data.full$Y)
# fit the second stage model and calculate the variance of the coefficients
plFit_list <- Dat_PLFit(data.full,P1OnlyFit_List$fit,dim_X,dim_Z,N)
PL_coef <- plFit_list$beta.pl
PL_coef
var_PL_coef <- plFit_list$var
var_PL_coef
plauc_list <- Dat_auc_PL(data.full,plFit_list)
PL_auc <- plauc_list$auc
PL_auc
PL_auc_var <- plauc_list$var
PL_auc_var

##########################################################################
###########  Fit two-Stage Benchmark Method and Calculate AUC ############
##########################################################################
two_stageXZ_fit <- Dat_twostage_XZ(data.full)
two_stageXZ_coef <- two_stageXZ_fit$beta
two_stageXZ_coef
two_stage_auc_res <- Dat_fullauc(two_stageXZ_fit$yhat)
two_stage_auc <- two_stage_auc_res$auc
two_stage_auc
