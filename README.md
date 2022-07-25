# TwoPhaseHighD

R code for the paper "A Two-Stage Modeling Approach to Building and Evaluating Risk Prediciton Models Using High-Dimensional Two-Phase Data"

This repository contains the functions and a simulated data example following the simulation scenarios in Section 4 of the paper. All functions can be found in the TwoPhaseHighD_Functions.R file and the simulated data example can be found in the TwoPhaseHighD_Example.R file. Each function in TwoPhaseHighD_Functions.R is described below.

data_gen: This is the main funciton used to generate the data set used in the simulation studies. The data set consists of 300 phase I predictors, $\mathbf{X}$, where the first phase I predictor $X_1$ is a binary variable, 3 phase II predictors $\mathbf{Z}$ and the case/control status $Y$. The characteristics of the data set can be modified by changing the arguments below.
* Input Arguments
  * N: The sample size of the phase I data
  * prob_y: The outcome prevalence
  * dim_X: The number of phase I predictors
  * dim_Z: The number of phase II predictors
  * nonzero_X: The log-odds ratios for the nonzero phase I predictors
  * Z_logOR: The log-odds ratios for the phase II predictors
  * rho: The correlation coefficient
 
func_A1gen: This function is used to create the AR1 covariance matrix for the simulad data set. It is a helper function called within data_gen.

Cont2Bin: This function is used to convert $X_1$ from a continuous to a binary variable based on a pre-specified quantile (freq). It is a helper function called within data_gen.

Dat_gensrs: This function generates the indicator variable $R$ using simple random sampling given a prevalence value (prob_nII) and sample size (n). $R_i=1$ indicats that $\mathbf{Z}$ is available for inidividual $i$ and $R_i=0$ indicates that $\mathbf{Z}$ is not available for individual $i$.

Dat_fullFit: This function fits the data generating logistic regression model which is used as the benchmark in Section 4.2. It takes the full data set and vector of data genenrating log-odds ratios as the input and outputs the vector of log-odds ratios for the model and predicted probabilities.

Dat_fullauc: This function calculates the model based AUC given a vector of predicted probabilities as the input. It can be used to calculate the AUC for the data generating model and the AUC for the two-stage model used as a benchmark in Section 4.1.

Dat_adaLasso: This function fits the stage I adaptive lasso model for the proposed two-stage method. It takes a matrix consisting of the phase I $\mathbf{X}$ variables and the case/control status $Y$ as the input and outputs the adaptive lasso model fit.

Dat_PLFit: This function fits the two-stage model using the proposed pseudo-score method and calculates the asymptotic variance for the coefficients $\mathbf{\beta}$. The output is a list containing the model coefficients (beta.pl), the asymptotic variance of the model coefficients (var), and several additional parameters needed to calculate the auc of the model. The inputs are described below.
* Input Arguments
  * data: The full data set. The data set should contain the phase I variables, followed by the phase II variables, the outcome, and finally the indicator variable $R$
  * ada_lasso_fit: The model fit object from cv.glmnet for the stage I adaptive lasso model. In the simulation, this is the fit object returned by Dat_adaLasso
  * dim_X: The number of phase I variables $\mathbf{X}$ 
  * dim_Z: The number of phase II variables $\mathbf{Z}$
  * N: The sample size of the phase I data

Dat_auc_PL: This function calculates the AUC for the proposed two-stage method and calculates the asymptotic variance for the AUC. It takes a list of inputs from the Dat_PLFit function (plFit_list) and the full data set. The output is a list containing the auc and the asymptotic variance of the auc.

Dat_twostage_XZ: This function fits the two-stage benchmark method from Section 4.1. It takes the data set as the input and returns the coefficient estiamtes and predicted probabilities as the outcome.
