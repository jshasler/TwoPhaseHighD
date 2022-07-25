##########################################################
##########################################################
## This file includes the functions used in the paper:
## "A Two-Stage Modeling Approach to Building and Evaluating Risk Prediction Models Using 
## High-Dimensional Two-Phase Data"
##########################################################
##########################################################

data_gen <- function(N,prob_y,dim_X,dim_Z,nonzero_X,Z_logOR,rho) {
  # generate full data, first variable is binary
  Cov <- func_A1gen(rho, (dim_X+dim_Z))
  XZ <- mvrnorm(N, rep(0,(dim_X+dim_Z)), Cov)
  X_binary <- Cont2Bin(XZ[,1],prev_X_efficient,1,N)
  X <- cbind(X_binary,XZ[,2:(dim_X)])
  Z <- XZ[,(dim_X+1):(dim_X+dim_Z)]
  Predictor_new <- cbind(X,Z)
  
  # organize vector of coefficients for model used to generate the outcome, Y
  X_logOR <- NULL
  for (i in 1:nonzero_X) {
    X_logOR <- c(X_logOR,nonzeroX_logOR[i],rep(log(1),dim_X/nonzero_X-1))
  }
  gamma_logOR <- c(X_logOR,Z_logOR)
  
  # solve for the intercept term
  find_gamma0 <- function(gamma0){
    exp_val = exp(gamma0 + Predictor_new%*%as.matrix(gamma_logOR))
    f <- mean(exp_val/(1 + exp_val)) - prob_y
  }
  gamma0 = uniroot(f = find_gamma0, c(-20, 20))$root
  exp_val = exp(gamma0 + Predictor_new%*%as.matrix(gamma_logOR))
  #epsilon_beta <- abs(c(gamma0,nonzeroX_logOR,Z_logOR))/numbeta
  
  # generate the outcome, Y
  prob <- exp_val/(1+exp_val)
  Y <- rbinom(N,1,prob)
  
  data_P1 <- data.frame(cbind(X,Z,Y))
  return(list(data_P1=data_P1,gamma_logOR=gamma_logOR,gamma0=gamma0))
}

func_A1gen <- function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  A1
}

Cont2Bin <- function(X,freq,p_bin,n) {
  X_Cont2Bin <- matrix(0,n,p_bin)
  for (i in 1:p_bin) {
    cutoff <- unname(quantile(X,1-freq[i]))
    for (j in 1:n) {
      X_Cont2Bin[j] <- X[j] >= cutoff
    }
    
  }
  X_Cont2Bin
}

Dat_gensrs <- function(prob_nII,n) {
  R <- rbinom(N,1,prob_nII)
  return(R)
}

Dat_fullFit <- function(data,gamma_logOR) {
  tmp <- cbind(data$Y,data[,which(gamma_logOR!=0)])
  names(tmp)[1] <- "Y"
  mod <- glm(Y~.,data=tmp,family=binomial)
  yhat.full <- predict(mod,tmp,type="response")
  beta.full <- mod$coefficients
  return(list(yhat.full=yhat.full, beta.full=beta.full))
}

Dat_adaLasso <- function(HighD_X,Y) {
  ridge <- cv.glmnet(HighD_X,Y,family="binomial",alpha=0)
  best_ridge_coef <- coef(ridge, s = ridge$lambda.min)
  best_ridge_coef <- as.numeric(best_ridge_coef)[-1]
  adaptive_lasso <- cv.glmnet(x=HighD_X,y=Y,family="binomial",alpha=1,penalty.factor=1/abs(best_ridge_coef))
  yhat.full <- predict(adaptive_lasso,s="lambda.min",HighD_X,type="response")
  return(list(fit=adaptive_lasso))
}

Dat_fullauc <- function(yhat) {
  # calculate AUC 
  c <- c(seq(0, 1, by=0.01)*max(yhat),1)
  denom_tpr <- sum(yhat)
  denom_fpr <- sum(1-yhat)
  num_tpr <- numeric(length(c))
  num_fpr <- numeric(length(c))
  for (j in 1:length(c)) {
    ind <- (yhat>=c[j])
    num_tpr[j] <- sum(yhat[ind==1])
    num_fpr[j] <- sum((1-yhat)[ind==1])
  }
  tpr <- num_tpr/denom_tpr
  fpr <- num_fpr/denom_fpr
  auc <- MESS::auc(fpr,tpr)
  tpr_05 <- sum(yhat[yhat>=0.5])/denom_tpr
  fpr_05 <- sum((1-yhat)[yhat>=0.5])/denom_fpr
  return(list(auc=auc,tpr_05=tpr_05,fpr_05=fpr_05))
}

Dat_PLFit <- function(data,ada_lasso_fit,dim_X,dim_Z,N) {
  HighD_X <- as.matrix(data[c(1:dim_X)])
  coef <- coef(ada_lasso_fit,s="lambda.min")
  coef_selected <- rownames(coef)[which(coef!=0)]
  coef_gamma <- coef[which(coef!=0)]
  X_gamma <- cbind(rep(1,N),HighD_X)%*%coef
  prob_initial_Y1 <- predict(ada_lasso_fit,HighD_X,type="response",s="lambda.min")
  prob_initial_Y0 <- 1-prob_initial_Y1
  data <- cbind(data,prob_initial_Y1,prob_initial_Y0)
  names(data)[306] <- "prob_initial_Y1"
  names(data)[307] <- "prob_initial_Y0"
  
  # score function for gamma
  X1 <- as.matrix(cbind(rep(1,nrow(data)),data[1:(dim_X)]))
  X1_A <- X1[,which(coef!=0)]
  dim_X_A <- ncol(X1_A)
  S_gamma <- apply(cbind(data$Y,X1_A, prob_initial_Y1), 1, function(x) x[2:(dim_X_A+1)]*(x[1]-x[dim_X_A+2]))
  
  # Fisher information matrix for gamma
  I_gamma <- solve(t(X1_A) %*% diag(as.vector(prob_initial_Y1*(prob_initial_Y0))) %*% X1_A/N)
  
  # influence function for gamma
  h_gamma <-  I_gamma %*% S_gamma
  
  # create dataset of phase I cases and phase I controls
  data_0 <- data[which(data$Y==0),]
  data_1 <- data[which(data$Y==1),]
  
  # fit missingness model for cases and controls
  f1_fit <- gam(R~bs(prob_initial_Y0,knots = (1-prob_y),Boundary.knots = c(0,1)),data=data_1,family=binomial)
  f0_fit <- gam(R~bs(prob_initial_Y1,knots = prob_y,Boundary.knots = c(0,1)),data=data_0,family=binomial)
  
  # calculate predicted probabilities for offset
  f1_hat <- predict(f1_fit,newdata=data,type="response")
  f0_hat <- predict(f0_fit,newdata=data,type="response")
  
  # calculate X matrix for f1 and f0
  Xmat_f1 <- cbind(rep(1,nrow(data)),bs(data$prob_initial_Y0,knots = (1-prob_y),Boundary.knots = c(0,1)))
  Xmat_f0 <- cbind(rep(1,nrow(data)),bs(data$prob_initial_Y1,knots = prob_y,Boundary.knots = c(0,1)))
  
  # calculate score function for tau1 and tau0
  S_tau1 <- apply(cbind(data$Y,data$R,Xmat_f1, f1_hat), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f1)-1)]*(x[2]-x[(3+ncol(Xmat_f1))]))
  S_tau0 <- apply(cbind((1-data$Y),data$R,Xmat_f0, f0_hat), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f0)-1)]*(x[2]-x[(3+ncol(Xmat_f0))]))
  
  # A matrix for tau (derivative of score with respect to tau)
  A_tau1 <- -solve(t(Xmat_f1) %*% diag(data$Y*as.vector(f1_hat*(1-f1_hat))) %*% (Xmat_f1)/N)
  A_tau0 <- -solve(t(Xmat_f0) %*% diag((1-data$Y)*as.vector(f0_hat*(1-f0_hat))) %*% (Xmat_f0)/N)
  
  # B matrix for tau (derivative of score with respect to gamma) - calculate numerically
  b1gammaforjaco <- function(gammader) {
    p0_jaco <- 1-exp(X1_A%*%gammader)/(1+exp(X1_A%*%gammader))
    
    # create dataset of phase I cases and phase I controls
    data_jaco <- cbind(data,p0_jaco)
    
    #calculate X matrix
    Xmat_f1_jaco <- cbind(rep(1,nrow(data_jaco)),bs(data_jaco$p0_jaco,knots=(1-prob_y),Boundary.knots = c(0,1)))
    f1_hat_jaco <- exp(Xmat_f1_jaco%*%f1_fit$coefficients)/(1+exp(Xmat_f1_jaco%*%f1_fit$coefficients))
    
    
    big_score <- apply(cbind(data_jaco$Y,data_jaco$R,Xmat_f1_jaco, f1_hat_jaco), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f1_jaco)-1)]*(x[2]-x[(3+ncol(Xmat_f1_jaco))]))
    mean_score <- rowMeans(big_score)
    return(mean_score)
  }
  
  b0gammaforjaco <- function(gammader) {
    p1_jaco <- exp(X1_A%*%gammader)/(1+exp(X1_A%*%gammader))
    
    # create dataset of phase I cases and phase I controls
    data_jaco <- cbind(data,p1_jaco)
    
    #calculate X matrix
    Xmat_f0_jaco <- cbind(rep(1,nrow(data_jaco)),bs(data_jaco$p1_jaco,knots=prob_y,Boundary.knots = c(0,1)))
    f0_hat_jaco <- exp(Xmat_f0_jaco%*%f0_fit$coefficients)/(1+exp(Xmat_f0_jaco%*%f0_fit$coefficients))
    
    
    big_score <- apply(cbind((1-data_jaco$Y),data_jaco$R,Xmat_f0_jaco, f0_hat_jaco), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f0_jaco)-1)]*(x[2]-x[(3+ncol(Xmat_f0_jaco))]))
    mean_score <- rowMeans(big_score)
    return(mean_score)
  }
  
  gammader <- (coef)[which(coef!=0)]
  B_tau1_gamma <- jacobian(b1gammaforjaco,gammader)
  B_tau0_gamma <- jacobian(b0gammaforjaco,gammader)
  
  # Calculate influence function for tau
  h_tau1 <- A_tau1 %*% (-B_tau1_gamma %*% h_gamma - S_tau1)
  h_tau0 <- A_tau0 %*% (-B_tau0_gamma %*% h_gamma - S_tau0)
  
  # calculate variance for tau
  coef_tau1 <- f1_fit$coefficients
  coef_tau0 <- f0_fit$coefficients
  
  # estimating equation for our two phase method
  XZ <- as.matrix(cbind(rep(1,nrow(data)),as.numeric(X_gamma),data[,(dim_X+1):(dim_X+dim_Z)]))
  estequal <- function(x){
    U <-  colSums(XZ*data$R*as.vector(data$Y - (exp(XZ%*%x)*as.vector(f1_hat/f0_hat))/(1 + exp(XZ%*%x)*as.vector(f1_hat/f0_hat))))
    return(U)
  }
  
  # solve the score function
  inibeta=rep(0,dim_Z+2)
  coef <- nleqslv(inibeta,estequal)$x
  
  # yscore.pl is equivalent to u_i
  yscore.pl <- exp(XZ%*%coef)*as.vector(f1_hat/f0_hat)/(1+exp(XZ%*%coef)*as.vector(f1_hat/f0_hat))
  
  # Calculate U for estimating equation
  U <- apply(cbind(data$R,data$Y,XZ, yscore.pl), 1, function(x) x[1]*x[3:(length(coef)+2)]*(x[2]-x[length(coef)+3]))
  
  # Calculate D_beta_inv
  neg_D_beta_inv <- solve(t(data$R*XZ) %*% diag(as.vector(yscore.pl*(1-yscore.pl))) %*% (XZ)/N)
  
  dgammaforjaco <- function(gammader) {
    p1_jaco <- exp(X1_A%*%gammader)/(1+exp(X1_A%*%gammader))
    p0_jaco <- 1-p1_jaco
    X_gamma_jaco <- X1_A%*%gammader
    
    # create dataset of phase I cases and phase I controls
    data_jaco <- cbind(data,p1_jaco,p0_jaco)
    
    #calculate X matrix
    Xmat_f0_jaco <- cbind(rep(1,nrow(data_jaco)),bs(data_jaco$p1_jaco,knots=prob_y,Boundary.knots = c(0,1)))
    Xmat_f1_jaco <- cbind(rep(1,nrow(data_jaco)),bs(data_jaco$p0_jaco,knots=(1-prob_y),Boundary.knots = c(0,1)))
    f0_hat_jaco <- exp(Xmat_f0_jaco%*%f0_fit$coefficients)/(1+exp(Xmat_f0_jaco%*%f0_fit$coefficients))
    f1_hat_jaco <- exp(Xmat_f1_jaco%*%f1_fit$coefficients)/(1+exp(Xmat_f1_jaco%*%f1_fit$coefficients))
    XZ_jaco <- as.matrix(cbind(rep(1,nrow(data)),as.numeric(X_gamma_jaco),data[,(dim_X+1):(dim_X+dim_Z)]))
    yhat_jaco <- exp(XZ_jaco%*%coef)*(as.vector(f1_hat_jaco/f0_hat_jaco))/(1+exp(XZ_jaco%*%coef)*as.vector(f1_hat_jaco/f0_hat_jaco))
    
    #calculate score
    big_U1 <- apply(cbind(data_jaco$R,data_jaco$Y,XZ_jaco,yhat_jaco), 1, function(x) x[1]*x[3:(length(coef)+2)]*(x[2]-x[length(coef)+3]))
    mean_score <- rowMeans(big_U1)
    return(mean_score)
  }
  D_gamma <- jacobian(dgammaforjaco,gammader)
  
  # ystar.pl is equivalent to xi_i
  ystar.pl <- exp(XZ%*%coef)/(1+exp(XZ%*%coef))
  
  # Calculate D_tau1 and D_tau0
  D_tau1 <- t(data$R*XZ) %*% diag(as.vector(-yscore.pl*(1-yscore.pl)*(1-f1_hat))) %*% (Xmat_f1)/N
  D_tau0 <- t(data$R*XZ) %*% diag(as.vector(yscore.pl*(1-yscore.pl)*(1-f0_hat))) %*% (Xmat_f0)/N
  
  # Claculate influence function for beta
  h_beta <- neg_D_beta_inv %*% (U + D_tau1 %*% h_tau1 + D_tau0 %*% h_tau0 + D_gamma %*% h_gamma)
  var_coef <- diag(h_beta %*% t(h_beta)/N^2)
  
  # estimate of f(x,z)
  fhat.pl <- (data$R)/(as.vector(f1_hat)*ystar.pl+as.vector(f0_hat)*(1-ystar.pl))/N
  
  
  return(list(beta.pl=coef,yhat.pl=ystar.pl,fhat.pl=fhat.pl,h_beta=h_beta,h_tau0=h_tau0,h_tau1=h_tau1,h_gamma=h_gamma,
              f1_hat=f1_hat,f0_hat=f0_hat,Xmat_f1=Xmat_f1,Xmat_f0=Xmat_f0,var=var_coef,
              coef_tau1=coef_tau1,coef_tau0=coef_tau0,coef_gamma=coef_gamma,XZ=XZ,X1_A=X1_A))
}

Dat_auc_PL <- function(data,plFit_list) {
  yhat.pl <- plFit_list$yhat.pl
  fhat.pl <- plFit_list$fhat.pl
  beta.pl <- plFit_list$beta.pl
  h_beta <- plFit_list$h_beta
  h_tau0 <- plFit_list$h_tau0
  h_tau1 <- plFit_list$h_tau1
  h_gamma <- plFit_list$h_gamma
  f1_hat <- as.vector(plFit_list$f1_hat)
  f0_hat <- as.vector(plFit_list$f0_hat)
  tau1.pl <- plFit_list$coef_tau1
  tau0.pl <- plFit_list$coef_tau0
  Xmat_f1 <- plFit_list$Xmat_f1
  Xmat_f0 <- plFit_list$Xmat_f0
  XZ <- plFit_list$XZ
  gamma <- plFit_list$coef_gamma
  h_gamma <- plFit_list$h_gamma
  X1_A <- plFit_list$X1_A
  
  # calculate the AUC
  denom_tpr.pl <- sum(yhat.pl*fhat.pl)
  denom_fpr.pl <- sum((1-yhat.pl)*fhat.pl)
  c <- c(seq(0, 1, by=0.01)*max(yhat.pl*data$R),1)
  num_tpr.pl <- numeric(length(c))
  num_fpr.pl <- numeric(length(c))
  for (j in 1:length(c)){
    ind <- (yhat.pl>=c[j])
    num_tpr.pl[j] <- sum((yhat.pl*fhat.pl)[ind==1])
    num_fpr.pl[j] <- sum(((1-yhat.pl)*fhat.pl)[ind==1])
  }
  tpr.pl <- num_tpr.pl/denom_tpr.pl
  fpr.pl <- num_fpr.pl/denom_fpr.pl
  auc <- MESS::auc(fpr.pl, tpr.pl)
  
  # numerical method for dFPR/dbeta and dTPR/dbeta as a function of c
  # Estimate df/dbeta
  df_beta <- apply(cbind(XZ, -fhat.pl^2*N*(f1_hat-f0_hat)*yhat.pl*(1-yhat.pl)), 1, function(x) x[1:(dim_Z+2)]*x[(dim_Z+2)+1])
  
  dFPR_beta <- matrix(0,nrow=length(c),ncol=(dim_Z+2))
  dTPR_beta <- matrix(0,nrow=length(c),ncol=(dim_Z+2))
  yfit1 <- numeric(nrow(data))
  yfit2 <- numeric(nrow(data))
  for (i in 1:(dim_Z+2)) {
    tmp1 <- beta.pl
    tmp2 <- beta.pl
    tmp1[i] <- beta.pl[i] + 0.001
    tmp2[i] <- beta.pl[i] - 0.001
    yfit1 <- exp(XZ %*% tmp1)/(1+exp(XZ %*% tmp1))
    yfit2 <- exp(XZ %*% tmp2)/(1+exp(XZ %*% tmp2))
    fhat1 <- fhat.pl + df_beta[i,]*0.001
    fhat2 <- fhat.pl - df_beta[i,]*0.001
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yfit1*fhat1)
    denom2_tpr <- sum(yfit2*fhat2)
    denom1_fpr <- sum((1-yfit1)*fhat1)
    denom2_fpr <- sum((1-yfit2)*fhat2)
    
    for (j in 1:length(c)){
      ind1 <- (yfit1>=c[j])
      ind2 <- (yfit2>=c[j])
      num1_tpr <- sum((yfit1*fhat1)[ind1==1])
      num2_tpr <- sum((yfit2*fhat2)[ind2==1])
      num1_fpr <- sum(((1-yfit1)*fhat1)[ind1==1])
      num2_fpr <- sum(((1-yfit2)*fhat2)[ind2==1])
      dTPR_beta[j,i] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*0.001)
      dFPR_beta[j,i] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*0.001)
    }
  }
  numtau1 <- length(tau1.pl)
  df_tau1 <- apply(cbind(Xmat_f1, -fhat.pl^2*N*f1_hat*(1-f1_hat)*yhat.pl), 1, function(x) x[1:numtau1]*x[numtau1+1])
  dFPR_tau1 <- matrix(0, nrow=length(c), ncol=numtau1)
  dTPR_tau1 <- matrix(0, nrow=length(c), ncol=numtau1)
  
  for (i in 1:numtau1) {
    fhat1 <- fhat.pl + df_tau1[i,]*0.001
    fhat2 <- fhat.pl - df_tau1[i,]*0.001
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yhat.pl*fhat1)
    denom2_tpr <- sum(yhat.pl*fhat2)
    denom1_fpr <- sum((1-yhat.pl)*fhat1)
    denom2_fpr <- sum((1-yhat.pl)*fhat2)
    
    for (j in 1:length(c)){
      ind <- (yhat.pl>=c[j])
      num1_tpr <- sum((yhat.pl*fhat1)[ind==1])
      num2_tpr <- sum((yhat.pl*fhat2)[ind==1])
      num1_fpr <- sum(((1-yhat.pl)*fhat1)[ind==1])
      num2_fpr <- sum(((1-yhat.pl)*fhat2)[ind==1])
      
      dTPR_tau1[j,i] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*0.001)
      dFPR_tau1[j,i] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*0.001)
    }
  }
  
  numtau0 <- length(tau0.pl)
  df_tau0 <- apply(cbind(Xmat_f0, -fhat.pl^2*N*f0_hat*(1-f0_hat)*(1-yhat.pl)), 1, function(x) x[1:numtau0]*x[numtau0+1])
  dFPR_tau0 <- matrix(0, nrow=length(c), ncol=numtau0)
  dTPR_tau0 <- matrix(0, nrow=length(c), ncol=numtau0)
  for (i in 1:numtau1) {
    tmp1 <- tau0.pl
    tmp2 <- tau0.pl
    fhat1 <- fhat.pl + df_tau0[i,]*0.001
    fhat2 <- fhat.pl - df_tau0[i,]*0.001
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yhat.pl*fhat1)
    denom2_tpr <- sum(yhat.pl*fhat2)
    denom1_fpr <- sum((1-yhat.pl)*fhat1)
    denom2_fpr <- sum((1-yhat.pl)*fhat2)
    
    for (j in 1:length(c)){
      ind <- (yhat.pl>=c[j])
      num1_tpr <- sum((yhat.pl*fhat1)[ind==1])
      num2_tpr <- sum((yhat.pl*fhat2)[ind==1])
      num1_fpr <- sum(((1-yhat.pl)*fhat1)[ind==1])
      num2_fpr <- sum(((1-yhat.pl)*fhat2)[ind==1])
      
      dTPR_tau0[j,i] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*0.001)
      dFPR_tau0[j,i] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*0.001)
    }
  }
  
  # calculate the derivative of f with respect to gamma numerically
  df_dgamma <- matrix(0,nrow=nrow(data),ncol=length(gamma))
  for (i in 1:length(gamma)) {
    tmp1 <- gamma
    tmp2 <- gamma
    tmp1[i] <- gamma[i] + 0.001
    tmp2[i] <- gamma[i] - 0.001
    X_gamma1 <- X1_A%*%tmp1
    X_gamma2 <- X1_A%*%tmp2
    X_mat1 <- as.matrix(XZ)
    X_mat1[,2] <- X_gamma1
    X_mat2 <- as.matrix(XZ)
    X_mat2[,2] <- X_gamma2
    yhat1 <- exp(X_mat1 %*% beta.pl)/(1+exp(X_mat1 %*% beta.pl))
    yhat2 <- exp(X_mat2 %*% beta.pl)/(1+exp(X_mat2 %*% beta.pl))
    p1_1 <- exp(X1_A%*%tmp1)/(1+exp(X1_A%*%tmp1))
    p1_2 <- exp(X1_A%*%tmp2)/(1+exp(X1_A%*%tmp2))
    p0_1 <- 1-p1_1
    p0_2 <- 1-p1_2
    Xmat_f0_1 <- cbind(rep(1,nrow(data)),bs(p1_1,knots=prob_y,Boundary.knots = c(0,1)))
    Xmat_f0_2 <- cbind(rep(1,nrow(data)),bs(p1_2,knots=prob_y,Boundary.knots = c(0,1)))
    Xmat_f1_1 <- cbind(rep(1,nrow(data)),bs(p0_1,knots=(1-prob_y),Boundary.knots = c(0,1)))
    Xmat_f1_2 <- cbind(rep(1,nrow(data)),bs(p0_2,knots=(1-prob_y),Boundary.knots = c(0,1)))
    f0_hat_1 <- exp(Xmat_f0_1%*%tau0.pl)/(1+exp(Xmat_f0_1%*%tau0.pl))
    f0_hat_2 <- exp(Xmat_f0_2%*%tau0.pl)/(1+exp(Xmat_f0_2%*%tau0.pl))
    f1_hat_1 <- exp(Xmat_f1_1%*%tau1.pl)/(1+exp(Xmat_f1_1%*%tau1.pl))
    f1_hat_2 <- exp(Xmat_f1_2%*%tau1.pl)/(1+exp(Xmat_f1_2%*%tau1.pl))
    fhat1 <- data$R/(f1_hat_1*yhat1+f0_hat_1*(1-yhat1))/N
    fhat2 <- data$R/(f1_hat_2*yhat2+f0_hat_2*(1-yhat2))/N
    df_dgamma[,i] <- (fhat1-fhat2)/(2*0.001)
  }
  
  # numerical method for dFPR/dgamma and dTPR/dgamma as a function of c
  dFPR_gamma <- matrix(0,nrow=length(c),ncol=length(gamma))
  dTPR_gamma <- matrix(0,nrow=length(c),ncol=length(gamma))
  for (i in 1:length(gamma)) {
    tmp1 <- gamma
    tmp2 <- gamma
    tmp1[i] <- gamma[i] + 0.001
    tmp2[i] <- gamma[i] - 0.001
    X_gamma1 <- X1_A%*%tmp1
    X_gamma2 <- X1_A%*%tmp2
    X_mat1 <- XZ
    X_mat1[,2] <- X_gamma1
    X_mat2 <- XZ
    X_mat2[,2] <- X_gamma2
    yfit1 <- exp(X_mat1 %*% beta.pl)/(1+exp(X_mat1 %*% beta.pl))
    yfit2 <- exp(X_mat2 %*% beta.pl)/(1+exp(X_mat2 %*% beta.pl))
    fhat1 <- fhat.pl + df_dgamma[,i]*0.001
    fhat2 <- fhat.pl - df_dgamma[,i]*0.001
    #fhat1 <- data$R/(f1_hat*yfit1+f0_hat*(1-yfit1))
    #fhat2 <- data$R/(f1_hat*yfit2+f0_hat*(1-yfit2))
    #fhat1 <- fhat.pl
    #fhat2 <- fhat.pl
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yfit1*fhat1)
    denom2_tpr <- sum(yfit2*fhat2)
    denom1_fpr <- sum((1-yfit1)*fhat1)
    denom2_fpr <- sum((1-yfit2)*fhat2)
    
    for (j in 1:length(c)){
      ind1 <- (yfit1>=c[j])
      ind2 <- (yfit2>=c[j])
      num1_tpr <- sum((yfit1*fhat1)[ind1==1])
      num2_tpr <- sum((yfit2*fhat2)[ind2==1])
      num1_fpr <- sum(((1-yfit1)*fhat1)[ind1==1])
      num2_fpr <- sum(((1-yfit2)*fhat2)[ind2==1])
      dTPR_gamma[j,i] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*0.001)
      dFPR_gamma[j,i] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*0.001)
    }
  }
  
  
  # calculate influence function for TPR/FPR as a function of c
  mat_tpr_pl <- matrix(0, nrow=nrow(data), ncol=length(c))
  mat_fpr_pl <- matrix(0, nrow=nrow(data), ncol=length(c))
  for (j in 1:length(c)){
    ind <- (yhat.pl>=c[j])
    mat_tpr_pl[,j] <- N*fhat.pl*yhat.pl*(ind-tpr.pl[j])/denom_tpr.pl # some values greater than 1 here
    mat_fpr_pl[,j] <- N*fhat.pl*(1-yhat.pl)*(ind-fpr.pl[j])/denom_fpr.pl
  }
  
  # calculate the influence function for tpr/fpr
  h_tpr <- t(dTPR_beta %*% h_beta + dTPR_tau1 %*% h_tau1 + dTPR_tau0 %*% h_tau0 + dTPR_gamma %*% h_gamma) + mat_tpr_pl
  h_fpr <- t(dFPR_beta %*% h_beta + dFPR_tau1 %*% h_tau1 + dFPR_tau0 %*% h_tau0 + dFPR_gamma %*% h_gamma) + mat_fpr_pl
  
  
  # calculate influence function for auc
  h_auc <- tryCatch(0-apply(h_fpr, 1, function(x) MESS::auc(tpr.pl, x))+apply(h_tpr, 1, function(x) MESS::auc(fpr.pl, x)),
                    error=function(e){message <- "Cannot solve";message})
  var_auc_pl <- tryCatch(t(h_auc) %*% h_auc/N^2,error=function(e){message <- "Could not solve";message})
  return(list(auc=auc,var=var_auc_pl))
}

Dat_twostage_XZ <- function(data) {
  # generate initial estimates of coefficients for model between Y and X using phase I samples
  HighD_X <- as.matrix(data[c(1:dim_X)])
  ridge <- cv.glmnet(HighD_X,data$Y,family="binomial",alpha=0)
  best_ridge_coef <- coef(ridge, s = ridge$lambda.min)
  best_ridge_coef <- as.numeric(best_ridge_coef)[-1]
  fit <- cv.glmnet(x=HighD_X,y=data$Y,family="binomial",alpha=1,penalty.factor=1/abs(best_ridge_coef))
  coef <- coef(fit,s="lambda.min")
  coef_selected <- rownames(coef)[which(coef!=0)]
  X_gamma <- cbind(rep(1,N),HighD_X)%*%coef
  
  
  # fit the stage II model
  Y <- data$Y
  X_gamma <- as.numeric(X_gamma)
  Z <- data[,(dim_X+1):(dim_X+dim_Z)]
  tmp <- cbind(Y,X_gamma,Z)
  mod_stageII <- glm(Y~.,family="binomial",data=tmp)
  coef_stageII <- mod_stageII$coefficients
  X_mat <- as.matrix(cbind(rep(1,nrow(data)),X_gamma,Z))
  
  # calculate predicted probabilities ystar
  ystar <- exp(X_mat%*%coef_stageII)/(1+exp(X_mat%*%coef_stageII))
  
  return(list(beta=coef_stageII,yhat=ystar))
}