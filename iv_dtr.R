library(bridgedist)
library(grf)
library(rgenoud)
library(modi)

expit <- function(x) exp(x)/(1 + exp(x))
scalar1 <- function(x) {x / sqrt(sum(x^2))}

Gen_Data <- function(N){
  L1 <- matrix(rnorm(5 * N, mean = 1.5, sd = 1), nrow = N, ncol = 5)
  Z1 <- rbinom(N, 1, 0.5)
  U1 <- rnorm(N)
  
  delta1 <- pnorm(L1 %*% c(0.8, 0, 0, 0, 0))
  prob_A1 <- pnorm(L1 %*% c(0.8, 0, 0, 0, 0) - 0.2*U1) * (1-delta1) + delta1 * Z1
  A1 <- rbinom(N, 1, prob_A1)
  q1 <- 0.2 + L1 %*% c(-0.6, -0.8, 0, 0, 0)
  e1 <- rnorm(N)
  
  R1 <- 0.5 + L1 %*% c(0.5, 0.8, 0.3, -0.5, 0.7) + A1 * q1 + 1.5*U1 + 1.4*e1
  
  L2 <- matrix(rnorm(5 * N, mean = 2, sd = 1), nrow = N, ncol = 5)
  Z2 <- rbinom(N, 1, 0.5)
  U2 <- rnorm(N)
  delta2 <- pnorm(L2 %*% c(0, 0.7, 0, 0, 0))
  prob_A2 <- pnorm(L2 %*% c(0, 0.7, 0, 0, 0) + 0.1*R1 + 0.5*U2) * (1-delta2) + delta2 * Z2
  A2 <- rbinom(N, 1, prob_A2)
  
  q2 <- L2 %*% c(1, 1, -1, 0, 0) - 0.2 * R1 + 0.4
  
  e2 <- rnorm(N)
  
  R2 <- 3 + L1 %*% c(0.5, 0.5, 0,0,0) + L2 %*% c(0, 0, 0,-1,1) + A2 * q2 + 2*U2 + 1.2*e2
  data <- data.frame(L1 = I(L1), Z1 = Z1, A1 = A1, R1 = R1, L2 = I(L2), Z2 = Z2, A2 = A2, R2 = R2)
  return(data)
}



log_likelihood <- function(params, X, Z, A) {
  n_X <- ncol(X) + 1
  alpha <- params[1:n_X]
  beta <- params[-c(1:n_X)]
  
  # 计算 Phi 的输入
  phi_beta <- as.vector(beta %*% t(cbind(1, X)))
  phi_alpha <- as.vector(alpha %*% t(cbind(1, X)))
  
  # 累积分布函数，数值稳定性调整
  Phi_beta <- pmin(pmax(pnorm(phi_beta), 1e-10), 1 - 1e-10)
  Phi_alpha <- pmin(pmax(pnorm(phi_alpha), 1e-10), 1 - 1e-10)
  
  # 对数似然项，添加保护性检查
  term1 <- log(pmax(Phi_beta * (1 - Phi_alpha) + Z * Phi_alpha, 1e-10)) * A
  term2 <- log(pmax(1 - Phi_beta * (1 - Phi_alpha) - Z * Phi_alpha, 1e-10)) * (1 - A)
  
  # 返回负对数似然
  return(-sum(term1 + term2))
}

learn_policy_mean <- function(data, model, regimeClass.stg1, regimeClass.stg2, print.level, pop.size, wait.generations){
  L1 <- data$L1
  Z1 <- data$Z1
  A1 <- data$A1
  R1 <- data$R1
  L2 <- data$L2
  Z2 <- data$Z2
  A2 <- data$A2
  R2 <- data$R2
  R <- R1 + R2
  s.tol <- diff(range(R))*1e-05
  
  if (model == "LR"){
    glmZ1 <- glm(Z1 ~ L1, family = binomial(link = "logit"))
    pi.Z1 <- Z1 * glmZ1$fitted.values + (1 - Z1) * (1 - glmZ1$fitted.values)
    
    glmZ2 <- glm(Z2 ~ L1+L2, family = binomial(link = "logit"))
    pi.Z2 <- Z2 * glmZ2$fitted.values + (1 - Z2) * (1 - glmZ2$fitted.values)
    
    
    a0 <- rep(0, 1 + ncol(L1))
    b0 <- rep(0, 1 + ncol(L1))
    init_params <- c(a0, b0)  

    result1 <- optim(
      par = init_params,          
      X=L1, Z=Z1, A=A1,
      fn = log_likelihood,        
      method = "BFGS",            
      control = list(maxit = 1000) 
    )
    alpha1 <- result1$par[1:(1 + ncol(L1))]
    delta.L1 <- pnorm(alpha1 %*% t(cbind(1, L1)))
    
    a0 <- rep(0, 1 + ncol(L2))
    b0 <- rep(0, 1 + ncol(L2))
    init_params <- c(a0, b0)  

    result2 <- optim(
      par = init_params,          
      X=L2, Z=Z2, A=A2,
      fn = log_likelihood,        
      method = "BFGS",           
      control = list(maxit = 1000)
    )
    alpha2 <- result2$par[1:(1 + ncol(L2))]
    delta.L2 <- pnorm(alpha2 %*% t(cbind(1, L2)))
  }else if (model == "RF") {
    grfZ1 <- probability_forest(X = cbind(L1), Y = as.factor(Z1))
    pi.Z1 <- Z1 * grfZ1$predictions[, 2] + (1 - Z1) * grfZ1$predictions[, 1]
    
    grfA1 <- probability_forest(X = cbind(L1,Z1), Y = as.factor(A1))
    delta.L1 <- predict(grfA1, cbind(L1, 1))$predictions[, "1"] - predict(grfA1, cbind(L1, 0))$predictions[, "1"]
    
    grfZ2 <- probability_forest(X = cbind(L1,L2), Y = as.factor(Z2))
    pi.Z2 <- Z2 * grfZ2$predictions[, 2] + (1 - Z2) * grfZ2$predictions[, 1]
    
    grfA2 <- probability_forest(X = cbind(L1,L2,Z2), Y = as.factor(A2))
    delta.L2 <- predict(grfA2, cbind(L1,L2, 1))$predictions[, "1"] - predict(grfA2, cbind(L1,L2, 0))$predictions[, "1"]
  } else {
    stop("model must be set to LR or RF")
  }
  
  regimeClass.stg1="A1~L1"
  regimeClass.stg2="A2~R1+L2"
  
  regimeClass.stg1 <- as.formula(regimeClass.stg1)
  regimeClass.stg2 <- as.formula(regimeClass.stg2)
  
  p.data1 <- model.matrix(regimeClass.stg1, data)
  p.data2 <- model.matrix(regimeClass.stg2, data)
  txname.stg1 <- as.character(regimeClass.stg1[[2]])
  txname.stg2 <- as.character(regimeClass.stg2[[2]])
  
  minVec1 <- apply(p.data1, MARGIN = 2, min)
  spanVec1 <- apply(p.data1, MARGIN = 2, FUN=function(x) max(x)-min(x))
  
  minVec2 <- apply(p.data2, MARGIN = 2, min)
  spanVec2 <- apply(p.data2, MARGIN = 2, FUN=function(x) max(x)-min(x))
  
  # Rescale each nonconstant variable in regimeClass to range between 0 and 1
  p.data.scale1 <- cbind(Intercept=1, apply(p.data1, MARGIN = 2, 
                                            FUN = function(x) (x-min(x))/(max(x)-min(x)))[,-1])
  p.data.scale2 <- cbind(Intercept=1, apply(p.data2, MARGIN = 2, 
                                            FUN = function(x) (x-min(x))/(max(x)-min(x)))[,-1])
  
  
  txVec1 <- try(data[, txname.stg1], silent = TRUE)
  txVec2 <- try(data[, txname.stg2], silent = TRUE)
  nvars.stg1 <- ncol(p.data1)
  nvars.stg2 <- ncol(p.data2)
  nvars.total <- nvars.stg1 + nvars.stg2
  
  # define the form of decision rules
  CovSpace_1 <- p.data.scale1
  CovSpace_2 <- p.data.scale2  
  
  
  
  V_mean_hat <- function(d){
    d1 <- d[1:nvars.stg1]
    A1d <- as.numeric(CovSpace_1%*%d1 > 0)
    d2 <- d[-c(1:nvars.stg1)]
    A2d <- as.numeric(CovSpace_2%*%d2 > 0)
    pic <- pi.Z1 * pi.Z2 * delta.L1 * delta.L2
    c <- as.numeric(A1d == A1) * as.numeric(A2d == A2)*(2 * A1 - 1) * (2 * Z1 - 1) *(2 * A2 - 1) * (2 * Z2 - 1) 
    wts <- c/pic
    result <- weighted.mean(wts * R)
    
    return(result)
  }
  
  est <- genoud(fn=V_mean_hat, nvars=nvars.total,
                print.level=print.level, max=TRUE,
                pop.size=pop.size,
                wait.generations=wait.generations,
                gradient.check=FALSE, BFGS=FALSE,
                P1=50, P2=50, P3=50, P4=50, P5=50,
                P6=50, P7=50, P8=50, P9=0,
                default.domains = 1,
                starting.values=NULL,
                hard.generation.limit=FALSE,
                solution.tolerance=s.tol,
                optim.method="Nelder-Mead")
  coef.1 <- scalar1(est$par[1:nvars.stg1])
  coef.2 <- scalar1(est$par[-c(1:nvars.stg1)])
  
  # parameter indexing the same estimated optimal treatment regime, where all
  # covariates are in the original scale
  coef.orgn.scale.1 <- rep(0,length(coef.1))
  coef.orgn.scale.2 <- rep(0,length(coef.2))
  
  coef.orgn.scale.1[1] <-coef.1[1]- sum(coef.1[-1]*minVec1[-1]/spanVec1[-1])
  coef.orgn.scale.2[1] <-coef.2[1]- sum(coef.2[-1]*minVec2[-1]/spanVec2[-1])
  
  coef.orgn.scale.1[-1] <- coef.1[-1]/spanVec1[-1]
  coef.orgn.scale.2[-1] <- coef.2[-1]/spanVec2[-1]
  
  coef.orgn.scale.1 <- scalar1(coef.orgn.scale.1)
  coef.orgn.scale.2 <- scalar1(coef.orgn.scale.2)
  
  names(coef.1) <- names(coef.orgn.scale.1)<- colnames(p.data.scale1) 
  names(coef.2) <- names(coef.orgn.scale.2)<- colnames(p.data.scale2)
  
  out <- list(coef.1 = coef.1, 
              coef.orgn.scale.1 =coef.orgn.scale.1,
              coef.2=coef.2,
              coef.orgn.scale.2 =coef.orgn.scale.2)
  return(out)
}

learn_policy_quan <- function(data, model, tau, regimeClass.stg1, regimeClass.stg2, print.level, pop.size, wait.generations){
  L1 <- data$L1
  Z1 <- data$Z1
  A1 <- data$A1
  R1 <- data$R1
  L2 <- data$L2
  Z2 <- data$Z2
  A2 <- data$A2
  R2 <- data$R2
  R <- R1 + R2
  s.tol <- diff(range(R))*1e-05
  
  if (model == "LR"){
    glmZ1 <- glm(Z1 ~ L1, family = binomial(link = "logit"))
    pi.Z1 <- Z1 * glmZ1$fitted.values + (1 - Z1) * (1 - glmZ1$fitted.values)
    
    glmZ2 <- glm(Z2 ~ L1+L2, family = binomial(link = "logit"))
    pi.Z2 <- Z2 * glmZ2$fitted.values + (1 - Z2) * (1 - glmZ2$fitted.values)
    
    
    a0 <- rep(0, 1 + ncol(L1))
    b0 <- rep(0, 1 + ncol(L1))
    # 设置初始值 (alpha_k, beta_k)
    init_params <- c(a0, b0) 
    
    # 优化
    result1 <- optim(
      par = init_params,         
      X=L1, Z=Z1, A=A1,
      fn = log_likelihood,        
      method = "BFGS",           
      control = list(maxit = 1000)
    )
    alpha1 <- result1$par[1:(1 + ncol(L1))]
    delta.L1 <- pnorm(alpha1 %*% t(cbind(1, L1)))
    
    a0 <- rep(0, 1 + ncol(L2))
    b0 <- rep(0, 1 + ncol(L2))
    init_params <- c(a0, b0) 

    result2 <- optim(
      par = init_params,         
      X=L2, Z=Z2, A=A2,
      fn = log_likelihood,      
      method = "BFGS",            
      control = list(maxit = 1000) 
    )
    alpha2 <- result2$par[1:(1 + ncol(L2))]
    delta.L2 <- pnorm(alpha2 %*% t(cbind(1, L2)))
  }else if (model == "RF") {
    grfZ1 <- probability_forest(X = cbind(L1), Y = as.factor(Z1))
    pi.Z1 <- Z1 * grfZ1$predictions[, 2] + (1 - Z1) * grfZ1$predictions[, 1]
    
    grfA1 <- probability_forest(X = cbind(L1,Z1), Y = as.factor(A1))
    delta.L1 <- predict(grfA1, cbind(L1, 1))$predictions[, "1"] - predict(grfA1, cbind(L1, 0))$predictions[, "1"]
    
    grfZ2 <- probability_forest(X = cbind(L2), Y = as.factor(Z2))
    pi.Z2 <- Z2 * grfZ2$predictions[, 2] + (1 - Z2) * grfZ2$predictions[, 1]
    
    grfA2 <- probability_forest(X = cbind(R1, L2,Z2), Y = as.factor(A2))
    delta.L2 <- predict(grfA2, cbind(R1, L2, 1))$predictions[, "1"] - predict(grfA2, cbind(R1, L2, 0))$predictions[, "1"]
  } else {
    stop("model must be set to LR or RF")
  }
  
  regimeClass.stg1="A1~L1"
  regimeClass.stg2="A2~R1+L2"
  
  regimeClass.stg1 <- as.formula(regimeClass.stg1)
  regimeClass.stg2 <- as.formula(regimeClass.stg2)
  
  p.data1 <- model.matrix(regimeClass.stg1, data)
  p.data2 <- model.matrix(regimeClass.stg2, data)
  txname.stg1 <- as.character(regimeClass.stg1[[2]])
  txname.stg2 <- as.character(regimeClass.stg2[[2]])
  
  minVec1 <- apply(p.data1, MARGIN = 2, min)
  spanVec1 <- apply(p.data1, MARGIN = 2, FUN=function(x) max(x)-min(x))
  
  minVec2 <- apply(p.data2, MARGIN = 2, min)
  spanVec2 <- apply(p.data2, MARGIN = 2, FUN=function(x) max(x)-min(x))
  
  # Rescale each nonconstant variable in regimeClass to range between 0 and 1
  p.data.scale1 <- cbind(Intercept=1, apply(p.data1, MARGIN = 2, 
                                            FUN = function(x) (x-min(x))/(max(x)-min(x)))[,-1])
  p.data.scale2 <- cbind(Intercept=1, apply(p.data2, MARGIN = 2, 
                                            FUN = function(x) (x-min(x))/(max(x)-min(x)))[,-1])
  
  
  txVec1 <- try(data[, txname.stg1], silent = TRUE)
  txVec2 <- try(data[, txname.stg2], silent = TRUE)
  nvars.stg1 <- ncol(p.data1)
  nvars.stg2 <- ncol(p.data2)
  nvars.total <- nvars.stg1 + nvars.stg2
  
  # define the form of decision rules
  CovSpace_1 <- p.data.scale1
  CovSpace_2 <- p.data.scale2  
  
  
  
  V_quan_hat <- function(d, tau){
    d1 <- d[1:nvars.stg1]
    A1d <- as.numeric(CovSpace_1%*%d1 > 0)
    d2 <- d[-c(1:nvars.stg1)]
    A2d <- as.numeric(CovSpace_2%*%d2 > 0)
    pic <- pi.Z1 * pi.Z2 * delta.L1 * delta.L2
    c <- as.numeric(A1d == A1) * as.numeric(A2d == A2)*(2 * A1 - 1) * (2 * Z1 - 1) *(2 * A2 - 1) * (2 * Z2 - 1) 
    wts <- c/pic
    result <- weighted.quantile(R, wts,  prob = tau, plot = FALSE)
    
    return(result)
  }
  
  est <- genoud(fn=V_quan_hat, nvars=nvars.total,tau =tau,
                print.level=print.level, max=TRUE,
                pop.size=pop.size,
                wait.generations=wait.generations,
                gradient.check=FALSE, BFGS=FALSE,
                P1=50, P2=50, P3=50, P4=50, P5=50,
                P6=50, P7=50, P8=50, P9=0,
                default.domains = 1,
                starting.values=NULL,
                hard.generation.limit=FALSE,
                solution.tolerance=s.tol,
                optim.method="Nelder-Mead")
  coef.1 <- scalar1(est$par[1:nvars.stg1])
  coef.2 <- scalar1(est$par[-c(1:nvars.stg1)])
  
  # parameter indexing the same estimated optimal treatment regime, where all
  # covariates are in the original scale
  coef.orgn.scale.1 <- rep(0,length(coef.1))
  coef.orgn.scale.2 <- rep(0,length(coef.2))
  
  coef.orgn.scale.1[1] <-coef.1[1]- sum(coef.1[-1]*minVec1[-1]/spanVec1[-1])
  coef.orgn.scale.2[1] <-coef.2[1]- sum(coef.2[-1]*minVec2[-1]/spanVec2[-1])
  
  coef.orgn.scale.1[-1] <- coef.1[-1]/spanVec1[-1]
  coef.orgn.scale.2[-1] <- coef.2[-1]/spanVec2[-1]
  
  coef.orgn.scale.1 <- scalar1(coef.orgn.scale.1)
  coef.orgn.scale.2 <- scalar1(coef.orgn.scale.2)
  
  names(coef.1) <- names(coef.orgn.scale.1)<- colnames(p.data.scale1) 
  names(coef.2) <- names(coef.orgn.scale.2)<- colnames(p.data.scale2)
  
  out <- list(coef.1 = coef.1, 
              coef.orgn.scale.1 =coef.orgn.scale.1,
              coef.2=coef.2,
              coef.orgn.scale.2 =coef.orgn.scale.2)
  return(out)
}


learn_policy <- function(data, model, setting, tau = NULL, regimeClass.stg= NULL, regimeClass.stg1=NULL, regimeClass.stg2=NULL, print.level, pop.size, wait.generations){
  if (setting == "Mean"){
    policy <- learn_policy_mean(data, model, regimeClass.stg1, regimeClass.stg2, print.level, pop.size, wait.generations)
  }else if (setting == "Quan"){
    policy <- learn_policy_quan(data, model, tau, regimeClass.stg1, regimeClass.stg2, print.level, pop.size, wait.generations)
  }else {
    stop("setting error")
  }
  return(policy)
}

