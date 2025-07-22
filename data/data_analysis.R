rm(list = ls())
library(readxl)
library(grf)
library(rgenoud)
library(modi)
library(MASS)
library(car)

scalar1 <- function(x) {x / sqrt(sum(x^2))}

# 定义转换函数
transformed_data <- function(data){
  label <- as.numeric(as.factor(data$"三级行业类别"))
  L1_original <- cbind(data$"2月交易数", data$"2月交易额")
  A1 <- data$"3月使用营销活动比例"
  R1_original <- data$"3月交易额" - data$"2月交易额"
  
  L2_original <- cbind(data$"5月交易数", data$"5月交易额")
  A2 <- data$"6月使用营销活动比例"
  R2_original <- data$"6月交易额" - data$"5月交易额"
  
  n <- nrow(data)
  Z1 <- rep(0, n)
  Z2 <- rep(0, n)
  tau_l <- 0.95
  tau_u <- 1
  
  for (i in 1:n){
    turnover_2 <- data[i, "2月交易额"]
    turnover_5 <- data[i, "5月交易额"]
    quan_turnover_2_l <- quantile(data$"2月交易额", tau_l)
    quan_turnover_2_u <- quantile(data$"2月交易额", tau_u)
    quan_turnover_5_l <- quantile(data$"5月交易额", tau_l)
    quan_turnover_5_u <- quantile(data$"5月交易额", tau_u)
    if (quan_turnover_2_l <= turnover_2 & turnover_2 <= quan_turnover_2_u){
      Z1[i] <- 1
    }
    if (quan_turnover_5_l <= turnover_5 & turnover_5 <= quan_turnover_5_u){
      Z2[i] <- 1
    }
  }
  
  # 计算最小值并进行偏移
  min_value <- min(cbind(R1_original, R2_original))
  shift <- abs(min_value) + 1
  R1 <- R1_original + shift
  R2 <- R2_original + shift
  
  # 对数转换
  L1 <- log(L1_original)
  R1 <- log(R1)
  L2 <- log(L2_original)
  R2 <- log(R2)
  R <- R1 + R2
  # 将转换后的数据组合成新的数据框
  data <- data.frame(label = label,
                     L1 = I(L1), Z1 = Z1, A1 = A1, R1 = R1, 
                     L2 = I(L2), Z2 = Z2, A2 = A2, R2 = R2, R=R
  )
  
  return(data)
}

mean_policy <- function(data, print.level, pop.size, wait.generations, domains){
  # data <- transformed_data(data)
  # label <- data$label
  L1 <- data$L1
  Z1 <- data$Z1
  A1 <- data$A1
  R1 <- data$R1
  L2 <- data$L2
  Z2 <- data$Z2
  A2 <- data$A2
  R2 <- data$R2
  R <- R1 + R2
  s.tol <- diff(range(R))*1e-03
  
  
  grfZ1 <- probability_forest(X = cbind(L1), Y = as.factor(Z1))
  pi.Z1 <- Z1 * grfZ1$predictions[, 2] + (1 - Z1) * grfZ1$predictions[, 1]
  delta.L1 <- mean(A1[Z1 == 1]) - mean(A1[Z1 == 0])
  grfZ2 <- probability_forest(X = cbind(L1,R1,L2), Y = as.factor(Z2))
  pi.Z2 <- Z2 * grfZ2$predictions[, 2] + (1 - Z2) * grfZ2$predictions[, 1]
  delta.L2 <- mean(A2[Z2 == 1]) - mean(A2[Z2 == 0])
  
  regimeClass.stg1="A1~L1"
  regimeClass.stg2="A2~L1+A1+R1+L2"
  
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
    result <- weighted.mean(R, wts)
    
    return(result)
  }
  
  est <- genoud(fn=V_mean_hat, nvars=nvars.total,
                print.level=print.level, max=TRUE,
                pop.size=pop.size,
                wait.generations=wait.generations,
                gradient.check=FALSE, BFGS=FALSE,
                P1=50, P2=50, P3=50, P4=50, P5=50,
                P6=50, P7=50, P8=50, P9=0,
                default.domains = domains,
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


quan_policy <- function(data, tau, print.level, pop.size, wait.generations, domains){
  # data <- transformed_data(data)
  # label <- data$label
  L1 <- data$L1
  Z1 <- data$Z1
  A1 <- data$A1
  R1 <- data$R1
  L2 <- data$L2
  Z2 <- data$Z2
  A2 <- data$A2
  R2 <- data$R2
  R <- R1 + R2
  s.tol <- diff(range(R))*1e-03
  
  
  grfZ1 <- probability_forest(X = cbind(L1), Y = as.factor(Z1))
  pi.Z1 <- Z1 * grfZ1$predictions[, 2] + (1 - Z1) * grfZ1$predictions[, 1]
  delta.L1 <- mean(A1[Z1 == 1]) - mean(A1[Z1 == 0])
  grfZ2 <- probability_forest(X = cbind(L1,R1,L2), Y = as.factor(Z2))
  pi.Z2 <- Z2 * grfZ2$predictions[, 2] + (1 - Z2) * grfZ2$predictions[, 1]
  delta.L2 <- mean(A2[Z2 == 1]) - mean(A2[Z2 == 0])
  
  regimeClass.stg1="A1~L1"
  regimeClass.stg2="A2~L1+A1+R1+L2"
  
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
  
  est <- genoud(fn=V_quan_hat,tau =tau, nvars=nvars.total,
                print.level=print.level, max=TRUE,
                pop.size=pop.size,
                wait.generations=wait.generations,
                gradient.check=FALSE, BFGS=FALSE,
                P1=50, P2=50, P3=50, P4=50, P5=50,
                P6=50, P7=50, P8=50, P9=0,
                default.domains = domains,
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

val_mean <- function(data, d){
  # data <- transformed_data(data)
  # label <- data$label
  L1 <- data$L1
  Z1 <- data$Z1
  A1 <- data$A1
  R1 <- data$R1
  L2 <- data$L2
  Z2 <- data$Z2
  A2 <- data$A2
  R2 <- data$R2
  R <- R1 + R2
  s.tol <- diff(range(R))*1e-03
  
  
  grfZ1 <- probability_forest(X = cbind(L1), Y = as.factor(Z1))
  pi.Z1 <- Z1 * grfZ1$predictions[, 2] + (1 - Z1) * grfZ1$predictions[, 1]
  delta.L1 <- mean(A1[Z1 == 1]) - mean(A1[Z1 == 0])
  grfZ2 <- probability_forest(X = cbind(L1,R1,L2), Y = as.factor(Z2))
  pi.Z2 <- Z2 * grfZ2$predictions[, 2] + (1 - Z2) * grfZ2$predictions[, 1]
  delta.L2 <- mean(A2[Z2 == 1]) - mean(A2[Z2 == 0])
  
  regimeClass.stg1="A1~L1"
  regimeClass.stg2="A2~L1+A1+R1+L2"
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
  CovSpace_1 <- cbind(1, L1)
  CovSpace_2 <- cbind(1, L1, A1, R1, L2)  
  
  
  
  V_mean_hat <- function(d){
    d1 <- d[1:nvars.stg1]
    A1d <- as.numeric(CovSpace_1%*%d1 > 0)
    d2 <- d[-c(1:nvars.stg1)]
    A2d <- as.numeric(CovSpace_2%*%d2 > 0)
    pic <- pi.Z1 * pi.Z2 * delta.L1 * delta.L2
    c <- as.numeric(A1d == A1) * as.numeric(A2d == A2)*(2 * A1 - 1) * (2 * Z1 - 1) *(2 * A2 - 1) * (2 * Z2 - 1) 
    wts <- c/pic
    result <- weighted.mean(R, wts)
    return(result)
  }
  return(V_mean_hat(d))
}


val_quan <- function(data, d, tau){
  # data <- transformed_data(data)
  # label <- data$label
  L1 <- data$L1
  Z1 <- data$Z1
  A1 <- data$A1
  R1 <- data$R1
  L2 <- data$L2
  Z2 <- data$Z2
  A2 <- data$A2
  R2 <- data$R2
  R <- R1 + R2
  s.tol <- diff(range(R))*1e-03
  
  
  grfZ1 <- probability_forest(X = cbind(L1), Y = as.factor(Z1))
  pi.Z1 <- Z1 * grfZ1$predictions[, 2] + (1 - Z1) * grfZ1$predictions[, 1]
  delta.L1 <- mean(A1[Z1 == 1]) - mean(A1[Z1 == 0])
  grfZ2 <- probability_forest(X = cbind(L1,R1,L2), Y = as.factor(Z2))
  pi.Z2 <- Z2 * grfZ2$predictions[, 2] + (1 - Z2) * grfZ2$predictions[, 1]
  delta.L2 <- mean(A2[Z2 == 1]) - mean(A2[Z2 == 0])
  
  regimeClass.stg1="A1~L1"
  regimeClass.stg2="A2~L1+A1+R1+L2"
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
  CovSpace_1 <- cbind(1, L1)
  CovSpace_2 <- cbind(1, L1,A1, R1, L2)  
  
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
  
  return(V_quan_hat(d, tau))
}

data <- read_excel("~/project/dtr/simulation/real data/results/0101开通营销活动比例_1.xlsx", sheet = 1, col_names = TRUE)
data <- data[,-1]

data <- subset(data,
               data$"2月交易数" > 100 &
                 data$"3月交易数" > 100 &
                 data$"5月交易数" > 100 &
                 data$"6月交易数" > 100 )

k <- 123
set.seed(k)
data <- subset(data, `三级行业类别` %in% c("小吃", "米线/米粉", "炒菜","麻辣烫"))
data <- transformed_data(data)

train_size <- floor(0.7 * nrow(data))

train_indices <- sample(seq_len(nrow(data)), size = train_size)

train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]


out <- mean_policy(train_data,print.level=0, pop.size=3000, wait.generations=8, domains=10)
set.seed(k)
out_0.5 <- quan_policy(train_data,tau = 0.5,print.level=0, pop.size=3000, wait.generations=8, domains=10)
