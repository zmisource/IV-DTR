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

data <- read_excel("/Users/zzzsy/project/dtr/simulation/simulation_final/data/0101开通营销活动比例_1.xlsx", sheet = 1, col_names = TRUE)
data <- data[,-1]

data <- subset(data,
               data$"2月交易数" > 100 &
                 data$"3月交易数" > 100 &
                 data$"5月交易数" > 100 &
                 data$"6月交易数" > 100 )

k <- 123
set.seed(k)
pop.size <- 3000
data <- subset(data, `三级行业类别` %in% c("小吃", "米线/米粉", "炒菜","麻辣烫"))
data <- transformed_data(data)

train_size <- floor(0.7 * nrow(data))

train_indices <- sample(seq_len(nrow(data)), size = train_size)

train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]


out <- mean_policy(train_data,print.level=0, pop.size=3000, wait.generations=8, domains=10)
set.seed(k)
out_0.5 <- quan_policy(train_data,tau = 0.5,print.level=0, pop.size=3000, wait.generations=8, domains=10)

library(quantreg)
library(quantoptr)

get_data <- function(data){
  # data <- transformed_data(data)
  # label <- data$label
  L1 <- data$L1
  A1 <- data$A1
  R1 <- data$R1
  L2 <- data$L2
  A2 <- data$A2
  R2 <- data$R2
  R <- R1 + R2
  return(data.frame(x1 = I(L1), a1 = A1, y1 = R1, x2 = I(L2), a2 = A2, y2 = R2, y = R))
}
train_data_qtr <- get_data(train_data)
s.tol <- diff(range(train_data_qtr$y))*1e-03
set.seed(k)

out1<- TwoStg_Mopt(data = train_data_qtr, regimeClass.stg1 = a1 ~ x1, regimeClass.stg2 = a2 ~ x1+a1 + y1 + x2,
                   moPropen1 = a1 ~ x1, moPropen2 = a2 ~ x1+a1 + y1 + x2,
                   p_level = 0,
                   cl.setup = 1, pop.size = pop.size, it.num = 8, s.tol = s.tol)

out1_0.5<- TwoStg_Qopt(data = train_data_qtr,tau=0.5, regimeClass.stg1 = a1 ~ x1, regimeClass.stg2 = a2 ~ x1+a1 + y1 + x2,
                       moPropen1 = a1 ~ x1, moPropen2 = a2 ~ x1+a1 + y1 + x2,
                       p_level = 0,
                       cl.setup = 1, pop.size = pop.size, it.num = 8, s.tol = s.tol)
library(DynTxRegime)
get_iq_data <- function(data){
  # data <- transformed_data(data)
  # label <- data$label
  L1 <- data$L1
  A1 <- data$A1
  R1 <- data$R1
  L2 <- data$L2
  A2 <- data$A2
  R2 <- data$R2
  R <- R1 + R2
  return(data.frame(X1 = I(L1), A1 = A1, Y1 = R1, X2 = I(L2), A2 = A2, Y2 = R2, Y = R))
}

train_data_iq <- get_iq_data(train_data)

moMain <- buildModelObj(model = ~X1+A1+ Y1 + X2,
                        solver.method = 'lm')

moCont <- buildModelObj(model = ~X1+A1+ Y1 + X2,
                        solver.method = 'lm')

fitSS <- iqLearnSS(moMain = moMain, moCont = moCont,
                   data = train_data_iq, response =train_data_iq$Y,  txName = 'A2')

# main effects model
moMain <- buildModelObj(model = ~X1,
                        solver.method = 'lm')

moCont <- buildModelObj(model = ~X1,
                        solver.method = 'lm')

fitFSC <- iqLearnFSM(moMain = moMain, moCont = moCont,
                     data = train_data_iq, response = fitSS,  txName = 'A1')

d_iq <- c(coef(object = fitFSC)$outcome$Combined[4:6], coef(object = fitSS)$outcome$Combined[8:14])

library(ITRSelect)

train_data_itr <- get_data(train_data)
X1 <- train_data_itr$x1
A1 <- train_data_itr$a1
Y1 <- train_data_itr$y1
X2 <- train_data_itr$x2
A2 <- train_data_itr$a2
Y <- train_data_itr$y

result <- PAL(Y~X1|A1|Y1 + X2|A2, , lambda.list = seq(0.1,1.1,0.1))


d_pal <- c(result$beta1.est, result$beta2.est)
d_iv <- c(out$coef.orgn.scale.1, out$coef.orgn.scale.2)
d_qtr <- c(out1$coef.orgn.scale.1, out1$coef.orgn.scale.2)
d_quan_iv <- c(out_0.5$coef.orgn.scale.1, out_0.5$coef.orgn.scale.2)
d_quan_qtr <- c(out1_0.5$coef.orgn.scale.1, out1_0.5$coef.orgn.scale.2)

tau <- 0.5
# 定义保留小数位数的参数
decimal_places <- 4  # 可以调整此值

# 动态生成格式字符串
format_string <- paste0("%.", decimal_places, "f")

# 计算 val_mean 结果，并格式化
mean_results <- paste(
  sprintf(format_string, val_mean(test_data, d_iv)),
  sprintf(format_string, val_mean(test_data, d_quan_iv)),
  sprintf(format_string, val_mean(test_data, d_qtr)),
  sprintf(format_string, val_mean(test_data, d_quan_qtr)),
  sprintf(format_string, val_mean(test_data, d_iq)),
  sprintf(format_string, val_mean(test_data, d_pal)),
  sep = " & "
)

# 计算 val_quan 结果，并格式化
quan_results <- paste(
  sprintf(format_string, val_quan(test_data, d_iv, tau)),
  sprintf(format_string, val_quan(test_data, d_quan_iv, tau)),
  sprintf(format_string, val_quan(test_data, d_qtr, tau)),
  sprintf(format_string, val_quan(test_data, d_quan_qtr, tau)),
  sprintf(format_string, val_quan(test_data, d_iq, tau = tau)),
  sprintf(format_string, val_quan(test_data, d_pal, tau = tau)),
  sep = " & "
)

# 输出结果
cat("Mean", "&", mean_results, "\\\\", "\n")
cat("Median", "&", quan_results, "\\\\", "\n")
