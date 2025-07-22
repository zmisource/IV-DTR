rm(list = ls())
library(bridgedist)
library(grf)
library(rgenoud)
library(modi)
expit <- function(x) exp(x)/(1 + exp(x))

Val_Data <- function(N, setting, d, tau = NULL){
  if (setting == 'Mean') {
    L1 <- matrix(rnorm(5 * N, mean = 1.5, sd = 1), nrow = N, ncol = 5)
    U1 <- rnorm(N)
    
    nvars.stg1 <- ncol(L1) + 1
    d1 <- as.numeric(d[1:nvars.stg1])
    A1 <- as.numeric(cbind(1,L1) %*% d1 > 0)
    q1 <- 0.2 + L1 %*% c(-0.6, -0.8, 0, 0, 0)
    e1 <- rnorm(N)

    R1 <- 0.5 + L1 %*% c(0.5, 0.8, 0.3, -0.5, 0.7) + A1 * q1 + 1.5*U1 + 1.4*e1
    
    L2 <- matrix(rnorm(5 * N, mean = 2, sd = 1), nrow = N, ncol = 5)
    U2 <- rnorm(N)
    d2 <- as.numeric(d[-c(1:nvars.stg1)])
    A2 <- as.numeric(cbind(1, R1, L2) %*% d2 > 0)
    
    q2 <- L2 %*% c(1, 1, -1, 0, 0) - 0.2 * R1 + 0.4
    
    e2 <- rnorm(N)

    R2 <- 3 + L1 %*% c(0.5, 0.5, 0,0,0) + L2 %*% c(0, 0, 0,-1,1) + A2 * q2 + 2*U2 + 1.2*e2
    out <- mean(R1 + R2)
  } else if (setting == 'Quan'){
    L1 <- matrix(runif(N * 5, min = 0, max = 1), nrow = N, ncol = 5)
    U1 <- rnorm(N)
    nvars.stg1 <- ncol(L1) + 1
    d1 <- as.numeric(d[1:nvars.stg1])
    A1 <- as.numeric(cbind(1,L1) %*% d1 > 0)
    q1 <- 0.2 + L1 %*% c(-0.6, -0.8, 0, 0, 0)
    e1 <- rnorm(N)
    
    R1 <- 0.5 + L1 %*% c(0.5, 0.8, 0.3, -0.5, 0.7) + A1 * q1 + 0.5*U1 + e1
    
    L2 <- matrix(runif(N * 5, min = 0, max = 1), nrow = N, ncol = 5)
    U2 <- rnorm(N)
    d2 <- as.numeric(d[-c(1:nvars.stg1)])
    A2 <- as.numeric(cbind(1, R1, L2) %*% d2 > 0)
    
    q2 <- L2 %*% c(1, 1, -1, 0, 0) - 0.2 * R1 + 0.4
    
    e2 <- rnorm(N)
    
    R2 <- 3 + L1 %*% c(0.5, 0.5, 0.5, 0.5, 0.5) + A2 * q2 + U2 + e2
    out <- quantile(R1+R2, probs = tau)
  }
  return(out)
}

N <- 100000
num_boot <- 500

set.seed(3)
load("Mean_1000.RData")

para <- as.data.frame(theta_ml)
Vmean_1000 <- rep(0,num_boot)

for (i in 1:num_boot) {
  dhat <- colMeans(para, na.rm = TRUE)
  Vmean_1000[i] <- Val_Data(N,'Mean', dhat, tau = NULL)
}


load("Mean_2000.RData")
para <- as.data.frame(theta_ml)

Vmean_2000 <- rep(0,num_boot)
for (i in 1:num_boot) {
  dhat <- colMeans(para, na.rm = TRUE)
  Vmean_2000[i] <- Val_Data(N,'Mean', dhat, tau = NULL)
}

load("Mean_5000.RData")
para <- as.data.frame(theta_ml)

Vmean_5000 <- rep(0,num_boot)
for (i in 1:num_boot) {
  dhat <- colMeans(para, na.rm = TRUE)
  Vmean_5000[i] <- Val_Data(N,'Mean', dhat, tau = NULL)
}


set.seed(3)
load("Quan_1000_0.5.RData")
para <- as.data.frame(theta_ml)

Vquan_1000_0.5 <- rep(0,num_boot)
for (i in 1:num_boot) {
  
  dhat <- colMeans(para, na.rm = TRUE)
  Vquan_1000_0.5[i] <- Val_Data(N,'Quan', dhat, tau = 0.5)
}

load("Quan_2000_0.5.RData")
para <- as.data.frame(theta_ml)

Vquan_2000_0.5 <- rep(0,num_boot)
for (i in 1:num_boot) {
  
  dhat <- colMeans(para, na.rm = TRUE)
  Vquan_2000_0.5[i] <- Val_Data(N,'Quan', dhat, tau = 0.5)
}

load("Quan_5000_0.5.RData")
para <- as.data.frame(theta_ml)

Vquan_5000_0.5 <- rep(0,num_boot)
for (i in 1:num_boot) {
  
  dhat <- colMeans(para, na.rm = TRUE)
  Vquan_5000_0.5[i] <- Val_Data(N,'Quan', dhat, tau = 0.5)
}

print_results <- function(data) {
  for (n in names(data)) {
    mean_values <- round(data[[n]]$mean, 2)
    sd_values <- round(data[[n]]$sd*100, 2)
    cat(n, "&", paste(mean_values, collapse = " & "), "\\\\", "\n")
    cat(" ", "&", paste(paste0("(", sd_values, ")"), collapse = " & "), "\\\\", "\n")
  }
}

results_mean <- list(
  "$n = 1000$" = list(mean = c(mean(Vmean_1000)),
                      sd = c(sd(Vmean_1000))),
  "$n = 2000$" = list(mean = c(mean(Vmean_2000)),
                      sd = c(sd(Vmean_2000))),
  "$n = 5000$" = list(mean = c(mean(Vmean_5000)),
                      sd = c(sd(Vmean_5000)))
)

results_quan <- list(
  "$n = 1000$" = list(mean = c(mean(Vquan_1000_0.5)),
                      sd = c(sd(Vquan_1000_0.5))),
  "$n = 2000$" = list(mean = c(mean(Vquan_2000_0.5)),
                      sd = c(sd(Vquan_2000_0.5))),
  "$n = 5000$" = list(mean = c(mean(Vquan_5000_0.5)),
                      sd = c(sd(Vquan_5000_0.5)))
)


print_results(results_mean)
print_results(results_quan)