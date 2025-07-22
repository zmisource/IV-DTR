source("iv_dtr.R")
N <- c(1000, 2000, 5000)
num_boot <- 500
model <- "RF"
setting <- "Mean"

print.level <- 0
pop.size <- 8000
wait.generations <- 8

for (i in N){
  set.seed(123)
  data <- Gen_Data(N)
  sample_output <- learn_policy(data, model, setting, tau = NULL, regimeClass.stg= NULL, regimeClass.stg1=NULL, regimeClass.stg2=NULL, print.level, pop.size, wait.generations)
  out_length <- length(c(sample_output$coef.orgn.scale.1, sample_output$coef.orgn.scale.2))
  
  theta_ml <- matrix(NA, nrow = num_boot, ncol = out_length)
  colnames(theta_ml) <- c(names(sample_output$coef.orgn.scale.1), names(sample_output$coef.orgn.scale.2))
  for (i in 1:num_boot) {
    cat(i, as.character(Sys.time()), "\n")
    data <- Gen_Data(N)
    res <- learn_policy(data, model, setting, tau = NULL, regimeClass.stg= NULL, regimeClass.stg1=NULL, regimeClass.stg2=NULL, print.level, pop.size, wait.generations)
    theta_ml[i, ] <- c(res$coef.orgn.scale.1, res$coef.orgn.scale.2)
  }
  
  file_name <- sprintf("%s_%d.RData", setting, N)
  save(theta_ml, file = file_name)
}



setting <- "Quan"
tau <- 0.5

for (i in N){
  set.seed(123)
  data <- Gen_Data(N)
  sample_output <- learn_policy(data, model, setting, tau = tau, regimeClass.stg= NULL, regimeClass.stg1=NULL, regimeClass.stg2=NULL, print.level, pop.size, wait.generations)
  out_length <- length(c(sample_output$coef.orgn.scale.1, sample_output$coef.orgn.scale.2))
  
  theta_ml <- matrix(NA, nrow = num_boot, ncol = out_length)
  colnames(theta_ml) <- c(names(sample_output$coef.orgn.scale.1), names(sample_output$coef.orgn.scale.2))
  for (i in 1:num_boot) {
    cat(i, as.character(Sys.time()), "\n")
    data <- Gen_Data(N)
    res <- learn_policy(data, model, setting, tau = tau, regimeClass.stg= NULL, regimeClass.stg1=NULL, regimeClass.stg2=NULL, print.level, pop.size, wait.generations)
    theta_ml[i, ] <- c(res$coef.orgn.scale.1, res$coef.orgn.scale.2)
  }
  
  file_name <- sprintf("%s_%d_%.2f.RData", setting, N,tau)
  save(theta_ml, file = file_name)
}

