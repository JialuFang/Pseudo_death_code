source("utilities.R")
library(copula)
library(survival)

S_prog_fct <- function(t){exp(-t)}
S_death_fct <- function(t){exp(-t/2)}
inv_S_prog_fct <- function(x){-log(x)}
inv_S_death_fct <- function(x){-2*log(x)}


seeds <- c(1:1000)
tau <- 2; 
samplesize <- c(100, 200, 500, 1000)
for (isample in 1:4){
  print(isample)
  kappa2_list <- numeric(1000)
  kappa3_list <- numeric(1000)
  true_kappa2_list <- numeric(1000)
  true_kappa3_list <- numeric(1000)
  rprog_list <- numeric(1000); rdeath_list <- numeric(1000)
  for(j in 1:1000){
    set.seed(seeds[j])
    res <- gen.cor.data(n=samplesize[isample], followup_para=1/4, inv_S_prog_fct, S_prog_fct, 
                     S_death_fct, inv_S_death_fct, tau, theta=1, theta_CR=0.1, A=1)
    rprog_list[j] <- res$rprog; rdeath_list[j] <- res$rdeath
    data <- res$data
    reskappa <- calc_kappa(data, tau)
    
    true_kappa2_list[j] <- reskappa$true_kappa2
    true_kappa3_list[j] <- reskappa$true_kappa3
    
    kappa2_list[j] <- reskappa$kappa2
    kappa3_list[j] <- reskappa$kappa3
  }
  results <- data.frame(kappa2=kappa2_list, true_kappa2=true_kappa2_list, 
                        kappa3=kappa3_list, true_kappa3=true_kappa3_list)
  file.name <- paste0("./test_results/","cor_kappa", "_size",  samplesize[isample], "_simu_1.Rdata")
  save(results, file=file.name)
}

