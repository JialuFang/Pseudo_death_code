#setwd('tmp/files/total_resp')
source("utilities.R")
library(survival)

# Scenario 1
S_prog_fct <- function(t){exp(-t)}
S_death_fct <- function(t){exp(-t/2)}
inv_S_prog_fct <- function(x){-log(x)}
inv_S_death_fct <- function(x){-2*log(x)}

# Scenario 2
S_prog_fct <- function(t){exp(-t/4*3)}
S_death_fct <- function(t){exp(-t/2)}
inv_S_prog_fct <- function(x){-4/3*log(x)}
inv_S_death_fct <- function(x){-2*log(x)}

# Scenario 4
shape_prog <- 0.5; shape_death <- 0.5
scale_prog <- 1; scale_death <- 2

# Scenario 3
# shape_prog <- 2; shape_death <- 2
# scale_prog <- 1; scale_death <- 2

#S_prog_fct <- function(t){exp(-(t / scale_prog)^shape_prog)}
#S_death_fct <- function(t){exp(-(t / scale_death)^shape_death)}
#inv_S_prog_fct <- function(x){scale_prog * (-log(x))^(1/shape_prog)}
#inv_S_death_fct <- function(x){scale_death * (-log(x))^(1/shape_death)}

curve(S_prog_fct, from = 0, to = 2, col = "blue", lwd = 2, 
      xlab = "Time", ylab = "Survival Probability", 
      main = "Survival Function Curves", ylim = c(0, 1))
curve(S_death_fct, from = 0, to = 2, col = "red", lwd = 2, add = TRUE)


seeds <- c(1:1000)
tau <- 2
samplesize <- c(100, 200, 500, 1000)
for (isample in 1:4){
  kappa2_list <- numeric(1000)
  kappa3_list <- numeric(1000)
  true_kappa2_list <- numeric(1000)
  true_kappa3_list <- numeric(1000)
  for(j in 1:1000){
    set.seed(seeds[j])
    data <- gen.data(n=samplesize[isample], followup_dist = "unif", followup_para=3, inv_S_prog_fct, S_prog_fct, 
                     S_death_fct, inv_S_death_fct, tau, theta=1, A=1)$data
    reskappa <- calc_kappa(data, tau)
      
    true_kappa2_list[j] <- reskappa$true_kappa2
    true_kappa3_list[j] <- reskappa$true_kappa3
    
    kappa2_list[j] <- reskappa$kappa2
    kappa3_list[j] <- reskappa$kappa3
  }
  results <- data.frame(kappa2=kappa2_list, true_kappa2=true_kappa2_list, 
                        kappa3=kappa3_list, true_kappa3=true_kappa3_list)
  file.name <- paste0("./kappa_accuracy_data/","kappa", "_size",  samplesize[isample], "_simu2.Rdata")
  print(file.name)
  save(results, file=file.name)
}
