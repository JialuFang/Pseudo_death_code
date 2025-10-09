source("utilities.R")
library(parallel)
library(survival)
library(WINS)
library(dplyr)
library(MASS)

S_prog_fct1 <- function(t,lambda=1/4){exp(-(1+lambda)*t*3/2)}
inv_S_prog_fct1 <- function(x,lambda=1/4){-2/3 * log(x) / (1 + lambda)}

S_death_fct1 <- function(t){exp(-5/8*t)}
inv_S_death_fct1 <- function(x){-8/5*log(x)}

S_prog_fct2 <- function(t){exp(-3/2*t)}
inv_S_prog_fct2 <- function(x){-2/3*log(x)}

S_death_fct2 <- function(t, lambda=2.5){
  lambda*(1/(1+3/4*t))+(1-lambda)*exp(-t/2)
}
inv_S_death_fct2 <- function(x, lambda=2.5){
  target_fct <- function(t, x, lambda){
    lambda*(1/(1+3/4*t))+(1-lambda)*exp(-t/2) -x
  }
  result <- uniroot(target_fct, lower = 0, upper = 1000, x = x, lambda = lambda)
  return(result$root)
}

#> mean(death_list1); mean(death_list2)
#0.63132,0.72017
#> mean(prog_list1); mean(prog_list2)
#0.22547,0.29668

followup_para1 <- 2; followup_para2 <- 2; tau <- 2
run.fn <- function(i){
  set.seed(seeds[i])
  if (i%%20 == 0){
    print(i)
  }
  ndots <- 4
  data1 <- gen_data_and_check_probs(n=samplesize, followup_dist='unif', followup_para=followup_para1, 
                                    inv_S_prog_fct1, S_prog_fct1, S_death_fct1, inv_S_death_fct1, 
                                    tau=tau, theta=1, check_ndots=ndots, A=1)$data
  data2 <- gen_data_and_check_probs(n=samplesize, followup_dist='unif', followup_para=followup_para2, 
                                    inv_S_prog_fct2, S_prog_fct2, S_death_fct2, inv_S_death_fct2, 
                                    tau=tau, theta=1, check_ndots=ndots, A=0)$data
  
  data_all <- rbind(data.frame(data1),data.frame(data2))
  
  #log-rank test
  logrank.res <- logrank.test(data_all)
  
  #wins
  data_use <- data.frame(id=data_all$id, arm=data_all$A, Delta_1=data_all$prog,
                         Delta_2=data_all$death, Y_1=data_all$t.oprog, Y_2=data_all$t.odeath)
  res_tte <- capture.output(win.stat(data = data_use, ep_type = "tte", arm.name = c("1","0"),
                                     priority = c(2:1), alpha = 0.05, digit = 5,
                                     stratum.weight = "unstratified", method = "unadjusted",
                                     pvalue = "two-sided"))
  wins.res <- extract_values(res_tte)
  
  #pseudo death welch Z
  pseu.res.Z <- pseudeath.test.Z(data1, data2, tau, B=200)
  
  #pseudo death permutation
  pseu.res.per <- pseudeath.test.per(data1, data2, tau, B=200)
  
  #pseudo death Wald 
  pseu.res.wald <- pseudeath.test.wald(data1, data2, tau, ndots=4, B=200)
  
  ress_pseu <- c(logrank.res, wins.res, pseu.res.Z, pseu.res.per, pseu.res.wald)
  return(ress_pseu)
}

samplesize_list <- c(50, 100, 200, 500)
for (ii in 1:4){
  nsimu <- 500
  samplesize <- samplesize_list[ii]
  test.res <- list()
  seeds <- c(1:nsimu)
  t <- system.time({
    test.res <- mclapply(1:nsimu, run.fn, mc.cores=40)
  })
  file.name <- paste0("./test_results/","Test_",samplesize,"_simu_use7", ".Rdata")
  save(test.res, file=file.name)
  print(paste("samplesize", samplesize, "finished"))
}

