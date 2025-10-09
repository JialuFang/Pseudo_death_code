source("utilities.R")
library(parallel)
library(survival)
library(WINS)
library(dplyr)

S_prog_fct <- function(t){exp(-t/4*3)}
S_death_fct <- function(t){exp(-t/2)}
inv_S_prog_fct <- function(x){-3/4*log(x)}
inv_S_death_fct <- function(x){-2*log(x)}

run.fn <- function(i){
  set.seed(seeds[i])
  if (i%%20 == 0){
    print(i)
  }
  ndots <- 4
  data1 <- gen_data_and_check_probs(n=samplesize, followup_dist='unif', followup_para=3, inv_S_prog_fct, 
                                    S_prog_fct, S_death_fct, inv_S_death_fct, 
                                    tau=2, theta=1, check_ndots=ndots, A=1)$data
  data2 <- gen_data_and_check_probs(n=samplesize, followup_dist='unif', followup_para=3, inv_S_prog_fct, 
                                    S_prog_fct, S_death_fct, inv_S_death_fct, 
                                    tau=2, theta=1, check_ndots=ndots, A=0)$data
  
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
  
  #pseudo death Z
  pseu.res.Z <- pseudeath.test.Z(data1, data2, tau, B=200)
  
  #pseudo death per
  pseu.res.per <- pseudeath.test.per(data1, data2, tau, B=200)
  
  # pseudo death wald
  pseu.res.wald <- pseudeath.test.wald(data1, data2, tau, ndots=4, B=200)
  
  ress_pseu <- c(logrank.res, wins.res, pseu.res.Z, pseu.res.per, pseu.res.wald)
  return(ress_pseu)
}

samplesize_list <- c(50, 100, 200, 500)
for (ii in 1:4){
  nsimu <- 1000
  tau <- 2; samplesize <- samplesize_list[ii]
  test.res <- list()
  seeds <- c(1:nsimu)
  t <- system.time({
    test.res <- mclapply(1:nsimu, run.fn, mc.cores=40)
  })
  file.name <- paste0("./test_results/","H0_Test_",samplesize,"_simu2", ".Rdata")
  save(test.res, file=file.name)
  print(paste("samplesize", samplesize, "finished"))
}



