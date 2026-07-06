source("utilities.R")
library(parallel)
library(survival)
library(WINS)
library(dplyr)

pseudeath.boot.diagnostic.strat <- function(data1, data2, tau, B = 200) {
  # observed estimates
  res_t1 <- pseudeath.mu.var(data1, tau)
  res_t2 <- pseudeath.mu.var(data2, tau)
  
  mu_t1 <- res_t1$mu
  mu_t2 <- res_t2$mu
  delta_hat <- mu_t1 - mu_t2
  
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  
  delta_boot <- numeric(B)
  valid_count <- 0
  
  while (valid_count < B) {
    idx1 <- sample(1:n1, size = n1, replace = TRUE)
    idx2 <- sample(1:n2, size = n2, replace = TRUE)
    
    boot1 <- data1[idx1, ]
    boot2 <- data2[idx2, ]
    
    mu1 <- tryCatch(pseudeath.mu.var(boot1, tau)$mu, error = function(e) NaN)
    mu2 <- tryCatch(pseudeath.mu.var(boot2, tau)$mu, error = function(e) NaN)
    
    if (!is.nan(mu1) && !is.nan(mu2) && length(mu1) > 0 && length(mu2) > 0) {
      valid_count <- valid_count + 1
      delta_boot[valid_count] <- mu1 - mu2
    }
  }
  
  vb_hat <- var(delta_boot)
  se_boot <- sqrt(vb_hat)
  z_boot <- delta_hat / se_boot
  p_boot <- 2 * (1 - pnorm(abs(z_boot)))
  
  return(list(
    mu1_hat = mu_t1,
    mu0_hat = mu_t2,
    delta_hat = delta_hat,
    vb_hat = vb_hat,
    se_boot = se_boot,
    z_boot = z_boot,
    p_boot = p_boot
  ))
}

run.validation.fn <- function(i) {
  set.seed(seeds[i])
  
  if (i %% 100 == 0) {
    print(paste("replicate", i))
  }
  
  ndots <- 4
  
  data1 <- gen_data_and_check_probs(
    n = samplesize,
    followup_dist = "unif",
    followup_para = 3,
    inv_S_prog_fct,
    S_prog_fct,
    S_death_fct,
    inv_S_death_fct,
    tau = tau,
    theta = 1,
    check_ndots = ndots,
    A = 1
  )$data
  
  data2 <- gen_data_and_check_probs(
    n = samplesize,
    followup_dist = "unif",
    followup_para = 3,
    inv_S_prog_fct,
    S_prog_fct,
    S_death_fct,
    inv_S_death_fct,
    tau = tau,
    theta = 1,
    check_ndots = ndots,
    A = 0
  )$data
  
  # choose one of the two functions
  diag_res <- pseudeath.boot.diagnostic.strat(data1, data2, tau = tau, B = Bboot)
  
  return(data.frame(
    rep = i,
    mu1_hat = diag_res$mu1_hat,
    mu0_hat = diag_res$mu0_hat,
    delta_hat = diag_res$delta_hat,
    vb_hat = diag_res$vb_hat,
    se_boot = diag_res$se_boot,
    z_boot = diag_res$z_boot,
    p_boot = diag_res$p_boot
  ))
}

S_prog_fct <- function(t) { exp(-t) }
S_death_fct <- function(t) { exp(-t/2) }
inv_S_prog_fct <- function(x) { -log(x) }
inv_S_death_fct <- function(x) { -2*log(x) }

samplesize_list <- c(50, 100, 200, 500)

for (ii in 1:4) {
  nsimu <- 1000
  Bboot <- 200
  tau <- 2
  samplesize <- samplesize_list[ii]
  seeds <- 1:nsimu
  
  cat("Start validation simulation for n =", samplesize, "\n")
  
  t_run <- system.time({
    validation.res <- mclapply(1:nsimu, run.validation.fn, mc.cores = 1)
  })
  
  validation.df <- bind_rows(validation.res)
  
  file.name <- paste0("./test_results/", "Validation_bootstrap_n", samplesize, "_simu1.Rdata")
  save(validation.df, t_run, file = file.name)
  
  cat("Finished validation simulation for n =", samplesize, "\n")
}
