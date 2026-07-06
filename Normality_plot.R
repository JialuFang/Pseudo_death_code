summarize_bootstrap_validation <- function(validation.df) {
  mc_mean_delta <- mean(validation.df$delta_hat, na.rm = TRUE)
  mc_var_delta  <- var(validation.df$delta_hat, na.rm = TRUE)
  mean_vb_hat   <- mean(validation.df$vb_hat, na.rm = TRUE)
  median_vb_hat <- median(validation.df$vb_hat, na.rm = TRUE)
  ratio_mean_mc <- mean_vb_hat / mc_var_delta
  diff_mean_mc  <- mean_vb_hat - mc_var_delta
  
  out <- data.frame(
    MC_mean_delta = mc_mean_delta,
    MC_var_delta = mc_var_delta,
    Mean_bootstrap_variance = mean_vb_hat,
    Median_bootstrap_variance = median_vb_hat,
    Ratio_meanBoot_to_MCvar = ratio_mean_mc,
    Diff_meanBoot_minus_MCvar = diff_mean_mc
  )
  
  return(out)
}

summarize_normal_approximation <- function(validation.df) {
  z <- validation.df$z_boot
  z <- z[is.finite(z)]
  
  emp_q <- quantile(z, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), na.rm = TRUE)
  theo_q <- qnorm(c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
  
  out1 <- data.frame(
    Mean_Z = mean(z),
    SD_Z = sd(z),
    Rej_5pct = mean(abs(z) > qnorm(0.975)),
    Rej_10pct = mean(abs(z) > qnorm(0.95))
  )
  
  out2 <- data.frame(
    Prob = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975),
    Empirical = as.numeric(emp_q),
    Theoretical = as.numeric(theo_q),
    Difference = as.numeric(emp_q - theo_q)
  )
  
  return(list(summary = out1, quantiles = out2))
}

samplesize_list <- c(50, 100, 200, 500)
res_table <- vector("list", length(samplesize_list))
Z_table <- vector("list", length(samplesize_list))
for (ii in 1:4) {
  samplesize <- samplesize_list[ii]
  file.name <- paste0("./Validation_bootstrap/", "Validation_bootstrap_n", samplesize, "_simu4.Rdata")
  load(file.name)
  tmp <- summarize_bootstrap_validation(validation.df)
  tmp <- cbind(samplesize = samplesize, tmp)
  res_table[[ii]] <- tmp
  
  Ztmp <- summarize_normal_approximation(validation.df)$summary
  Q5<- summarize_normal_approximation(validation.df)$quantiles[2,2]
  Q95<- summarize_normal_approximation(validation.df)$quantiles[6,2]
  Ztmp <- cbind(Ztmp, Q5, Q95)
  Ztmp <- cbind(samplesize = samplesize, Ztmp)
  Z_table[[ii]] <- Ztmp
}
res_table <- bind_rows(res_table)
Z_table <- bind_rows(Z_table)
