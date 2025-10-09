alpha <- 0.05
samplesize_list <- c(50, 100, 200, 500)

pvalue_df <- data.frame(matrix(ncol = 9, nrow = 4))
colnames(pvalue_df) <- c("sample_size", "log_rank_prog", "log_rank_death", "win_ratio", 
                         "net_benefit", "win_odds", "pseudo_death_Z", "pseudo_death_per",
                         "pseudo_death_wald")

theta_CR <- 1
for (j in c(1:4)){
  samplesize <- samplesize_list[j]
  #file.name <- paste0("./test_results/","H0_Test_",samplesize,"_simu_cor",theta_CR*10,".Rdata")
  file.name <- paste0("./typeIerror_data/","H0_Test_",samplesize,"_simu1", ".Rdata")
  print(file.name)
  load(file.name)
  
  log_rank_prog_p_values <- numeric(length(test.res))
  log_rank_death_p_values <- numeric(length(test.res))
  win_ratio_p_values <- numeric(length(test.res))
  net_benefit_p_values <- numeric(length(test.res))
  win_odds_p_values <- numeric(length(test.res))
  pseudo_death_Z_p_values <- numeric(length(test.res))
  pseudo_death_per_p_values <- numeric(length(test.res))
  pseudo_death_wald_p_values <- numeric(length(test.res))
  
  for (i in seq_along(test.res)) {
    log_rank_prog_p_values[i] <- test.res[[i]]$log_rank_prog_p
    log_rank_death_p_values[i] <- test.res[[i]]$log_rank_death_p
    win_ratio_p <- test.res[[i]]$win_ratio_p
    if (grepl("<", win_ratio_p)) {
      win_ratio_p_values[i] <- 0
    } else {
      win_ratio_p_values[i] <- as.numeric(win_ratio_p)
    }
    net_benefit_p <- test.res[[i]]$net_benefit_p
    if (grepl("<", net_benefit_p)) {
      net_benefit_p_values[i] <- 0
    } else {
      net_benefit_p_values[i] <- as.numeric(net_benefit_p)
    }
    win_odds_p <- test.res[[i]]$win_odds_p
    if (grepl("<", win_odds_p)) {
      win_odds_p_values[i] <- 0
    } else {
      win_odds_p_values[i] <- as.numeric(win_odds_p)
    }
    pseudo_death_Z_p_values[i] <- test.res[[i]]$pseudo_death_B_p
    pseudo_death_per_p_values[i] <- test.res[[i]]$pseudo_death_per_p
    pseudo_death_wald_p_values[i] <- test.res[[i]]$pseudo_death_wal_n4_p
  }
  
  pvalue_df[j,1] <- samplesize
  pvalue_df[j,2] <- mean(log_rank_prog_p_values < alpha)
  pvalue_df[j,3] <- mean(log_rank_death_p_values < alpha)
  pvalue_df[j,4] <- mean(win_ratio_p_values < alpha)
  pvalue_df[j,5] <- mean(net_benefit_p_values < alpha)
  pvalue_df[j,6] <- mean(win_odds_p_values < alpha)
  pvalue_df[j,7] <-  mean(pseudo_death_Z_p_values < alpha)
  pvalue_df[j,8] <-  mean(pseudo_death_per_p_values < alpha)
  pvalue_df[j,9] <-  mean(pseudo_death_wald_p_values < alpha)
} 
