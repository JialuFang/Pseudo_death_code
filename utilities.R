########generate data
gen.data <- function(n, followup_dist, followup_para, inv_S_prog_fct, S_prog_fct, 
                     S_death_fct, inv_S_death_fct, tau, theta=1, A){
  #args:
  #   n: num of subjects to generate
  #   followup_para: parameter for generating follow-up period 
  #   death_para: parameter for generating death time
  #   prog_para: parameter for generating progression time
  #   tau: an arbitrary time point within the study period at which the survival outcome is of interest
  #Return:
  #   t.prog: actual time from enrollment to progression
  #   prog: binary data that indicates whether progression occurs or not. if progression occurs, prog=1
  #   t.oprog: observed time from enrollment to progression. if no progression, t.eprog=min{t.followup, tau}
  #   rprog: censoring rate for progression event
  #   t.death: actual time from enrollment to death
  #   death: binary data that indicates whether death occurs or not. if death occurs, death=1
  #   t.odeath: observed time from enrollment to death. if no death, t.edeath=min{t.followup, tau}
  #   rdeath: censoring rate for death event
  #   t.followup: follow up time for each subject
  #   t.censor: censoring time for each subject
  
  x_value_fct <- function(u, v, theta){
    a <- (u * v^(theta+1))^(-theta/(theta+1)) - v^{-theta} + 1
    return(1-a**(-1/theta))
  }
  
  final_data <- data.frame()
  
  while (nrow(final_data) < n) {
    n_need <- n - nrow(final_data)  
    
    U.prog <- runif(n_need, 0.01, 0.99)
    t.prog <- sapply(U.prog, inv_S_prog_fct)
    
    U.death <- runif(n_need, 0.01, 0.99)
    v_value <- 1 - S_prog_fct(t.prog)
    x_value <- x_value_fct(U.death, v_value, theta)
    t.death <- sapply(x_value, inv_S_death_fct)
    
    prog <- rep(0, n_need)
    death <- rep(0, n_need)
    
    if (followup_dist == 'unif') {
      t.followup <- runif(n_need, 0, followup_para)
    } else if (followup_dist == 'exp') {
      t.followup <- rexp(n_need, followup_para)
    }
    
    t.censor <- pmin(tau, t.followup)
    
    prog[t.prog <= t.censor] <- 1
    t.oprog <- t.prog
    t.oprog[prog == 0] <- t.censor[prog == 0]
    
    death[t.death <= t.censor] <- 1
    t.odeath <- t.death
    t.odeath[death == 0] <- t.censor[death == 0]
    
    temp_data <- data.frame(prog=prog, t.oprog=t.oprog, t.prog=t.prog,
                            death=death, t.odeath=t.odeath, t.death=t.death, 
                            t.followup=t.followup, t.censor=t.censor, A=A)
    
    temp_data <- subset(temp_data, t.death > t.prog)
    
    final_data <- rbind(final_data, temp_data)
  }
  
  final_data <- final_data[1:n, ]
  final_data$id <- 1:n
  
  nprog = sum(final_data$prog); ndeath = sum(final_data$death)
  
  res <- list(rprog=1-nprog/n, rdeath=1-ndeath/n,
              data=final_data)
  
  return(res)
}

gen.cor.data <- function(n, followup_para, inv_S_prog_fct, S_prog_fct, 
                         S_death_fct, inv_S_death_fct, tau, theta=1, theta_CR, A){
  # by default, followup_dist='unif'
  #args:
  #   n: num of subjects to generate
  #   death_para: parameter for generating death time
  #   prog_para: parameter for generating progression time
  #   tau: an arbitrary time point within the study period at which the survival outcome is of interest
  #   followup_para: parameter for generating follow-up period
  #   follow-up period follows the exponential distribution
  #Return:
  #   t.prog: actual time from enrollment to progression
  #   prog: binary data that indicates whether progression occurs or not. if progression occurs, prog=1
  #   t.oprog: observed time from enrollment to progression. if no progression, t.eprog=min{t.followup, tau}
  #   rprog: censoring rate for progression event
  #   t.death: actual time from enrollment to death
  #   death: binary data that indicates whether death occurs or not. if death occurs, death=1
  #   t.odeath: observed time from enrollment to death. if no death, t.edeath=min{t.followup, tau}
  #   rdeath: censoring rate for death event
  #   t.followup: follow up time for each subject
  #   t.censor: censoring time for each subject
  
  x_value_fct <- function(u, v, theta){
    a <- (u * v^(theta+1))^(-theta/(theta+1)) - v^{-theta} + 1
    return(1-a**(-1/theta))
  }
  
  final_data <- data.frame()
  
  while (nrow(final_data) < n) {
    n_need <- n - nrow(final_data)
    
    clayton_cop <- claytonCopula(param = theta_CR, dim = 2)
    u <- rCopula(n_need, clayton_cop)
    
    U.prog <- u[,1]
    t.prog <- sapply(U.prog, inv_S_prog_fct)
    U.death <- runif(n_need, 0, 1)
    v_value <- 1 - S_prog_fct(t.prog)
    x_value <- x_value_fct(U.death, v_value, theta)
    t.death <- sapply(x_value, inv_S_death_fct)
    
    prog <- rep(0, n_need)
    death <- rep(0, n_need)
    
    followup = rep(0, n_need);
    t.followup = -log(u[,2]) / followup_para
    
    t.censor = pmin(tau, t.followup)
    
    prog[t.prog <= t.censor] <- 1
    t.oprog <- t.prog
    t.oprog[prog == 0] <- t.censor[prog == 0]
    
    death[t.death <= t.censor] <- 1
    t.odeath <- t.death
    t.odeath[death == 0] <- t.censor[death == 0]
    
    temp_data <- data.frame(prog=prog, t.oprog=t.oprog, t.prog=t.prog,
                            death=death, t.odeath=t.odeath, t.death=t.death, 
                            t.followup=t.followup, t.censor=t.censor, A=A)
    
    temp_data <- subset(temp_data, t.death > t.prog)
    
    final_data <- rbind(final_data, temp_data)
  }
  
  final_data <- final_data[1:n, ]
  final_data$id <- 1:n
  
  nprog = sum(final_data$prog); ndeath = sum(final_data$death)
  
  res <- list(rprog=1-nprog/n, rdeath=1-ndeath/n,
              data=final_data)
  
  return(res)
}

gen_data_and_check_probs <- function(n,followup_dist,  followup_para, inv_S_prog_fct, 
                                     S_prog_fct, S_death_fct, inv_S_death_fct, 
                                     tau, theta=1, A, check_ndots=4) {
  ndots <- check_ndots
  tau_list <- seq(tau / ndots, tau, by = tau / ndots)
  
  while (TRUE) {
    res <- gen.data(n = n, followup_dist = followup_dist, followup_para = followup_para, 
                    inv_S_prog_fct, S_prog_fct, S_death_fct, inv_S_death_fct, tau, theta, A)
    data <- res$data
    
    Pr1 <- rep(0, ndots)
    Pr2 <- rep(0, ndots)
    Pr3 <- rep(0, ndots)
    Pr1[ndots] <- mean(data$death == 1)
    Pr2[ndots] <- mean(data$death == 0 & data$prog == 1)
    Pr3[ndots] <- mean(data$death == 0 & data$prog == 0)
    
    for (d in 1:(ndots - 1)) {
      data_test <- trunc_data(data, arb.t = tau_list[d])
      Pr1[d] <- mean(data_test$death == 1)
      Pr2[d] <- mean(data_test$death == 0 & data_test$prog == 1)
      Pr3[d] <- mean(data_test$death == 0 & data_test$prog == 0)
    }
    
    if (all(Pr1 != 0) && all(Pr2 != 0) && all(Pr3 != 0)) {
      return(res)
    }
  }
}

########truncate data
trunc_data <- function(data, arb.t){
  n <- dim(data)[1]
  prog = data$prog
  death = data$death
  
  t.censor <- data$t.odeath
  t.censor = pmin(arb.t, t.censor)
  
  prog[arb.t < data$t.oprog & data$prog == 1]=0
  t.oprog <- data$t.oprog
  t.oprog[prog==0]=t.censor[prog==0];
  
  death[arb.t < data$t.odeath & data$death == 1]=0
  t.odeath <- data$t.odeath
  t.odeath[death==0]=t.censor[death==0];
  
  res <- data.frame(prog=prog, t.oprog=t.oprog,
                     death=death, t.odeath=t.odeath, A = data$A)
  rownames(res) <- rownames(data)
  return(res)
}

########calculate data
calc_kappa <- function(data, tau){
  ##generate the censor event and censor time
  data$censor <- ifelse(data$death == 0, 1, 0)
  data$t.ocensor <- data$t.odeath
  
  ##calculate the cumulative hazard function of C
  surv_obj <- Surv(time = data$t.ocensor, event = data$censor)
  fit <- survfit(surv_obj ~ 1)
  
  #The observed data fall into one of the three categories with probability,
  #(1) Pr1: A patient first experiences progression and then death 
  #(2) Pr2: A patient experiences progression but not death
  #(3) Pr3: A patient experiences neither progression nor death
  Pr1 <- mean(data$death == 1)
  Pr2 <- mean(data$death == 0 & data$prog == 1)
  Pr3 <- mean(data$death == 0 & data$prog == 0)
  
  ##calculate kappa2 and kappa3
  results21 <- numeric(length(fit$time)-1)
  results22 <- numeric(length(fit$time)-1)
  results31 <- numeric(length(fit$time)-1)
  results32 <- numeric(length(fit$time)-1)
  
  cumhaz.diff <- c(fit$cumhaz[1],diff(fit$cumhaz))
  index <- length(fit$time)-1
  tau.eps <- fit$time[index]
  surv.tau.eps <- fit$surv[index]
  surv.diff <- c(fit$surv[1]-1, diff(fit$surv))
  
  for (i in 1:(length(fit$time)-1)) {
    t <- fit$time[i]
    results21[i] <- mean(t < data$t.odeath & data$t.oprog <= t) * cumhaz.diff[i]
    results22[i] <- mean(tau.eps < data$t.odeath & data$t.oprog <= t) * surv.diff[i]
  }
  
  for (i in 1:(length(fit$time)-1)) {
    t <- fit$time[i]
    results31[i] <- mean(t < data$t.odeath & data$t.oprog > t) * cumhaz.diff[i]
    results32[i] <- mean(tau.eps < data$t.odeath & data$t.oprog > t) * surv.diff[i]
  }
  
  kappa2 <- (sum(results21)+sum(results22)/surv.tau.eps)/Pr2
  kappa3 <- (sum(results31)+sum(results32)/surv.tau.eps)/Pr3
  
  true_kappa2 <- sum(data$t.death < tau & data$death == 0 & data$prog == 1)/sum(data$death == 0 & data$prog == 1)
  true_kappa3 <- sum(data$t.death < tau & data$death == 0 & data$prog == 0)/sum(data$death == 0 & data$prog == 0)
  
  return(list(kappa2=kappa2, true_kappa2=true_kappa2,
              kappa3=kappa3, true_kappa3=true_kappa3))
}

###########Win Statistics#########
logrank.test <- function(data){
  surv_object_prog <- Surv(time = data$t.oprog, event = data$prog)
  log_rank_test_prog <- survdiff(surv_object_prog ~ data$A)
  p.prog <- 1 - pchisq(log_rank_test_prog$chisq, df = 1)
  surv_object_death <- Surv(time = data$t.odeath, event = data$death)
  log_rank_test_death <- survdiff(surv_object_death ~ data$A)
  p.death <- 1 - pchisq(log_rank_test_death$chisq, df = 1)
  return(list(log_rank_prog=log_rank_test_prog$chisq, log_rank_prog_p=p.prog, 
              log_rank_death=log_rank_test_death$chisq, log_rank_death_p=p.death))
}

extract_values <- function(res_tte) {
  win_ratio_pattern <- "Win Ratio :\\s*(-?[0-9\\.]+)"
  net_benefit_pattern <- "Net Benefit :\\s*(-?[0-9\\.]+)"
  win_odds_pattern <- "Win Odds :\\s*(-?[0-9\\.]+)"
  p_value_pattern <- "two-sided p-value is:\\s*(-?[<0-9\\. ]+)"
  
  win_ratio <- NULL
  win_ratio_p <- NULL
  net_benefit <- NULL
  net_benefit_p <- NULL
  win_odds <- NULL
  win_odds_p <- NULL
  count <- 0
  
  for (line in res_tte) {
    if (grepl(win_ratio_pattern, line)) {
      win_ratio <- sub(win_ratio_pattern, "\\1", line)
      count <- 1
    }
    if (grepl(net_benefit_pattern, line)) {
      net_benefit <- sub(net_benefit_pattern, "\\1", line)
      count <- 2
    }
    if (grepl(win_odds_pattern, line)) {
      win_odds <- sub(win_odds_pattern, "\\1", line)
      count <- 3
    }
    if (grepl(p_value_pattern, line)) {
      p_value <- sub(p_value_pattern, "\\1", line)
      if (count == 1) {
        win_ratio_p <- p_value
      } else if (count == 2) {
        net_benefit_p <- p_value
      } else if (count == 3) {
        win_odds_p <- p_value
      }
    }
  }
  
  return(list(
    win_ratio = as.numeric(win_ratio), win_ratio_p = win_ratio_p,
    net_benefit = as.numeric(net_benefit), net_benefit_p = net_benefit_p,
    win_odds = as.numeric(win_odds), win_odds_p = win_odds_p
  ))
}

###########Bootstrap (transform Welch's t into Z)#########
pseudeath.mu.var <- function(data, tau){
  cal.mean <- function(data, kappa1=1, kappa2, kappa3){
    data$NT1 = data$death == 1
    data$NT0_NS1 = data$death == 0 & data$prog == 1
    data$NT0_NS0 = data$death == 0 & data$prog == 0
    
    group_means <- data %>%
      summarise(
        mean_NT1 = mean(NT1),
        mean_NT0_NS1 = mean(NT0_NS1),
        mean_NT0_NS0 = mean(NT0_NS0)
      )
    
    use_means <- as.vector(group_means)
    mu <- kappa1*use_means$mean_NT1+kappa2*use_means$mean_NT0_NS1+kappa3*use_means$mean_NT0_NS0
    
    return(mu)
  }
  
  cal.var <- function(data, kappa1=1, kappa2, kappa3){
    data$NT1 = data$death == 1
    data$NT0_NS1 = data$death == 0 & data$prog == 1
    data$NT0_NS0 = data$death == 0 & data$prog == 0
    
    group_vars <- data %>%
      summarise(
        var_NT1 = var(NT1),
        var_NT0_NS1 = var(NT0_NS1),
        var_NT0_NS0 = var(NT0_NS0)
      )
    
    group_covs <- data %>%
      summarise(
        cov_NT1_NT0_NS1 = cov(NT1, NT0_NS1),
        cov_NT0_NS1_NT0_NS0 = cov(NT0_NS1, NT0_NS0),
        cov_NT1_NT0_NS0 = cov(NT1, NT0_NS0)
      )
    
    use_vars <- as.vector(group_vars)
    use_covs <- as.vector(group_covs)
    sigma2 <- kappa1^2 * use_vars$var_NT1 + kappa2^2 * use_vars$var_NT0_NS1 +
      kappa3^2 * use_vars$var_NT0_NS0 + 2*kappa1*kappa2*use_covs$cov_NT1_NT0_NS1 +
      2*kappa2*kappa3*use_covs$cov_NT0_NS1_NT0_NS0 + 2*kappa1*kappa3*use_covs$cov_NT1_NT0_NS0
    
    return(sigma2)
  }
  
  res_kappa <- calc_kappa(data, tau)
  kappa2 <- res_kappa$kappa2; kappa3 <- res_kappa$kappa3
  mu <- cal.mean(data, kappa1=1, kappa2, kappa3)
  sigma2 <- cal.var(data, kappa1=1, kappa2, kappa3)
  return(list(mu=mu, sigma2=sigma2, kappa2=kappa2, kappa3=kappa3))
}

pseudeath.test.Z <- function(data1, data2, tau, B){
  data_all <- rbind(data.frame(data1),data.frame(data2))
  res_t1 <- pseudeath.mu.var(data1, tau)
  res_t2 <- pseudeath.mu.var(data2, tau)
  mu_t1 <- res_t1$mu; mu_t2 <- res_t2$mu
  
  mu1_list <- numeric(B)
  mu2_list <- numeric(B)
  valid_count <- 0
  
  n <- dim(data1)[1]
  while (valid_count < B) {
    indices <- sample(1:(2*n), size = 2*n, replace = TRUE)
    bootstrap_sample <- data_all[indices, ]
    bootstrap_sample1 <- subset(bootstrap_sample, A == 1)
    bootstrap_sample2 <- subset(bootstrap_sample, A == 0)
    mu1 <- pseudeath.mu.var(bootstrap_sample1, tau)$mu
    mu2 <- pseudeath.mu.var(bootstrap_sample2, tau)$mu

    
    if (!is.nan(mu1) && !is.nan(mu2) && length(mu1) > 0 && length(mu2) > 0) {
      valid_count <- valid_count + 1
      mu1_list[valid_count] <- mu1
      mu2_list[valid_count] <- mu2
    } 
  }
  sigma2_dif <- var(mu1_list-mu2_list)
  Z_score <- (mu_t1-mu_t2)/sqrt(sigma2_dif)
  pvalue_Z <- 2 * (1 - pnorm(abs(Z_score)))
  return(list(pseudo_death_B=Z_score, pseudo_death_B_p=pvalue_Z))
}

###########permutation test#########
pseudeath.test.per <- function(data1, data2, tau, B){
  data_all <- rbind(data.frame(data1),data.frame(data2))
  res_t1 <- pseudeath.mu.var(data1, tau)
  res_t2 <- pseudeath.mu.var(data2, tau)
  mu_t1 <- res_t1$mu; mu_t2 <- res_t2$mu
  mu_dif <- mu_t1-mu_t2
  
  mu1_list <- numeric(B)
  mu2_list <- numeric(B)
  valid_count <- 0
  
  n <- dim(data1)[1]
  while (valid_count < B) {
    
    bootstrap_sample <- data_all
    bootstrap_sample$A <- rbinom(n = nrow(data_all), size = 1, prob = 0.5)
    
    bootstrap_sample1 <- subset(bootstrap_sample, A == 1)
    bootstrap_sample2 <- subset(bootstrap_sample, A == 0)
    
    mu1 <- pseudeath.mu.var(bootstrap_sample1, tau)$mu
    mu2 <- pseudeath.mu.var(bootstrap_sample2, tau)$mu
    
    if (!is.nan(mu1) && !is.nan(mu2) && length(mu1) > 0 && length(mu2) > 0) {
      valid_count <- valid_count + 1
      mu1_list[valid_count] <- mu1
      mu2_list[valid_count] <- mu2
    } 
  }
  mu_dif_list <- mu1_list - mu2_list
  pvalue_per <- mean(abs(mu_dif_list) > abs(mu_dif))
  return(list(pseudo_death_per_p=pvalue_per))
}

###########Wald test#########
pseudeath.mu.single <- function(data, tau){
  cal.mean <- function(data, kappa1=1, kappa2, kappa3){
    data$NT1 = data$death == 1
    data$NT0_NS1 = data$death == 0 & data$prog == 1
    data$NT0_NS0 = data$death == 0 & data$prog == 0
    
    group_means <- data %>%
      summarise(
        mean_NT1 = mean(NT1),
        mean_NT0_NS1 = mean(NT0_NS1),
        mean_NT0_NS0 = mean(NT0_NS0)
      )
    
    use_means <- as.vector(group_means)
    mu <- kappa1*use_means$mean_NT1+kappa2*use_means$mean_NT0_NS1+kappa3*use_means$mean_NT0_NS0
    
    return(mu)
  }
  
  res_kappa <- calc_kappa(data, tau)
  kappa2 <- res_kappa$kappa2; kappa3 <- res_kappa$kappa3
  mu <- cal.mean(data, kappa1=1, kappa2, kappa3)
  return(mu)
}

pseudeath.mu.vec <- function(data, tau, ndots=4){
  tau_list <- seq(tau/ndots, tau, by = tau/ndots)
  data_list <- list()
  for (i in 1:(ndots-1)){
    data_tru <- trunc_data(data, arb.t=tau_list[i])
    data_list[[paste0("data", i)]] <- data_tru
  }
  data_list[[paste0("data", ndots)]] <- data
  
  mu.vec <- rep(0, ndots)
  for (i in 1:ndots){
    res <- pseudeath.mu.single(data_list[[paste0("data", i)]], tau_list[i])
    if (length(res) == 0 || is.null(res) || is.na(res)) {
      mu.vec[i] <- NA  
    } else {
      mu.vec[i] <- res
    }
  }
  
  return(mu.vec)
}

pseudeath.test.wald <- function(data1, data2, tau, ndots, B){
  data_all <- rbind(data.frame(data1),data.frame(data2))
  mean1 <- pseudeath.mu.vec(data1, tau, ndots)
  mean2 <- pseudeath.mu.vec(data2, tau, ndots)
  mean_diff <- mean1 - mean2
  
  delta_boot_list <- list()
  valid_count <- 0
  n <- nrow(data1)
  
  while (valid_count < B) {
    indices <- sample(1:(2*n), size = 2*n, replace = TRUE)
    bootstrap_sample <- data_all[indices, ]
    bootstrap_sample1 <- subset(bootstrap_sample, A == 1)
    bootstrap_sample2 <- subset(bootstrap_sample, A == 0)
    
    mu1 <- pseudeath.mu.vec(bootstrap_sample1, tau, ndots)
    mu2 <- pseudeath.mu.vec(bootstrap_sample2, tau, ndots)
    
    mu_dif <- mu1-mu2
    if (all(!is.na(mu_dif))) {
      valid_count <- valid_count + 1
      delta_boot_list[[valid_count]] <- mu_dif
    }
  }
  
  delta_mat <- do.call(rbind, delta_boot_list)
  cov_pooled <- cov(delta_mat)
  
  # Compute the test statistic and p value
  chi_stat <- t(mean_diff) %*% ginv(cov_pooled) %*% mean_diff
  p_value <- 1 - pchisq(chi_stat, df = length(mean_diff))
  
  stat_name <- paste0("pseudo_death_wal_n", ndots)
  pvalue_name <- paste0("pseudo_death_wal_n", ndots, "_p")
  
  result <- list()
  result[[stat_name]] <- chi_stat
  result[[pvalue_name]] <- p_value
  
  return(result)
}