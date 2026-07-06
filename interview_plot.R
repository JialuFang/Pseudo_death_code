library(ggplot2)
library(gridExtra)
alpha <- 0.05
samplesize_list <- c(50, 100, 200, 500)

pvalue_df <- data.frame(matrix(ncol = 7, nrow = 4))
colnames(pvalue_df) <- c("sample_size", "log_rank_prog", "log_rank_death", "win_statistic", 
                         "pseudo_death_Z", "pseudo_death_per",
                         "pseudo_death_wald")

res.ls <- list()
for (k in c(1:8)){
  print(k)
  for (j in c(1:4)){
    samplesize <- samplesize_list[j]
    file.name <- paste0("./power_data/","Test_",samplesize,"_simu_use",k, ".Rdata")
    load(file.name)
    print(file.name)
    
    log_rank_prog_p_values <- numeric(length(test.res))
    log_rank_death_p_values <- numeric(length(test.res))
    win_ratio_p_values <- numeric(length(test.res))
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
      pseudo_death_Z_p_values[i] <- test.res[[i]]$pseudo_death_B_p
      pseudo_death_per_p_values[i] <- test.res[[i]]$pseudo_death_per_p
      pseudo_death_wald_p_values[i] <- test.res[[i]]$pseudo_death_wal_n4_p
    }
    
    pvalue_df[j,1] <- samplesize
    pvalue_df[j,2] <- mean(log_rank_prog_p_values < alpha)
    pvalue_df[j,3] <- mean(log_rank_death_p_values < alpha)
    pvalue_df[j,4] <- mean(win_ratio_p_values < alpha)
    pvalue_df[j,5] <-  mean(pseudo_death_Z_p_values < alpha)
    pvalue_df[j,6] <-  mean(pseudo_death_per_p_values < alpha)
    pvalue_df[j,7] <-  mean(pseudo_death_wald_p_values < alpha)
  } 
  res.ls[[k]] <- pvalue_df
}

plot_list <- list()

title<- c("Proportional hazards", "Crossing hazards", "Crossing survival functions", 
          "Early difference", "Middle difference", "Late difference",
          "Counteracting effect", "Absent effect")

pdf(file = "plots/power_seminar.pdf", width = 12, height = 4.5)

pal <- setNames(
  c("#D62728", # 红 red
    "#FF7F0E", # 橙 orange
    "#FFD700", # 黄 yellow (gold)
    "#2CA02C", # 绿 green
    "#1F77B4", # 蓝 blue
    "#9467BD"  # 紫 purple
  ),
  m.names
)

for (i in 1:8){
  grp.names <- samplesize_list
  m.names <-   c("Log-rank TTP", "Log-rank OS", "Win statistics", "Pseudo-death Z", "Pseudo-death\n Permutation",
                 "Pseudo-death\n Wald")
  g.var <- rep(grp.names, each=length(m.names))
  m.var <- rep(m.names, times=length(grp.names))
  data <- res.ls[[i]][, -1]
  v.var <- as.vector(t(data))
  data <- data.frame(g=factor(g.var, levels=grp.names), m=factor(m.var, levels=m.names), v=v.var)
  
  p <- ggplot(data = data, mapping = aes(x = g, y = v, color = m, group = m)) +
    geom_point(size = 1.2) +
    geom_line() +
    scale_color_manual(values = pal, breaks = m.names, limits = m.names) +
    theme(legend.key.size = unit(2, 'mm'), 
          legend.text = element_text(size = 6),
          axis.title.x  = element_text(size = 10),
          axis.title.y  = element_text(size = 10),
          legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
          legend.position = 'right',
          plot.title = element_text(size = 10,hjust=0.5)) +
    xlab("n") + ylab("Power") + 
    guides(color=guide_legend(title = "")) + 
    ggtitle(title[i])
  
  plot_list[[i]] <- p
}

grid.arrange(grobs = plot_list, nrow = 2, ncol = 4)

dev.off()
