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

plot_validation_figures <- function(validation.df, main_suffix = "") {
  z <- validation.df$z_boot
  z <- z[is.finite(z)]
  
  hist(
    z,
    breaks = 30,
    probability = TRUE,
    main = bquote("Histogram of " * Z[B] ~ .(main_suffix)),
    xlab = expression(Z[B] * " statistic"),
    col = "lightgray",
    border = "white",
    cex.main=0.9, cex.lab=0.8
  )
  curve(dnorm(x, mean = 0, sd = 1), add = TRUE, col = "red", lwd = 2)
  
  qqnorm(z, main = bquote("Q-Q plot of " * Z[B] ~ .(main_suffix)), 
         pch = 19, cex = 0.5, cex.main=0.9, cex.lab=0.8)
  qqline(z, col = "red", lwd = 2)
  
}

files <- paste0("./test_results/Validation_bootstrap_n500_simu", 1:4, ".Rdata")

pdf("Normality_check.pdf", width = 5, height = 7)

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(4, 2), mar = c(4, 4, 3, 1))

for (i in 1:4) {
  load(files[i])   # load validation.df
  plot_validation_figures(validation.df, main_suffix = paste("under Scenario", i))
}
par(oldpar)
dev.off()

