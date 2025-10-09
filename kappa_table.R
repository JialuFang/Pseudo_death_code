if (!require(ggplot2)) {
  install.packages("ggplot2")
}
library(ggplot2)
if (!require(gridExtra)) {
  install.packages("gridExtra")
}
library(gridExtra)


matrix_kappa2 <- matrix(0, nrow=4, ncol=5)
matrix_kappa3 <- matrix(0, nrow=4, ncol=5)

col_names <- c("estimated kappa mean", "estimated kappa sd", 
               "true kappa mean", "true kappa sd", "MSE")
row_names <- c("n=100", "n=200", "n=500", "n=1000")


dimnames(matrix_kappa2) <- list(row_names, col_names)
dimnames(matrix_kappa3) <- list(row_names, col_names)

theta_CR <- 1

samplesize <- c(100, 200, 500, 1000)
for (i in c(1:4)){
  #file.name <- paste0("./kappa_accuracy_data/","kappa", "_size",  samplesize[i], "_simu1",".Rdata")
  file.name <- paste0("./kappa_accuracy_cor_data/","cor_kappa", "_size",  samplesize[i], 
                        "_simu_",theta_CR*10,".Rdata")
  print(file.name)
  load(file.name) 
  matrix_kappa2[i,1] <- mean(results$kappa2)
  matrix_kappa2[i,2] <- sd(results$kappa2)
  matrix_kappa2[i,3] <- mean(results$true_kappa2)
  matrix_kappa2[i,4] <- sd(results$true_kappa2)
  matrix_kappa2[i,5] <- mean((results$true_kappa2 - results$kappa2)^2)
  
  
  matrix_kappa3[i,1] <- mean(results$kappa3)
  matrix_kappa3[i,2] <- sd(results$kappa3)
  matrix_kappa3[i,3] <- mean(results$true_kappa3)
  matrix_kappa3[i,4] <- sd(results$true_kappa3)
  matrix_kappa3[i,5] <- mean((results$true_kappa3 - results$kappa3)^2)
}

matrix_kappa2
matrix_kappa3

