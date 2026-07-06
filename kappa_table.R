if (!require(ggplot2)) {
  install.packages("ggplot2")
}
library(ggplot2)
if (!require(gridExtra)) {
  install.packages("gridExtra")
}
library(gridExtra)


matrix_kappa2 <- matrix(0, nrow=4, ncol=3)
matrix_kappa3 <- matrix(0, nrow=4, ncol=3)

col_names <- c("estimated kappa mean", "estimated kappa sd", "MSE")
row_names <- c("n=100", "n=200", "n=500", "n=1000")


dimnames(matrix_kappa2) <- list(row_names, col_names)
dimnames(matrix_kappa3) <- list(row_names, col_names)

theta_CR <- 2

true_kappa2 <- 0.252; true_kappa3 <- 0.277
true_kappa2 <- 0.247; true_kappa3 <- 0.239
true_kappa2 <- 0.303; true_kappa3 <- 0.471
true_kappa2 <- 0.196; true_kappa3 <- 0.129


true_kappa2 <- 0.151; true_kappa3 <- 0.253
true_kappa2 <- 0.157; true_kappa3 <- 0.256
true_kappa2 <- 0.174; true_kappa3 <- 0.276
true_kappa2 <- 0.187; true_kappa3 <- 0.292
true_kappa2 <- 0.202; true_kappa3 <- 0.309


samplesize <- c(100, 200, 500, 1000)
for (i in c(1:4)){
  #file.name <- paste0("./kappa_accuracy_data/","kappa", "_size",  samplesize[i], "_simu4",".Rdata")
  file.name <- paste0("./kappa_accuracy_cor_data/","cor_kappa", "_size",  samplesize[i], 
                        "_simu_",theta_CR*10,".Rdata")
  print(file.name)
  load(file.name) 
  matrix_kappa2[i,1] <- mean(results$kappa2)
  matrix_kappa2[i,2] <- sd(results$kappa2)
  matrix_kappa2[i,3] <- mean((true_kappa2 - results$kappa2)^2)
  
  
  matrix_kappa3[i,1] <- mean(results$kappa3)
  matrix_kappa3[i,2] <- sd(results$kappa3)
  matrix_kappa3[i,3] <- mean((true_kappa3 - results$kappa3)^2)
}

matrix_kappa2
matrix_kappa3

