source("utilities.R")
S_prog_fct1 <- function(t){exp(-3/2*t)}
inv_S_prog_fct1 <- function(x){-2/3*log(x)}

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

followup_para1 <- 2; followup_para2 <- 2; tau <- 2; samplesize <- 200; ndots <- 4
nsample <- 500
death_list1 <- rep(0, nsample); death_list2 <- rep(0, nsample)
prog_list1 <- rep(0, nsample); prog_list2 <- rep(0, nsample)

for (i in 1:nsample){
  set.seed(i)
  res1 <- gen.data(n=samplesize, followup_dist='unif', followup_para=followup_para1, 
                                    inv_S_prog_fct1, S_prog_fct1, S_death_fct1, inv_S_death_fct1, 
                                    tau=tau, theta=1, A=1)
  res2 <- gen.data(n=samplesize, followup_dist='unif', followup_para=followup_para2, 
                                    inv_S_prog_fct2, S_prog_fct2, S_death_fct2, inv_S_death_fct2, 
                                    tau=tau, theta=1, A=0)
  death_list1[i]<-res1$rdeath; death_list2[i]<-res2$rdeath
  prog_list1[i]<-res1$rprog; prog_list2[i]<- res2$rprog
}

mean(death_list1); mean(death_list2)
mean(prog_list1); mean(prog_list2)
