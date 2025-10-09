source("utilities.R")
source("LoadData.R")
library(surrosurv)
library(survminer)
library(survival)
library(gridExtra)
library(grid)
library(WINS)
library(dplyr)


mydata=loading("gastadv")

# checking censoring rates
ndots <- 4
rprog_list1 <- matrix(NA, nrow = 1, ncol = ndots)
rprog_list2 <- matrix(NA, nrow = 1, ncol = ndots)
rdeath_list1 <- matrix(NA, nrow = 1, ncol = ndots)
rdeath_list2 <- matrix(NA, nrow = 1, ncol = ndots)

data1 <- mydata[mydata$A==1,]
data2 <- mydata[mydata$A==0,]

rprog_list1[1,ndots] <- 1-mean(data1$prog)
rdeath_list1[1,ndots] <- 1-mean(data1$death)
rprog_list2[1,ndots] <- 1-mean(data2$prog)
rdeath_list2[1,ndots] <- 1-mean(data2$death)

tau <- max(mydata$t.odeath)
tau_list <- seq(tau/ndots, tau, by = tau/ndots)
for (j in 1:(ndots-1)){
  sample_tru1 <- trunc_data(data1, arb.t=tau_list[j])
  rprog_list1[1,j] <- 1-mean(sample_tru1$prog)
  rdeath_list1[1,j] <- 1-mean(sample_tru1$death)
  sample_tru2 <- trunc_data(data2, arb.t=tau_list[j])
  rprog_list2[1,j] <- 1-mean(sample_tru2$prog)
  rdeath_list2[1,j] <- 1-mean(sample_tru2$death)
}

# plot the data
png(filename = "plots/realdata_12th.png", width = 2000, height = 1000, res = 300)

fit1 <- survfit(Surv(t.oprog, prog) ~ A, data = mydata)
plot1 <- ggsurvplot(fit1,data = mydata,
           xlab = "Time(days)", ylab = "Survival probability for PFS", 
           title = "", legend.title = "", 
           legend.labs = c("Control Group", "Treatment Group"),
           palette = c("#E7B800", "#2E9FDF"), 
           conf.int = TRUE,                   
           lwd = 0.8,                         
           pval = TRUE,
           pval.coord = c(700, 0.75),
           censor.shape = 1,                  
           censor.size = 0.1,
           ggtheme = theme_minimal(),
           legend = "bottom"
)

fit2 <- survfit(Surv(t.odeath, death) ~ A, data = mydata)
plot2 <- ggsurvplot(fit2,data = mydata,
           xlab = "Time(days)", ylab = "Survival probability for OS", 
           title = "", legend.title = "", 
           legend.labs = c("Control Group", "Treatment Group"), 
           palette = c("#E7B800", "#2E9FDF"),  
           conf.int = TRUE,                    
           lwd = 0.8,                         
           pval = TRUE,
           pval.coord = c(700, 0.75),
           censor.shape = 1,                   
           censor.size = 0.1,
           ggtheme = theme_minimal(),
           legend = "bottom"
)

grid.arrange(
  arrangeGrob(plot1$plot, plot2$plot, ncol = 2)
)

dev.off()

# log rank
logrank.res <- logrank.test(mydata)

# win statistic
data_wins <- data.frame(id=c(1:dim(mydata)[1]), arm=mydata$A, Delta_1=mydata$prog,
                       Delta_2=mydata$death, Y_1=mydata$t.oprog, Y_2=mydata$t.odeath)
res_tte <- capture.output(win.stat(data = data_wins, ep_type = "tte", arm.name = c("1","0"),
                                   priority = c(2:1), alpha = 0.05, digit = 5,
                                   stratum.weight = "unstratified", method = "unadjusted",
                                   pvalue = "two-sided"))
wins.res <- extract_values(res_tte)

# pseudo-death #####check dif tau
tau <- max(mydata$t.odeath)
ndots <- 4
tau_list <- seq(tau/ndots, tau, by = tau/ndots)
pseu.res.wel <- pseudeath.test.T(mydata[mydata$A==1,], mydata[mydata$A==0,], tau)

pseu.res.T2 <- pseudeath.test.T2(mydata[mydata$A==1,], mydata[mydata$A==0,], tau, ndots=4)

pseu.res.Z <- pseudeath.test.Z(mydata[mydata$A==1,], mydata[mydata$A==0,], tau, B=200)

set.seed(14)
for (j in 1:(ndots)){
  print(tau_list[j])
  sample_tru <- trunc_data(mydata,tau_list[j])
  #data_wins <- data.frame(id=c(1:dim(sample_tru)[1]), arm=sample_tru$A, Delta_1=sample_tru$prog,
  #                        Delta_2=sample_tru$death, Y_1=sample_tru$t.oprog, Y_2=sample_tru$t.odeath)
  #res_tte <- capture.output(win.stat(data = data_wins, ep_type = "tte", arm.name = c("1","0"),
  #                                   priority = c(2:1), alpha = 0.05, digit = 5,
  #                                   stratum.weight = "unstratified", method = "unadjusted",
  #                                   pvalue = "two-sided"))
  #wins.res <- extract_values(res_tte)
  #print(wins.res$win_ratio_p)
  pseu.res.Z <- pseudeath.test.Z(sample_tru[sample_tru$A==1,], sample_tru[sample_tru$A==0,], tau, B=500)
  print(pseu.res.Z$pseudo_death_B_p)
}

set.seed(1)
for (j in 1:(ndots)){
  print(tau_list[j])
  sample_tru <- trunc_data(mydata,tau_list[j])
  pseu.res.per <- pseudeath.test.per(sample_tru[sample_tru$A==1,], sample_tru[sample_tru$A==0,], tau, B=500)
  print(pseu.res.per$pseudo_death_per_p)
}

for (i in c(1:3)){
  pseu.res.wald <- pseudeath.test.wald(mydata[mydata$A==1,], mydata[mydata$A==0,], tau, ndots=4, B=5000)
  print(pseu.res.wald$pseudo_death_wal_n4_p)
}




# multi-state modelling
library(mstate)

tmat <- transMat(list(c(2, 3), c(3), c()), names = c("Start", "Progression", "Death"))

msdata <- msprep(
  data = mydata,
  trans = tmat,
  time = c(NA, "t.oprog", "t.odeath"),
  status = c(NA, "prog", "death"),
  keep = c("A"),
  start = list(state = 1, time = 0)
)

#sum(mydata$prog == 1 & mydata$death == 1 & mydata$t.oprog == mydata$t.odeath)

head(msdata)
msdata_expanded <- expand.covs(msdata, covs = "A", longnames = FALSE)

# Start → Progression（trans == 1）
cox1 <- coxph(Surv(Tstart, Tstop, status) ~ A, data = subset(msdata_expanded, trans == 1))
summary(cox1)

# Start → Death（trans == 2）
cox2 <- coxph(Surv(Tstart, Tstop, status) ~ A, data = subset(msdata_expanded, trans == 2))
summary(cox2)

# Progression → Death（trans == 3）
cox3 <- coxph(Surv(Tstart, Tstop, status) ~ A, data = subset(msdata_expanded, trans == 3))
summary(cox3)
