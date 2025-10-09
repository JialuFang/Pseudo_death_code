library(sos)
#findFn("colondeath")

library(frailtypack)
library(surrosurv)

loading<-function(dataname){
  if(dataname=="gastadv"){
    data("gastadv", package = "surrosurv")
    mydata=gastadv
    mydata=mydata[mydata$trialref==15,] #1  16 19 20 2  3  4  5  6  7  8  9  10 11 12 13 14 15 17 18
    myname=c("prog","t.oprog","death","t.odeath","A")
    mydata=mydata[,c(3,4,2,1,6)]
    colnames(mydata)=myname
    mydata$A=as.factor(mydata$A+0.5)
  }
  
  if(dataname=="dataOvarian"){
    data(dataOvarian, package = "frailtypack")
    mydata=dataOvarian
    mydata=mydata[mydata$trialID==1,]
    myname=c("prog","t.oprog","death","t.odeath","A")
    mydata=mydata[,c(5,4,7,6,3)]
    colnames(mydata)=myname
    mydata$A=as.factor(mydata$A)
  }
  
  if(dataname=="colon"){
    data(cancer, package = "survival")
    recurrence_data <- colon %>%
      filter(etype == 1) %>%
      select(id, rx, status, time) %>%
      rename(prog = status, t.oprog = time)
    
    death_data <- colon %>%
      filter(etype == 2) %>%
      select(id, status, time) %>%
      rename(death = status, t.odeath = time)
    
    mydata <- recurrence_data %>%
      inner_join(death_data, by = "id") %>%
      mutate(A = ifelse(rx == "Obs", 0, 1)) %>%
      select(id, prog, t.oprog, death, t.odeath, A)
    mydata$A=as.factor(mydata$A)
  }
  
  return(mydata)
}

#mean(mydata[mydata$A==1,]$death == 1)
#mean(mydata[mydata$A==0,]$death == 1)
#mean(mydata[mydata$A==1,]$death == 0 & mydata[mydata$A==1,]$prog == 1)
#mean(mydata[mydata$A==0,]$death == 0 & mydata[mydata$A==0,]$prog == 1)
#mean(mydata[mydata$A==1,]$death == 0 & mydata[mydata$A==1,]$prog == 0)
#mean(mydata[mydata$A==0,]$death == 0 & mydata[mydata$A==0,]$prog == 0)

