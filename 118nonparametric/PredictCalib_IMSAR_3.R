setwd("F:/DataSUC/code11.15/Wind power prediction codes (2)/Wind power prediction codes 013020")
source("IMS_AR_Lib.R")
#source("~/Desktop/Time Dependence/IMS_AR_Lib.R")
library(readr)
library('NHMSAR')
library('msm')
library('pscl')
#WindPowerData <- read.table("~/Desktop/Time Dependence/Wind_Empirical/WindPowerData.txt")
WindPowerData <- read.table("WindPowerData3.txt")
WindPowerData = WindPowerData/4
colnames(WindPowerData)<-c("Real")
set.seed(34)
#######################################
## Get predictions for future wind power by IMSAR
## data: historical wind power data
## step: step ahead
## p: lag of order
## Return: Empirical cdf of the predictive distribution
#######################################
getPosterior_IMSAR = function(data,step,p){
  res = IMS_AR_P(data,100,100,1,1,0,p)
  IMSAR_by_model = matrix(0,100,step)
  IMSAR_all = c()
  #IMSAR_predictions = c()
  for(i in 1:100){
    sample_model = res[[i]]
    count = sample_model[[1]]
    Psi = sample_model[[2]]
    sigma = sample_model[[3]]
    S = sample_model[[4]]
    N = sample_model[[5]]
    K = dim(N)[1]
    count[1:p] = count[1:p]-1
    S = S[(p+1):length(S)]
    removed = c()
    remove = 0
    K_removed = c()
    for(k in 1:K){
      if(count[k]==0){
        K = K-1
        K_removed = c(K_removed,k)
        remove = remove+1
      }
      removed = c(removed,remove)
    }
    if(length(K_removed) > 0){
      N = N[-K_removed,-K_removed]
      if(!is.matrix(N)){
        N = as.matrix(N)
      }
      Psi = matrix(Psi[-K_removed, ],K,(p+1))
      sigma = sigma[-K_removed]
      count = count[-K_removed]
    }
    for(tt in 1:length(S)){
      S[tt] = S[tt]-removed[S[tt]]
    }
    P = N
    for(k in 1:K){
      P[k,] = N[k,]/sum(N[k,])
    }
    
    prev = data[length(data):(length(data)-p+1)]
    prev_S = S[length(S)]
    IMSAR_pred = matrix(0,100,step)
    for(rep in 1:100){
      IMSAR_pred[rep,] = getPrediction_IMSAR(prev,p,prev_S,P,Psi,sigma,step)
      IMSAR_all = rbind(IMSAR_all,IMSAR_pred[rep,])
    }
    IMSAR_by_model[i,] = apply(IMSAR_pred,2,mean)
    #IMSAR_predictions= rbind(IMSAR_predictions,data-resid)
  }
#  IMSAR_all = exp(IMSAR_all)/(1+exp(IMSAR_all))
  res = list()
  for (i in 1:step){
    res=c(res,ecdf(IMSAR_all[,i]))
    samp = sample(IMSAR_all[,i],50,replace = T)
    write.table(t(samp),"pred3IMSARL3.csv",sep=",",append=TRUE,row.names=FALSE,quote = FALSE,col.names=FALSE)
    
  }
  return(c(res,apply(IMSAR_all,2,mean)))
}

TT = 100
start = 1
for(start in seq(from=101, to=864, by=4)){
  print(start)
  step = 4
  data = WindPowerData$Real[(start-TT):(start+TT+step-1)]
  train = data[1:TT]
  ramp = 0
  if(isRamp(train)){
    ramp = 1
  }
  p = 1
 # curTime = proc.time()[1]
  res = getPosterior_IMSAR(train,step,p)
 # cputime = proc.time()[1]-curTime
 # pits = matrix(0,1,step)
 # pi_width_99 = matrix(0,1,step)
  #pi_cover_99 = matrix(0,1,step)
  #pi_width_90 = matrix(0,1,step)
 # pi_cover_90 = matrix(0,1,step)
 # error = matrix(0,1,step)
 # skillScores = matrix(0,1,step)
 # alphas = seq(0.025,0.975,0.025)
  for(i in 1:step){
 #   cdf = res[[i]]
 #   pred = res[[step+i]]
  #  test_point = data[TT+i]
 #   pi_width_99[i] = quantile(cdf,0.005)-quantile(cdf,0.995)
  #  pi_cover_99[i] = (quantile(cdf,0.995)>=test_point & test_point >= quantile(cdf,0.005))*1
  #  pi_width_90[i] = quantile(cdf,0.95)-quantile(cdf,0.05)
  #  pi_cover_90[i] = (quantile(cdf,0.95)>=test_point & test_point >= quantile(cdf,0.05))*1
    #pi_width_50[i] = quantile(cdf,0.75)-quantile(cdf,0.25)
    #pi_cover_50[i] = (quantile(cdf,0.75)>=test_point & test_point >= quantile(cdf,0.25))*1
 #   error[i] =  pred-test_point
 #   skillScores[i] = sum(((test_point<=quantile(cdf,alphas))*1-alphas)*(test_point-quantile(cdf,alphas)))
  }
 # skillScores = round(skillScores,2)
 # pi_width_99 = round(pi_width_99,3)
 # pi_width_90 = round(pi_width_90,3)
  #pi_width_50 = round(pi_width_50,3)
 # error = round(error,3)
#  write.table(skillScores,"M100_IMSAR_skill.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
 # write.table(pi_width_99,"M100_IMSAR_piWid_99.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
 # write.table(pi_cover_99,"M100_IMSAR_Cover_99.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
 # write.table(pi_width_90,"M100_IMSAR_piWid_90.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
 # write.table(pi_cover_90,"M100_IMSAR_Cover_90.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
  #write.table(pi_width_50,"M100_IMSAR_piWid_50.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
  #write.table(pi_cover_50,"M100_IMSAR_Cover_50.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
  #write.table(error,"M10_IMSAR_Error.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
#  write.table(cputime,"M100_IMSAR_Time.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
#  write.table(ramp,"M100_IMSAR_Ramp.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
}

