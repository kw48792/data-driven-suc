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

#######################################
## Get predictions for future wind power by Persistence Model
## data: historical wind power data
## step: step ahead
## h: # of historical obs
## Return: Empirical cdf of the predictive distribution
#######################################
getPosterior_Persistence = function(data,step,h){
  T = length(data)
  for(i in 1:step){
    samples = data[T]-(data[T-seq(0,h-i,1)]-data[T-seq(i,h,1)])
    samp = sample(samples,50,replace = T)
    write.table(t(samp),"pred3.csv",sep=",",append=TRUE,row.names=FALSE,quote = FALSE,col.names=FALSE)
  }
  return()
}


TT = 100
start = 1
for(start in seq(from=101, to=864, by=4)){
  print(start)
  step = 4
  data = WindPowerData$Real[(start-TT):(start+TT+step-1)]
  #data = log(data/(1-data))
  train = data[1:TT]
  test_point = data[TT+1]
  h = TT-step+1
  res = getPosterior_Persistence(train,step,h)
 # pi_width_99 = matrix(0,1,step)
 # pi_cover_99 = matrix(0,1,step)
 # pi_width_90 = matrix(0,1,step)
 # pi_cover_90 = matrix(0,1,step)
  #pi_width_50 = matrix(0,1,step)
  #pi_cover_50 = matrix(0,1,step)
 # error = matrix(0,1,step)
 # skillScores = matrix(0,1,step)
 # alphas = seq(0.025,0.975,0.025)
  for(i in 1:step){
  #  samp = sample(res[[i]],50,replace = T)
    
  #  pred = data[TT]
  #  test_point = data[TT+i]
  # pi_width_99[i] = quantile(cdf,0.005)-quantile(cdf,0.995)
  #  pi_cover_99[i] = (quantile(cdf,0.995)>=test_point & test_point >= quantile(cdf,0.005))*1
  #  pi_width_90[i] = quantile(cdf,0.95)-quantile(cdf,0.05)
  #  pi_cover_90[i] = (quantile(cdf,0.95)>=test_point & test_point >= quantile(cdf,0.05))*1
    #pi_width_50[i] = quantile(cdf,0.75)-quantile(cdf,0.25)
    #pi_cover_50[i] = (quantile(cdf,0.75)>=test_point & test_point >= quantile(cdf,0.25))*1
   # error[i] =  pred-test_point
  #  skillScores[i] = sum(((test_point<=quantile(cdf,alphas))*1-alphas)*(test_point-quantile(cdf,alphas)))
  }
  #write.table(skillScores,"M100_Peris_skill.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
  #write.table(pi_width_99,"M100_Peris_piWid_99.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
  #write.table(pi_cover_99,"M100_Peris_Cover_99.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
  #write.table(pi_width_90,"M100_Peris_piWid_90.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
  #write.table(pi_cover_90,"M100_Peris_Cover_90.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
  #write.table(pi_width_50,"M100_IMSAR_piWid_50.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
  #write.table(pi_cover_50,"M100_IMSAR_Cover_50.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)
  #write.table(error,"M100_Peris_Error.txt",sep=",",append=TRUE, eol="\n",row.names=FALSE,col.names=FALSE)

  }

