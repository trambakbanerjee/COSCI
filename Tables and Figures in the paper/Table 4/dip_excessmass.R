
wd<-'Enter your working directory here'

#-------------------------------------------------------
require(diptest)
excess_mass<-function(data){
  
  p<- ncol(data)
  pvals<- matrix(0,p,1)
  decision<- matrix(0,1,p)
  for (i in 1:p){
    pvals[i]<- dip.test(data[,i])$p.value
    
  }
  index<- order(pvals)
  t<- as.matrix(p.adjust(pvals[index],method = "BH"))
  pvals<- t[order(index)]
  decision[,which(pvals<=0.05)]<-1
  return(decision)
}

setwd(wd)
datalocation <- wd
source("coscilibrary.R")

reps<-50
p<-5000
s<-7
n<- 200
temp<- paste(datalocation,toString(n),sep="/")
out<- matrix(0,reps,p)
for (r in 1:reps){
  
  name<- paste('data_',toString(r),'.RData',sep="")
  load(paste(temp,name,sep="/"))
  d <- d[,1:p]
  out[r,] <- excess_mass(d)
}

FN.excessmass<- rowSums(out[,1:s]==0)# Signal identified as Noise
FP.excessmass<- rowSums(out[,(s+1):p]>0)# Noise identified as signal
tab_200<- cbind(sum(FN.excessmass)/reps,sd(FN.excessmass)/sqrt(reps),
                 sum(FP.excessmass)/reps,sd(FP.excessmass)/sqrt(reps))

n<- 1000
temp<- paste(datalocation,toString(n),sep="/")
out<- matrix(0,reps,p)
for (r in 1:reps){
  
  name<- paste('data_',toString(r),'.RData',sep="")
  load(paste(temp,name,sep="/"))
  d <- d[,1:p]
  out[r,] <- excess_mass(d)
}

FN.excessmass<- rowSums(out[,1:s]==0)# Signal identified as Noise
FP.excessmass<- rowSums(out[,(s+1):p]>0)# Noise identified as signal
tab_1000<- cbind(sum(FN.excessmass)/reps,sd(FN.excessmass)/sqrt(reps),
                 sum(FP.excessmass)/reps,sd(FP.excessmass)/sqrt(reps))

n<- 2500
temp<- paste(datalocation,toString(n),sep="/")
out<- matrix(0,reps,p)
for (r in 1:reps){
  
  name<- paste('data_',toString(r),'.RData',sep="")
  load(paste(temp,name,sep="/"))
  d <- d[,1:p]
  out[r,] <- excess_mass(d)
}

FN.excessmass<- rowSums(out[,1:s]==0)# Signal identified as Noise
FP.excessmass<- rowSums(out[,(s+1):p]>0)# Noise identified as signal
tab_2500<- cbind(sum(FN.excessmass)/reps,sd(FN.excessmass)/sqrt(reps),
                 sum(FP.excessmass)/reps,sd(FP.excessmass)/sqrt(reps))



