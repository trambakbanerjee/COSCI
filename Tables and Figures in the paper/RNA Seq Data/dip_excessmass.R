
require(diptest)
require(optimbase)

wd<- 'Enter your working directory here'
#-----------------------------------------------------------

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
#------------------------------------------------------------

setwd(wd)
load("RNASeq.RData")
genes<-colnames(d)
genes<-as.data.frame(genes[-c(1,2)])
true.class<- c(as.numeric(d[,2]))
nclass<- length(unique(true.class))
tdata<- d[,-c(1,2)]
tdata<- matrix(as.numeric(unlist(tdata)),nrow=nrow(tdata))
tdata[tdata>1]<- 1
tdata<-log(tdata)
n<-nrow(tdata)
p<-ncol(tdata)
out <- excess_mass(tdata)

print(sum(out))
idx<- (out==1)
d<- tdata[,idx]
D<- scale(d,TRUE,TRUE)
CER.EMKM<- matrix(0,30,1)

for (i in 1:30){
  set.seed(i^2)
  vkm<- kmeans(D,nclass, iter.max = 50, nstart = 30,algorithm = "Hartigan-Wong")
  predclass.km<- vkm$cluster
  
  Q<- as.matrix(true.class)
  m<-n
  num<-0
  den<- m*(m-1)/2
  for (ii in 1:(m-1)){
    
    tt<- predclass.km[ii]-predclass.km[(ii+1):m]
    tt[tt!=0]=1
    tt <- 1-tt
    dd<- Q[ii]-Q[(ii+1):m]
    dd[dd!=0]=1
    dd <- 1-dd
    num<- sum(abs(tt-dd))+num
    
  }
  
  CER.EMVKM[i]<- num/den
  print(i)
}
avgCER.BMTKM<- c(mean(CER.BMTVKM),sd(CER.BMTVKM)/sqrt(30))

