
require(diptest)
require(optimbase)

wd<- 'Enter your working directory here'
#----------------------------------------------
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
d<- as.matrix(read.csv('lung2_x.txt',header=F,sep='\t'))
true.class<- as.matrix(c(rep(0,139),rep(1,64)))
nclass<- length(unique(true.class))
d<- as.matrix(transpose(d))
d<- d[1:203,]
d<- scale(d,TRUE,TRUE)
p<-ncol(d)
n<- nrow(d)
tdata<- d
out <- excess_mass(tdata)
print(sum(out))

idx<- (out==1)
D<- scale(tdata[,idx],TRUE,TRUE)

CER.VKM<- matrix(0,30,1)
for (i in 1:30){
  set.seed(i^2)
  vkm<- kmeans(D,nclass, iter.max = 50, nstart = 30,algorithm = "Hartigan-Wong")
  predclass.full.km<- vkm$cluster
  
  Q<- as.matrix(true.class)
  m<-n
  num<-0
  den<- m*(m-1)/2
  for (ii in 1:(m-1)){
    
    tt<- predclass.full.km[ii]-predclass.full.km[(ii+1):m]
    tt[tt!=0]=1
    tt <- 1-tt
    dd<- Q[ii]-Q[(ii+1):m]
    dd[dd!=0]=1
    dd <- 1-dd
    num<- sum(abs(tt-dd))+num
    
  }
  
  CER.VKM[i]<- num/den
  print(i)
}
avgCER.KM<- c(mean(CER.VKM),sd(CER.VKM)/sqrt(30))

