
require(diptest)
require(optimbase)

wd<- 'Enter your working directory here'
#------------------------------------------------------------------------------
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
load(paste(wd,"/cardiodata.RData",sep=""))
d<- cardiodata
true.class<- as.matrix(c(rep(1,44),rep(2,19)))
nclass<- length(unique(true.class))
d<-as.matrix(transpose(d))
d<- scale(d,TRUE,TRUE)
p<-ncol(d)
n<- nrow(d)
out <- excess_mass(d)

print(sum(out))



