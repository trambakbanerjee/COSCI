
require(diptest)

wd<-'Enter your working directory here'
#----------------------------------------------------
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
load(paste(wd,"/su.RData",sep=""))

d<- matrix(unlist(su$x),ncol=5565,byrow=TRUE)
true.class<-as.matrix(as.numeric(su$y))
nclass<- length(unique(true.class))
d<- scale(d,TRUE,TRUE)
p<-ncol(d)
n<- nrow(d)
out <- excess_mass(d)

print(sum(out))



