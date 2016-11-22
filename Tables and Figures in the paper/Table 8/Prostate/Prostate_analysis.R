

require(ggplot2)
require(optimbase)
require(bbmle)
require(sparcl)
require(gmodels)
require(gridExtra)

wd<- 'Enter your working directory here'
#----------------------------------------------------------------
cosci<-function(d,min.alpha,small.perturbation){
  
  p <- ncol(d)
  dim.score<- matrix(0,nrow=p,ncol=2)
  for (j in 1:p){
    #remove ties
    x<- d[,j]
    x<-x+rnorm(length(x))*small.perturbation
    x<-sort(x);
    n=length(x);
    #initializiations 
    C=n;#representing number of clusters by C
    ind<-1:C;#cluster index of observations
    ind.along.path<-c();  
    l.old<-0;
    lambda<-0;
    merge.points<-c();
    no.of.clusters<-c();
    probability<-c();
    boundaries<-c();
    freq<-rep(1,n);
    centroids<-x;
    unique.ind<-ind;
    counter=3;
    state<- 1
    while(state==1){
      l.new<-(centroids[2:C]-centroids[1:(C-1)])/(freq[2:C]+freq[1:(C-1)]);
      l.new.mod<-l.new[l.new>=l.old];
      l.old<-min(l.new.mod);
      i=match(l.old,l.new);
      merge.points<-c(merge.points,C);
      a.2<-min(x[(ind==unique.ind[i+1])]);
      a.1<-max(x[(ind==unique.ind[i])])
      probability<-rbind(probability,c(freq[i],freq[i+1])/n);
      centroids[i]<-(centroids[i]*freq[i]+centroids[i+1]*freq[i+1])/(freq[i]+freq[i+1]);
      centroids<-centroids[-(i+1)];
      freq[i]=freq[i]+freq[i+1];
      freq<-freq[-(i+1)];
      ind[(ind==unique.ind[i+1])]<-unique.ind[i];
      unique.ind<-unique.ind[-(i+1)];
      lambda<-c(lambda,l.old);
      C=C-1;
      if (C==1){
        #paste('We only have one cluster if the penalization parameter exceeds');
        l<-list(prob=probability);
        state<-0;
      }
    }
    prob1<- as.matrix(l$prob)
    score<- as.matrix(cbind(apply(prob1,1,min),apply(prob1,1,sum)))
    grp<- matrix(0,nrow=nrow(prob1),ncol=1)
    idx <- c(score[,1]>min.alpha & score[,2]>=0.5)
    grp[idx]<-score[idx,1]
    m<- max(grp)
    l<- min(grp)
    #prob1<- cbind(prob1,grp)
    dim.score[j,1:2]<- cbind(max(grp[idx]),mean(grp[idx]))
    
  }
  dim.score[is.nan(dim.score)]<-0
  
  return(dim.score)
}

setwd(wd)
d<- as.matrix(read.csv('prostate_x.txt',header=F,sep=' '))
true.class<- as.matrix(read.csv('prostate_y.txt',header=F,sep=' '))
d<- as.matrix(transpose(d))
d<- scale(d,TRUE,TRUE)
p<-ncol(d)
n<- nrow(d)
tdata<- d
score<-as.matrix(cosci(tdata,0.001,10^{-6}))
score.raw<- score[,1]
x<-score.raw/0.5

#-------------------------------------------------------------------
require(bbmle)
ux <- unique(x)
mode<- ux[which.max(tabulate(match(x, ux)))]
see<- order(x)
interval<- c(0,x[see[ceiling(0.9*p)]])
xi<- x[x<=interval[2]]
N0<- length(xi)
N<- length(x)
theta<- (N0/N)

nll<-function(a,b,n0,u){
  
  -sum(log(dbeta(xi,a,b)))+n0*log(pbeta(u,a,b))
}

out <- mle2(nll,start=list(a=2,b=5),data=list(n0=N0,u=interval[2]))

ahat<-as.vector(out@fullcoef[1])
bhat<- as.vector(out@fullcoef[2])
pi0<- theta/pbeta(interval[2],ahat,bhat)

f.null<-matrix(0,p,1)
index<- x<=interval[2]
f.null<- dbeta(x,ahat,bhat)

histinfo<-hist(x,breaks=100)
xx<- histinfo$mids
X<- cbind(xx,xx^{2},xx^{3},xx^{4},xx^{5},xx^{6})
y<- histinfo$counts
m<-glm(y~X,family="poisson")
beta<-as.vector(m$coefficients)
beta<- matrix(beta,ncol=1,byrow=TRUE)
X<-cbind(matrix(1,p,1),x,x^{2},x^{3},x^{4},x^{5},x^{6})
v<- exp(X%*%beta)
fhat<- v/(length(x)*0.01)

d1<-as.data.frame(cbind(x,(pi0*f.null)))
colnames(d1)<-c("X","f0")
d2<-as.data.frame(cbind(x,fhat))
colnames(d2)<-c("X","fhat")

ggplot() + geom_line(data=d1,(aes(x=X,y=f0)),color='blue')+geom_line(data=d2,(aes(x=X,y=fhat)),color='red')

fdr<- (pi0*f.null)/fhat
# ---------------- We will do MDR screening now -----------------------------------------

Ti<- fdr[order(fdr)]
Ti[Ti>1]<- 1
Yi<- 1-Ti
eps<- 1-pi0
nom.level<-(1/log(p))
cutoff.1<- p*eps*nom.level
ks<- min(which((rev(cumsum(rev(Yi))))<=cutoff.1))
idx.1<- (fdr<=Ti[ks])
alpha.1<- score.raw[idx.1]
select.d1<- d[,idx.1]
r.s<- sum(idx.1)

# ---------------- Now do FPR signal discovery -----------------------------------------
fdr.1<- fdr[idx.1]
Ti<- fdr.1[order(fdr.1)]
Ti[Ti>1]<- 1
cutoff.2<- min(nom.level,0.1)
kd<-  max(which((cumsum(Ti)/(1:r.s))<=cutoff.2))
idx.2<- (fdr.1<=Ti[kd])
alpha.2<- alpha.1[idx.2]
select.d2<- select.d1[,idx.2]
write.csv(select.d2,paste(wd,"/selected.csv",sep=""))

#------------------------- COSCI+SpKM and COSCI+KM --------------------------

D<- scale(select.d2,TRUE,TRUE)

km.perm <- KMeansSparseCluster.permute(D,K=2,nperms=25)
km.out <- KMeansSparseCluster(D,K=2,wbounds=km.perm$bestw)
predclass.spkm<- km.out[[1]]$Cs

Q<- as.matrix(true.class)
m<-n
num<-0
den<- m*(m-1)/2
for (ii in 1:(m-1)){
  
  tt<- predclass.spkm[ii]-predclass.spkm[(ii+1):m]
  tt[tt!=0]=1
  tt <- 1-tt
  dd<- Q[ii]-Q[(ii+1):m]
  dd[dd!=0]=1
  dd <- 1-dd
  num<- sum(abs(tt-dd))+num
  
}
CER.COSCISpKM<- num/den

CER.COSCIVKM<- matrix(0,30,1)
for (i in 1:30){
  set.seed(i^2)
  vkm<- kmeans(D,2, iter.max = 50, nstart = 30,algorithm = "Hartigan-Wong")
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
  
  CER.COSCIVKM[i]<- num/den
  print(i)
}
avgCER.COSCIKM<- c(mean(CER.COSCIVKM),sd(CER.COSCIVKM)/sqrt(30))

#------------------------- Full SpKM and KM --------------------------

D<- scale(d,TRUE,TRUE)
nclass<-2

km.perm <- KMeansSparseCluster.permute(D,K=nclass,nperms=25)
km.out <- KMeansSparseCluster(D,K=nclass,wbounds=km.perm$bestw)
predclass.full.spkm<- km.out[[1]]$Cs

Q<- as.matrix(true.class)
m<-n
num<-0
den<- m*(m-1)/2
for (ii in 1:(m-1)){
  
  tt<- predclass.full.spkm[ii]-predclass.full.spkm[(ii+1):m]
  tt[tt!=0]=1
  tt <- 1-tt
  dd<- Q[ii]-Q[(ii+1):m]
  dd[dd!=0]=1
  dd <- 1-dd
  num<- sum(abs(tt-dd))+num
  
}
CER.SpKM<- num/den


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
