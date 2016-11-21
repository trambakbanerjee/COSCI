
require(optimbase)
require(MASS)

#------------------------------------------------------
# This script uses snowfall package
#------------------------------------------------------

wd<- 'Enter your working directory here'
ncpu<- 8 #Enter the number of workers here
#----------------------------------------------------
setwd(wd)
cosci<-function(dat,min.alpha,small.perturbation){
  
  p <- ncol(dat)
  dim.score<- matrix(0,nrow=p,ncol=3)
  for (j in 1:p){
    #remove ties
    x<- dat[,j]
    set.seed(1)
    nn<- rnorm(length(x))*small.perturbation
    x<-x+nn
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
    dim.score[j,]<- cbind(max(grp[idx]),mean(grp[idx]),sum(idx))
    print(j)
    
  }
  dim.score[is.nan(dim.score)]<-0
  
  return(dim.score[,1])
}
cosci_2wayselect<-function(score,gamma){
  
  require(ggplot2)
  require(bbmle)
  require(gridExtra)
  
  p<- length(score)
  x<- score/0.5
  see<- order(x)
  interval<- c(0,x[see[ceiling(gamma*p)]])
  xi<- x[x<=interval[2]]
  N0<- length(xi)
  N<- length(x)
  theta<- (N0/N)
  
  nll<-function(a,b,n0,u){
    
    -sum(log(dbeta(xi,a,b)))+n0*log(pbeta(u,a,b))
  }
  
  out <- mle2(nll,start=list(a=0.2,b=5),data=list(n0=N0,u=interval[2]))
  
  ahat<-as.vector(out@fullcoef[1])
  bhat<- as.vector(out@fullcoef[2])
  pi0<- min(theta/pbeta(interval[2],ahat,bhat),0.99)
  
  f.null<-matrix(0,p,1)
  index<- x<=interval[2]
  f.null<- dbeta(x,ahat,bhat)
  
  histinfo<-hist(x,breaks=150,plot=FALSE)
  xx<- histinfo$mids
  X<- cbind(xx,xx^{2},xx^{3})
  y<- histinfo$counts
  m<-glm(y~X,family="poisson")
  beta<-as.vector(m$coefficients)
  beta<- matrix(beta,ncol=1,byrow=TRUE)
  X<-cbind(matrix(1,p,1),x,x^{2},x^{3})
  v<- exp(X%*%beta)
  gap<-histinfo$breaks[2]-histinfo$breaks[1]
  fhat<- v/(length(x)*gap)
  
  fdr<- (pi0*f.null)/fhat
  fdr[fdr>1]=1
  
  d1<-as.data.frame(cbind(x,(pi0*f.null)))
  colnames(d1)<-c("X","f0")
  d2<-as.data.frame(cbind(x,fhat))
  colnames(d2)<-c("X","fhat")
  ggplot() + geom_line(data=d1,(aes(x=X,y=f0)),color='blue')+geom_line(data=d2,(aes(x=X,y=fhat)),color='red')
  
  # ---------------- We will do MDR screening now -----------------------------------------
  
  Ti<- fdr[order(fdr)]
  Ti[Ti>1]<- 1
  Yi<- 1-Ti
  eps<- 1-pi0
  nom.level<-(1/log(p))
  cutoff.1<- p*eps*nom.level
  ks<- min(which((rev(cumsum(rev(Yi))))<=cutoff.1))
  idx.1<- (fdr<=Ti[ks])
  alpha.1<- score[idx.1]
  select.d1<- (1:p)*idx.1
  select.d1<- select.d1[select.d1>0]
  r.s<- sum(idx.1)
  
  # ---------------- Now do FPR signal discovery -----------------------------------------
  fdr.1<- fdr[idx.1]
  Ti<- fdr.1[order(fdr.1)]
  Ti[Ti>1]<- 1
  cutoff.2<- min(nom.level,0.1)
  kd<-  max(which((cumsum(Ti)/(1:r.s))<=cutoff.2))
  idx.2<- (fdr.1<=Ti[kd])
  alpha.2<- alpha.1[idx.2]
  select.d2<- select.d1[idx.2]
  
  return(list("selected"=select.d2))
  
}

#--------------------------------------------------------------------------
# Example 1 - Non Product structure (Totally 4 dims)

reps <- 10
p<- 25
idx<- combn(p,2,simplify=TRUE)
noise<- p-4   
n1 <- 1000
n2 <- 1000
n<-n1+n2
mu1 <- c(0.9, -0.9)
mu2 <- c(-0.9, 0.9)
mu3<- c(5,5)
rho1 <- 0.9
rho2<- -0.9
sigma1 = matrix(c(1, rho1,rho1,1),2)
sigma2 = matrix(c(1, rho2,rho2,1),2)

# Generate points on the unit circle
set.seed(30^3)
u1 <- transpose(runif(10,min=-1,max=1))
u2 <- sqrt(1-u1^2)
R = rbind(cbind(u1,u2),cbind(u1,-u2),cbind(1,0),cbind(-1,0),cbind(0,1),cbind(0,-1))
out<-matrix(0,24,reps)
outdir<-matrix(0,reps,1)
S<- matrix(0,reps,p+ncol(idx))
dir.opt<-matrix(0,ncol(idx),2*reps)
dir.threshold<- 0.9

library(snowfall)  
sfInit(parallel = TRUE,cpu = ncpu)
sfLibrary(optimbase)

for (r in 1:reps){
  #0. First generate a suitable data
  set.seed(10*r)
  x1 <- mvrnorm(n1, mu1, sigma1)
  set.seed(10*r)
  y1 <- mvrnorm(n2, mu2, sigma1)
  X1<- rbind(x1,y1)
  set.seed(10*r)
  X2 <- c(rbeta(n/2, 4, 6, ncp = 0),rbeta(n/2, 7,3, ncp = 0))
  set.seed(10*r)
  X3<- c(rlnorm(n/2,0.2,0.35),rnorm(n/2,4,0.5))
  set.seed(10*r)
  y = matrix(rnorm(n*noise,0,1),n,noise)
  D <- cbind(X1[sample(nrow(X1)),],X2,X3,y)
  P<-scale(D,center=TRUE,scale=TRUE)
  
  outcosci.1<-cosci(P,0,10^{-6})
   
  wrapper1<- function(dir){
    s<- transpose(R[dir,])
    XX<- matrix(0,n,ncol(idx))
    for (i in 1:ncol(idx)){
      
      XX[,i]<- P[,c(idx[,i])]%*%s
    }
    data <- XX
    out<-bmt_is(data,0,10^{-6})
    return(out)
  }
  sfExport('n','P','R','cosci','idx')
  out.cosci.2<- sfClusterApplyLB(1:24,wrapper1)
  outcosci.2<- matrix(unlist(out.cosci.2),ncol=ncol(idx),byrow=TRUE)
  sfRemoveAll()
  S[r,]<- c(outcosci.1,apply(outcosci.2,2,max))
  tempdir<- apply(outcosci.2,2,which.max)
  st<-2*(r-1)+1
  en<-2*r
  dir.opt[,st:en]<- R[tempdir,]
  print(r)
}
sfStop()
save.image(paste(wd,'/cosci_2way.RData',sep=""))

heads<- matrix("",1,325)
for (i in 1:325){
  if(i<=25){
    heads[i]<-toString(i)
  }
  if(i>25){
    s<-i-25
  heads[i]<- paste(idx[1,s],idx[2,s],sep="_")
  }
}

alpha<- c(0.05,0.08,0.1,0.12,0.15,0.20,0.25)
sig<- c(1,2,3,4)
noise<- 5:25
FN<- matrix(0,length(alpha),2)
FP<- FN
for (i in 1:length(alpha)){
  results<- matrix(0,reps,2)
  for (r in 1:reps){
    
    st<-2*(r-1)+1
    en<-2*r
    selected<- which(S[r,]>alpha[i])
    features<- heads[selected]
    combs<-selected[selected>25] -25
    combs.features<- features[which(selected>25)]
    dirs<- as.matrix(dir.opt[combs,st:en])
    dirs[abs(dirs[,1])>=dir.threshold,1]<-1
    dirs[dirs[,1]==1,2]<-0
    dirs[abs(dirs[,2])>=dir.threshold,2]<-1
    dirs[dirs[,2]==1,1]<-0
    combs.info<-as.data.frame(cbind(combs.features,transpose(idx[,combs]),
                                    combs+25,dirs))
    colnames(combs.info)<-c("features","feature.x","feature.y","index","dir.x",
                            "dir.y")
    a<- unique(as.numeric(as.character(combs.info$feature.x[abs(dirs[,1])>0])))
    b<- unique(as.numeric(as.character(combs.info$feature.y[abs(dirs[,2])>0])))
    c<- unique(c(a,b))
    selected.final<- unique(c(selected[selected<=25],c))
    results[r,]<- cbind(sum(!(sig%in%selected.final)),
                        sum(noise%in%selected.final))
  }
  
  FN[i,]<- cbind(mean(results[,1]),sd(results[,1])/sqrt(reps))
  FP[i,]<- cbind(mean(results[,2]),sd(results[,2])/sqrt(reps))  
}

tab_cosci<- as.data.frame(cbind(FN,FP))

#---- Data Driven COSCI ----------
scoremat<- S
out.coscidd<- matrix(0,reps,p)
id.coscidd<-out.coscidd
results<- matrix(0,reps,2)
sig<- c(1,2,3,4)
noise<- 5:25
for (i in 1:reps){
  
  score <- scoremat[i,]
  temp<-cosci_2wayselect(score,0.9)
  selected<-temp$selected
  features<- heads[selected]
  combs<-selected[selected>25] -25
  combs.features<- features[which(selected>25)]
  dirs<- as.matrix(dir.opt[combs,st:en])
  dirs[abs(dirs[,1])>=dir.threshold,1]<-1
  dirs[dirs[,1]==1,2]<-0
  dirs[abs(dirs[,2])>=dir.threshold,2]<-1
  dirs[dirs[,2]==1,1]<-0
  combs.info<-as.data.frame(cbind(combs.features,transpose(idx[,combs]),
                                  combs+25,dirs))
  colnames(combs.info)<-c("features","feature.x","feature.y","index","dir.x",
                          "dir.y")
  a<- unique(as.numeric(as.character(combs.info$feature.x[abs(dirs[,1])>0])))
  b<- unique(as.numeric(as.character(combs.info$feature.y[abs(dirs[,2])>0])))
  c<- unique(c(a,b))
  selected.final<- unique(c(selected[selected<=25],c))
  results[i,]<- cbind(sum(!(sig%in%selected.final)),
                      sum(noise%in%selected.final))
}
FN.dd<- cbind(mean(results[,1]),sd(results[,1])/sqrt(reps))
FP.dd<- cbind(mean(results[,2]),sd(results[,2])/sqrt(reps))  



