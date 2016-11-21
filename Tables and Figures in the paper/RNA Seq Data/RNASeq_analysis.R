require(ggplot2)
require(bbmle)
library(gridExtra)

wd<- 'Enter your working directory here'
#-----------------------------------------------------------
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
    idx <- c(score[,1]>min.alpha & score[,2]>0)
    grp[idx]<-score[idx,1]
    m<- max(grp)
    l<- min(grp)
    dim.score[j,]<- cbind(max(grp[idx]),mean(grp[idx]),sum(idx))
    print(j)
    
  }
  dim.score[is.nan(dim.score)]<-0
  
  return(dim.score)
}

setwd(wd)
load(paste(wd,"/RNASeq.RData",sep=""))
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
true.index<-c(1760,5937,6001,4628,1901,2201,1900,6803,4038,4317,2788,6005,
              4768,3179,3922,4362,3628,1827,5582,5479,6770,1032,3298,3297,3336,1631,
              1635,4196,8634,2157,1455,3672,3671)

score<-as.matrix(cosci(tdata,10^{-5},10^{-4}))
#-------------------------------------------
score.raw<- score[,1]
x<-score.raw/0.5
save.image(paste(wd,'/RNASeq_cosci.RData',sep=""))
#-------------------------------------------------------------------
load(paste(wd,'/RNASeq_cosci.RData',sep=""))

ux <- unique(x)
mode<- ux[which.max(tabulate(match(x, ux)))]
see<- order(x)
interval<- c(0,x[see[ceiling(0.7*p)]])
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

histinfo<-hist(x,breaks=150)
xx<- histinfo$mids
X<- cbind(xx,xx^{2},xx^{3},xx^{4},xx^{5},xx^{6},xx^{7})
y<- histinfo$counts
m<-glm(y~X,family="poisson")
beta<-as.vector(m$coefficients)
beta<- matrix(beta,ncol=1,byrow=TRUE)
X<-cbind(matrix(1,p,1),x,x^{2},x^{3},x^{4},x^{5},x^{6},x^{7})
v<- exp(X%*%beta)
fhat<- v/(length(x)*0.01)

fdr<- (pi0*f.null)/fhat
fdr[fdr>1]=1
plot(x,fdr)
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
select.d1<- tdata[,idx.1]
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
write.csv('select.d2',paste(wd,'/selected.csv',sep=""))

fdr.2<- fdr.1[order(fdr.1)]
cutoff.3<- fdr.2[70]
idx.3<- (fdr.1<=cutoff.3)
select.d3<- select.d1[,idx.3]
alpha.3<-alpha.1[idx.3]

#------------------------- get some plots --------------------------

d0 <- data.frame(x)
d1<-as.data.frame(cbind(x,(pi0*f.null)))
colnames(d1)<-c("X","f0")
d2<-as.data.frame(cbind(x,fhat))
colnames(d2)<-c("X","fhat")
g1<-ggplot(d0, aes(x = x)) + geom_histogram(aes(y = ..density..),binwidth=0.05,color='black',fill='white') +
  geom_line(data=d1,(aes(x=X,y=f0)),color='blue',size=1) + 
  geom_line(data=d2,(aes(x=X,y=fhat)),color='red',size=1)+
  xlab("2*Sj")+ylab("density of 2*Sj")

dframe1<- as.data.frame(cbind(1:p,score.raw[order(score.raw)]))
colnames(dframe1)<- c("x","S_j")
idx<-sum(score.raw<=min(alpha.2))
dframe31<- as.data.frame(cbind(1:idx,min(alpha.2)*rep(1,idx)))
colnames(dframe31)<-c("X","Y")
dframe32<- as.data.frame(cbind(idx*rep(1,2),c(0,min(alpha.2))))
colnames(dframe32)<-c("X","Y")
idx<- score.raw[score.raw<min(alpha.2)]
y<- 1:length(idx)
dframe4<-as.data.frame(cbind(y,idx[order(idx)]))
colnames(dframe4)<-c("X","Y")

lll<- order(score.raw)
ggg<- matrix(0,33,1)
for (i in 1:33){
  ggg[i]=which(true.index[i]==lll)
}
dframe5<- as.data.frame(cbind(ggg,score.raw[lll[ggg]]))
colnames(dframe5)<- c("X","Y")

g2<-ggplot(dframe1,aes(x=x,y=S_j)) + geom_point(color='gray60',size=0.5)+
  geom_line(data=dframe31,aes(x=X,y=Y),color='black',size=0.3)+
  geom_line(data=dframe32,aes(x=X,y=Y),color='black',size=0.3)+
  geom_point(data=dframe4,aes(x=X,y=Y),color='red',size=0.5)+
  geom_point(data=dframe5,aes(x=X,y=Y),color='blue',size=3)+
  xlab("Coordinates")+ylab("Sj")+
  theme(panel.background=element_rect(fill="white"),
        axis.line=element_line(color="black"))

grid.arrange(g1,g2,ncol=2)

#-------------------------COSCI + SpKM and COSCI + KM-----------------

require(sparcl)
require(gmodels)

D<- scale(select.d2,TRUE,TRUE)

km.perm <- KMeansSparseCluster.permute(D,K=nclass,nperms=25)
km.out <- KMeansSparseCluster(D,K=nclass,wbounds=km.perm$bestw)
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
  
  CER.COSCIVKM[i]<- num/den
  print(i)
}
avgCER.COSCIKM<- c(mean(CER.COSCIVKM),sd(CER.COSCIVKM)/sqrt(30))

#------- Heatmap ----------------------
set.seed(42)
vkm<- kmeans(D,nclass, iter.max = 50, nstart = 30,algorithm = "Hartigan-Wong")
predclass.km<- vkm$cluster
orderpredclass.km<- predclass.km[order(predclass.km)]
ordereddframe5<- dframe5[order(dframe5$Y),]
ordereddframe5<- as.numeric(row.names(ordereddframe5))
lineage<- as.character(genes[true.index[ordereddframe5],])
lineage[16]<-"Serpina3f"
lineage[17]<- "Pbx1"
classticks<- matrix(0,18,2)
for (iticks in 1:18){
  classticks[iticks,1]<- min(which(orderpredclass.km==(iticks+1)))
}
for (i in 1:18){
  if (i>1){
    classticks[i,2]<- ceiling(0.5*(classticks[i,1]+classticks[(i-1),1]))
  }
  if (i==1){
    classticks[i,2]<- ceiling(0.5*(classticks[i,1]+1))
  }
}
selected<- cbind(select.d2,predclass.km)
selected<-selected[order(selected[,2305]),]
colors = c(seq(-5,-4,length=2),seq(-4,-2,length=4),seq(-2,0,length=4))
mycolor<- colorRampPalette(c("bisque","gold","red"))(9);
image(1:dim(selected)[1],1:33,tdata[orderpredclass.km,true.index[ordereddframe5]],
      xlab='Cells',ylab='',col=mycolor,xaxt='n',yaxt='n',breaks=colors);
vv<- as.character(c(rep("black",21),rep("blue",12)))
axis(2, at=seq(1,21,1),labels=lineage[1:21],line = -0.5, tick = FALSE,cex.axis=0.8,
     las=2,col.axis='black')
axis(2, at=seq(22,33,1),labels=lineage[22:33],line = -0.5, tick = FALSE,cex.axis=0.8,
     las=2,col.axis='blue')
axis(1, at=c(classticks[,2],0.5*(n+classticks[18,1])),
     labels=seq(1,19,1),
     line = -0.5, tick = FALSE,cex.axis=0.8)
abline(v=classticks[,1],col='black',lwd=2)

image(1:dim(selected)[1],1:8716,tdata[,],
      xlab='Cells',ylab='',col=mycolor,xaxt='n',yaxt='n',breaks=colors);
vv<- as.character(c(rep("black",21),rep("blue",12)))
axis(2, at=seq(1,21,1),labels=lineage[1:21],line = -0.5, tick = FALSE,cex.axis=0.8,
     las=2,col.axis='black')
axis(2, at=seq(22,33,1),labels=lineage[22:33],line = -0.5, tick = FALSE,cex.axis=0.8,
     las=2,col.axis='blue')
axis(1, at=c(classticks[,2],0.5*(n+classticks[18,1])),
     labels=seq(1,19,1),
     line = -0.5, tick = FALSE,cex.axis=0.8)
abline(v=classticks[,1],col='black',lwd=2)


#------------------------- Full SpKM and KM --------------------------
# This is a time intensive process -------------------------------------

D<- scale(tdata,TRUE,TRUE)

km.perm <- KMeansSparseCluster.permute(D,K=nclass,nperms=2,nvals=2)
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
# --------------------------------------------------



