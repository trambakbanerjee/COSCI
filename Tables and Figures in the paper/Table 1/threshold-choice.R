
#############################################################################
# This script generates Table 1 in the paper.                               #
#                                                                           #
# Please provide the working directory name in line 8 below.                #
#############################################################################

setwd('/home/sukumar/Documents/Study/Research/BMT Extension/Simulations in paper/simtab1 in paper')

#-----------------------------------------------------------------------------------------------------
cluster.path.faster<-function(x,alpha=0.10,small.perturbation=10^{-6})
{
#inputs:- observation vector: x; 
#alpha :  merging proportion level threshold; merges with prob of both clusters above alpha are detected; 
#returns the sequence of lambda where clusters merge along with cluster index at that point
    
    #remove ties
    x<-x+rnorm(length(x))*small.perturbation
    
	x<-sort(x);
	n=length(x);
	#initializiations 
	C=n;#representing number of clusters by C
	ind<-1:C;#cluster index of observations
	ind.along.path<-c();  
	l.old<-0;
	lambda<-0;
	split.points<-c();
	merge.points<-c();
	split.points<-c();
	no.of.clusters<-c();
	probability<-c();
	boundaries<-c();
	freq<-rep(1,n);
    centroids<-x;
 	unique.ind<-ind;
 	counter=3;
 	while(TRUE){
	 	l.new<-(centroids[2:C]-centroids[1:(C-1)])/(freq[2:C]+freq[1:(C-1)]);
	 	l.new.mod<-l.new[l.new>=l.old];
		l.old<-min(l.new.mod);
		i=match(l.old,l.new);
		if (min(freq[i],freq[i+1])>alpha*n){
			#print(paste('BEWARE: there is a big merge when total no. of clusters are ',C));
			merge.points<-c(merge.points,C);
			a.2<-min(x[(ind==unique.ind[i+1])]);
			a.1<-max(x[(ind==unique.ind[i])])
			split.points<-c(split.points,(freq[i+1]*a.2+freq[i]*a.1)/(freq[i+1]+freq[i]));
			probability<-rbind(probability,c(freq[i],freq[i+1])/n);
			boundaries<-rbind(boundaries,c(min(x[(ind==unique.ind[i])]),max(x[(ind==unique.ind[i])]),min(x[(ind==unique.ind[i+1])]),max(x[(ind==unique.ind[i+1])])));
			counter<-0;
		}
		counter<-counter+1;
		if (counter < 3){
			ind.along.path<-rbind(ind.along.path,ind);
			no.of.clusters<-c(no.of.clusters,C);
		}		
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
				l<-list(path=no.of.clusters,lambda.path=lambda,index=ind.along.path,merge=merge.points,splits=split.points,prob=probability,boundaries=boundaries);
				return(l);
		}
	}
}

cluster.path.faster.1<-function(x1,alpha)
{
	c<-cluster.path.faster(x1,alpha);
	l<-c$splits;
	if(length(l)<1){return(1);}
	kk<-dim(c$prob)[1];
	if (sum(c$prob[kk,])<0.5){return(1);}
	return(length(l)+1)	
}

gen.mix.normal<-function(n,m=c(-20,0,20),p=rep(1,length(m))/length(m))
{
	x<-c();
	n1<-floor(n*p);
	if(n-sum(n1)>0){n1[1]<-n1[1]+n-sum(n1);}
	for (i in 1:length(m)){x<-c(x,rnorm(n1[i])+m[i]);}
	return(x);
}


#------------------------------------------------------------------------------------------------------
n1=c(100,500,1000,2000,5000,10000)

I=100;
result=matrix(0,I,length(n1))
result.mod=matrix(0,I,length(n1))
answer=c();
answer.mod=c();	

for (j in 1:6){
	n=n1[j]
	for (i in 1:I)
	{
		x1=gen.mix.normal(n,m=c(0,0,0),p=c(0.3,0.35,0.35));
		ans=cluster.path.faster(x1,0.001)
		result[i,j]=max(apply(ans$prob,1,min));
		temp=(apply(ans$prob,1,sum)>0.5)
		result.mod[i,j]=max(apply(ans$prob,1,min)*temp);
		print(i)
	}
	answer=cbind(answer,summary(result[,j]));
	answer.mod=cbind(answer.mod,summary(result.mod[,j]));
	
	print(j)
}

save.image("answer-normal.RData")

n1=c(100,500,1000,2000,5000,10000)
result=matrix(0,I,length(n1))
result.mod=matrix(0,I,length(n1))
answer=c();
answer.mod=c();	

for (j in 1:6){
	n=n1[j]
	for (i in 1:I)
	{
		x1=rt(n,df=1);
		ans=cluster.path.faster(x1,0.001)
		result[i,j]=max(apply(ans$prob,1,min));
		temp=(apply(ans$prob,1,sum)>0.5)
		result.mod[i,j]=max(apply(ans$prob,1,min)*temp);
		print(i)
	}
	answer=cbind(answer,summary(result[,j]));
	answer.mod=cbind(answer.mod,summary(result.mod[,j]));	
	print(j)
}

save.image("answer-t.RData")

n1=c(100,500,1000,2000,5000,10000)
result=matrix(0,I,length(n1))
result.mod=matrix(0,I,length(n1))
answer=c();
answer.mod=c();	

for (j in 1:6){
	n=n1[j]
	for (i in 1:I)
	{
		x1=rexp(n);
		ans=cluster.path.faster(x1,0.001)
		result[i,j]=max(apply(ans$prob,1,min));
		temp=(apply(ans$prob,1,sum)>0.5)
		result.mod[i,j]=max(apply(ans$prob,1,min)*temp);
		print(i)
	}
	answer=cbind(answer,summary(result[,j]));
	answer.mod=cbind(answer.mod,summary(result.mod[,j]));	
	print(j)
}

save.image("answer-exp.RData")


n1=c(100,500,1000,2000,5000,10000)
result=matrix(0,I,length(n1))
result.mod=matrix(0,I,length(n1))
answer=c();
answer.mod=c();	

for (j in 1:6){
	n=n1[j]
	for (i in 1:I)
	{
		x1=rcauchy(n);
		ans=cluster.path.faster(x1,0.001)
		result[i,j]=max(apply(ans$prob,1,min));
		temp=(apply(ans$prob,1,sum)>0.5)
		result.mod[i,j]=max(apply(ans$prob,1,min)*temp);
		print(i)
	}
	answer=cbind(answer,summary(result[,j]));
	answer.mod=cbind(answer.mod,summary(result.mod[,j]));	
	print(j)
}

save.image("answer-cauchy.RData")

require(smoothmest)
n1=c(100,500,1000,2000,5000,10000)
result=matrix(0,I,length(n1))
result.mod=matrix(0,I,length(n1))
answer=c();
answer.mod=c();	

for (j in 1:6){
  n=n1[j]
  for (i in 1:I)
  {
    x1=rdoublex(n,0,1);
    ans=cluster.path.faster(x1,0.001)
    result[i,j]=max(apply(ans$prob,1,min));
    temp=(apply(ans$prob,1,sum)>0.5)
    result.mod[i,j]=max(apply(ans$prob,1,min)*temp);
    print(i)
  }
  answer=cbind(answer,summary(result[,j]));
  answer.mod=cbind(answer.mod,summary(result.mod[,j]));	
  print(j)
}
save.image("answer-dexp.RData")

require(fExtremes)
n1=c(100,500,1000,2000,5000,10000)
result=matrix(0,I,length(n1))
result.mod=matrix(0,I,length(n1))
answer=c();
answer.mod=c();	

for (j in 1:6){
  n=n1[j]
  for (i in 1:I)
  {
    x1=rgev(n,xi = 0.8,mu = 0,beta = 1);
    ans=cluster.path.faster(x1,0.001)
    result[i,j]=max(apply(ans$prob,1,min));
    temp=(apply(ans$prob,1,sum)>0.5)
    result.mod[i,j]=max(apply(ans$prob,1,min)*temp);
    print(i)
  }
  answer=cbind(answer,summary(result[,j]));
  answer.mod=cbind(answer.mod,summary(result.mod[,j]));	
  print(j)
}
save.image("answer-gev.RData")

n1=c(100,500,1000,2000,5000,10000)
result=matrix(0,I,length(n1))
result.mod=matrix(0,I,length(n1))
answer=c();
answer.mod=c();	

for (j in 1:6){
  n=n1[j]
  for (i in 1:I)
  {
    x1=rbeta(n,1,3);
    ans=cluster.path.faster(x1,0.001)
    result[i,j]=max(apply(ans$prob,1,min));
    temp=(apply(ans$prob,1,sum)>0.5)
    result.mod[i,j]=max(apply(ans$prob,1,min)*temp);
    print(i)
  }
  answer=cbind(answer,summary(result[,j]));
  answer.mod=cbind(answer.mod,summary(result.mod[,j]));	
  print(j)
}
save.image("answer-beta13.RData")

require(triangle)
n1=c(100,500,1000,2000,5000,10000)
result=matrix(0,I,length(n1))
result.mod=matrix(0,I,length(n1))
answer=c();
answer.mod=c();	

for (j in 1:6){
  n=n1[j]
  for (i in 1:I)
  {
    x1=rtriangle(n,0,1,c=0.8);
    ans=cluster.path.faster(x1,0.001)
    result[i,j]=max(apply(ans$prob,1,min));
    temp=(apply(ans$prob,1,sum)>0.5)
    result.mod[i,j]=max(apply(ans$prob,1,min)*temp);
    print(i)
  }
  answer=cbind(answer,summary(result[,j]));
  answer.mod=cbind(answer.mod,summary(result.mod[,j]));	
  print(j)
}
save.image("answer-triangle.RData")

#---------------------------------------------------------


load("answer-normal.RData")

alpha=c(1,2,5,7.5,10,15,20,25)/100;
alpha.ans.2=c();

for (i in 1:length(alpha)){
alpha1=alpha[i];
alpha.ans.2=cbind(alpha.ans.2,apply((result.mod>alpha1),2,sum)/100*100)
}

alpha.ans.norm<- alpha.ans.2

load("answer-t.RData")

alpha=c(1,2,5,7.5,10,15,20,25)/100;
alpha.ans.2=c();

for (i in 1:length(alpha)){
alpha1=alpha[i];
alpha.ans.2=cbind(alpha.ans.2,apply((result.mod>alpha1),2,sum)/100*100)
}

alpha.ans.t<-alpha.ans.2

load("answer-exp.RData")

alpha=c(1,2,5,7.5,10,15,20,25)/100;
alpha.ans.2=c();

for (i in 1:length(alpha)){
alpha1=alpha[i];
alpha.ans.2=cbind(alpha.ans.2,apply((result.mod>alpha1),2,sum)/100*100)
}

alpha.ans.exp<-alpha.ans.2

load("answer-cauchy.RData")

alpha=c(1,2,5,7.5,10,15,20,25)/100;
alpha.ans.2=c();

for (i in 1:length(alpha)){
  alpha1=alpha[i];
  alpha.ans.2=cbind(alpha.ans.2,apply((result.mod>alpha1),2,sum)/100*100)
}

alpha.ans.cauchy<-alpha.ans.2

load("answer-dexp.RData")

alpha=c(1,2,5,7.5,10,15,20,25)/100;
alpha.ans.2=c();

for (i in 1:length(alpha)){
alpha1=alpha[i];
alpha.ans.2=cbind(alpha.ans.2,apply((result.mod>alpha1),2,sum)/100*100)
}

alpha.ans.dexp<-alpha.ans.2

load("answer-gev.RData")

alpha=c(1,2,5,7.5,10,15,20,25)/100;
alpha.ans.2=c();

for (i in 1:length(alpha)){
  alpha1=alpha[i];
  alpha.ans.2=cbind(alpha.ans.2,apply((result.mod>alpha1),2,sum)/100*100)
}

alpha.ans.gev<-alpha.ans.2

load("answer-beta13.RData")

alpha=c(1,2,5,7.5,10,15,20,25)/100;
alpha.ans.2=c();

for (i in 1:length(alpha)){
  alpha1=alpha[i];
  alpha.ans.2=cbind(alpha.ans.2,apply((result.mod>alpha1),2,sum)/100*100)
}

alpha.ans.beta13<-alpha.ans.2

load("answer-triangle.RData")

alpha=c(1,2,5,7.5,10,15,20,25)/100;
alpha.ans.2=c();

for (i in 1:length(alpha)){
  alpha1=alpha[i];
  alpha.ans.2=cbind(alpha.ans.2,apply((result.mod>alpha1),2,sum)/100*100)
}

alpha.ans.triang<-alpha.ans.2
