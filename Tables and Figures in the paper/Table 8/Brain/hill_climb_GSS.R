#This version is the same as hill_climb_GSS_revised.R in the folder of "Discarded"
require(sparcl)
require(mclust)
require(optimbase)
require(phyclust)

wd<-'Enter your working directory here'
#-------------------------------------------------------

withinss = function(X,G,K){
  if(is.vector(X) && is.atomic(X)){
    group = list()
    wcss = 0
    for(j in 1:K){
      group[[j]] = which(G == j)
      l = length(group[[j]])
      if(l>1){
        wcss = wcss + var(X[group[[j]]])*(l-1)
      }
    }
  } else if (is.matrix(X)){
    group = list()
    wcss = 0
    for(j in 1:K){
      group[[j]] = which(G == j)
      l = length(group[[j]])
      if(l>1){
        wcss = wcss + sum(apply(X[group[[j]],],2,var)*(l-1))
      }   
    }
  } else {
    cat("X is niether a vector nor a matrix! \n")
    return(NULL)
  } 
  
  return(wcss)
}

Alternate= function(X, k,tot, initial_set, s, itermax, threshold){
  p = dim(X)[2]
  set0 = initial_set
  set1 = c(0)
  iternum = 0
  while(iternum<= itermax && length(setdiff(set1,set0)) + length(setdiff(set0,set1)) > threshold ){
    clustering = kmeans(X[,set0],iter.max = 20, centers=k,algorithm = "Hartigan-Wong",trace = 0,nstart=2)
    result = clustering$cluster
    
    wcss = apply(X,2,withinss,G = result, K = k)
    iternum = iternum + 1
    
    set1 = set0
    set0 = which(rank((tot-wcss)/tot,ties.method = "random") > p-s)
  }
  out = list(final_set = set0, iternum = iternum, result = result, betweenss = clustering$betweenss)
  return(out)
}

#compute within-cluster distance by clustering feature by feature, select S of size s based on this
hill_climb_GSS = function(X,k,nperms=20,itermax,threshold,tolerance){
  X = scale(X)
  n = dim(X)[1]
  p = dim(X)[2]
  tot = apply(X,2,var)*(n-1)
  permx <- list()
  for(i in 1:nperms){
    permx[[i]] <- matrix(NA, nrow=n, ncol=p)
    for(j in 1:p) permx[[i]][,j] <- sample(X[,j])
  }
  wcss = rep(0,p)
  for(j in 1:p){
    clustering = kmeans(X[,j],iter.max = 10, centers = k,algorithm = "Hartigan-Wong", trace = 0)
    wcss[j] = clustering$tot.withinss
  }
  rank0 = rank((tot-wcss)/tot,ties.method = "random")
  golden.ratio = 2/(sqrt(5) +1)
  iteration = 0
  upper.bound = p
  lower.bound = 1
  p1 = floor(upper.bound - golden.ratio*(upper.bound-lower.bound))
  p2 = floor(lower.bound + golden.ratio*(upper.bound-lower.bound))
  #evaluate the gap statistics using p1 and p2
  initial_set = which(rank0 > p-p1)
  out1 = Alternate(X, k,tot, initial_set, p1, itermax, threshold)
  permtots = rep(0,nperms)
  for(t in 1:nperms){
    permresult = kmeans(permx[[t]][,out1$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
    permtots[t] <- permresult$betweenss
  }
  gap1 = log(out1$betweenss) - mean(log(permtots))
  initial_set = which(rank0 > p-p2)
  out2 = Alternate(X, k,tot, initial_set, p2, itermax, threshold)
  permtots = rep(0,nperms)
  for(t in 1:nperms){
    permresult = kmeans(permx[[t]][,out2$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
    permtots[t] <- permresult$betweenss
  }
  gap2 = log(out2$betweenss) - mean(log(permtots))  
  
  while(abs(upper.bound - lower.bound) > tolerance)
  {
    iteration = iteration + 1
    if(gap2 < gap1) # then the maximum is to the left of x2
    {
      upper.bound = p2
      p2 = p1
      gap2 = gap1
      p1 = floor(upper.bound - golden.ratio*(upper.bound - lower.bound))
      #evaluate gaps for p1
      initial_set = which(rank0 > p-p1)
      out1 = Alternate(X, k,tot, initial_set, p1, itermax, threshold)
      
      permtots = rep(0,nperms)
      for(t in 1:nperms){
        permresult = kmeans(permx[[t]][,out1$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
        permtots[t] <- permresult$betweenss
      }
      gap1 = log(out1$betweenss) - mean(log(permtots))
    } else {
      # the minimum is to the right of x1
      lower.bound = p1
      p1 = p2
      gap1 = gap2
      p2 = floor(lower.bound + golden.ratio * (upper.bound - lower.bound))
      #evaluate gaps for p2
      initial_set = which(rank0 > p-p2)
      out2 = Alternate(X, k,tot, initial_set, p2, itermax, threshold)
      
      permtots = rep(0,nperms)
      for(t in 1:nperms){
        permresult = kmeans(permx[[t]][,out2$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
        permtots[t] <- permresult$betweenss
      }
      gap2 = log(out2$betweenss) - mean(log(permtots))
    }   
  }
  s = floor((lower.bound + upper.bound)/2)
  initial_set = which(rank0 > s)
  out = Alternate(X, k,tot, initial_set, s, 2*itermax, threshold)
  output = list(final_set = out$final_set, iternum = iteration, result = out$result, s = s)
  return(output)
}
#------------------------------------------------------
setwd(wd)
d<- as.matrix(read.csv('brain_x.txt',header=F,sep=' '))
true.class<- as.matrix(read.csv('brain_y.txt',header=F,sep=' '))
nclass<- length(unique(true.class))
d<- as.matrix(transpose(d))
d<- scale(d,TRUE,TRUE)
p<-ncol(d)
n<- nrow(d)

out <- matrix(0,1,p)
fit2 = hill_climb_GSS(d,k=nclass,nperms=25,itermax=20,threshold=0,tolerance=0)
out[fit2$final_set]<- 1

predclass<- fit2$result
Q<- as.matrix(true.class)
m<-n
num<-0
den<- m*(m-1)/2
for (ii in 1:(m-1)){
  
  tt<- predclass[ii]-predclass[(ii+1):m]
  tt[tt!=0]=1
  tt <- 1-tt
  dd<- Q[ii]-Q[(ii+1):m]
  dd[dd!=0]=1
  dd <- 1-dd
  num<- sum(abs(tt-dd))+num
  
}
CER.SAS<- num/den
