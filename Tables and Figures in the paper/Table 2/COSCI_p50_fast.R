
wd<- 'Enter your working directory here'
ncpu<-8 #Number of workers
generate.data.first<- TRUE #FALSE if data already exists at datalocation

#-----------------------------------------------------------------------
setwd(wd)
datalocation <- wd
source("coscilibrary.R")
library(sparcl)

sim.grp<- function(x,r){
  n<- nrow(x)
  p<- ncol(x)
  sig<- matrix(0,n,(p+1))
  c1<- c(rep(0,n/2),rep(1,n/2))
  c2<- c1
  c3<- c(rep(0,n/4),rep(1,n/4),rep(2,n/4),rep(3,n/4))
  c4<- c1
  set.seed(r)
  index<-cbind(sample(n),sample(n),sample(n),sample(n))
  c1<- factor(c1[index[,1]])
  c2<- factor(c2[index[,2]])
  c3<- factor(c3[index[,3]])
  c4<- factor(c4[index[,4]])
  grp<- interaction(c1,c2,c3,c4)
  grp<- factor(as.numeric(grp)-1)
  
  sig[,1]<- x[index[,1],1]
  sig[,2]<- x[index[,2],2]
  sig[,3:4]<- x[index[,3],3:4]
  sig[,5]<- x[index[,4],5]
  sig[,6]<- grp
  
  return(sig)
}

generate.data<- function(n,datalocation){
  
  temp<- paste(datalocation,toString(n),sep="/")
  dir.create(temp)
  p<- 50
  sigs<- 5
  noise<- p-sigs
  reps<- 50
  
  # signal parameters
  a1 = 4
  b1 = 6
  a2 = 7
  b2 = 3
  mu1 <- c(0, 0)
  mu2 <- c(0,-4)
  mu3 <- c(4,0)
  mu4 <- c(4,-4)
  rho1 <- -0.85
  rho2 <- 0.85
  sigma1 = matrix(c(1, rho1,rho1,1),2)
  sigma2 = matrix(c(1,rho2,rho2,1),2)
  
  for (r in 1:reps){
    
    name1<- paste('data_',toString(r),'.RData',sep="")
    name2<- paste('data_',toString(r),'.csv',sep="")
    set.seed(10*r)
    x1 <- c(rbeta(n/2, a1, b1, ncp = 0),rbeta(n/2, a2,b2, ncp = 0))
    set.seed(10*r)
    x2<- c(rlnorm(n/2,0.2,0.35),rnorm(n/2,4,0.5))
    set.seed(10*r)
    x31 <- mvrnorm(n/4, mu1, sigma1)
    set.seed(10*r)
    x32 <- mvrnorm(n/4, mu2, sigma2)
    set.seed(10*r)
    y31 <- mvrnorm(n/4, mu3, sigma2)
    set.seed(10*r)
    y32 <- mvrnorm(n/4, mu4, sigma1)
    x3 <- rbind(x31,x32,y31,y32)
    set.seed(10*r)
    x4 <- c(rdoublex(n/2,mu=3,lambda=1.5),rdoublex(n/2,mu=5,lambda=1.5))
    sig<- sim.grp(cbind(x1,x2,x3,x4),r)
    grp<- sig[,6]
    set.seed(10*r)
    y = matrix(rnorm(n*noise,0,1),n,noise)
    d <- cbind(sig[,1:5],y,grp)
    save(d,file=paste(temp,name1,sep="/"))
    write.csv(d, file = paste(temp,name2,sep="/"))
    
  }
}

COSCI_p50<- function(n,datalocation,kmhc){
  
  library(snowfall)  
  sfInit(parallel = TRUE,cpu = ncpu)
  #sfSource('coscilibrary.R')
  
  p<- 50
  sigs<- 5
  reps<- 50
  
  wrapper1<- function(r){
    temp<- paste(datalocation,toString(n),sep="/")
    name<- paste('data_',toString(r),'.RData',sep="")
    load(paste(temp,name,sep="/"))
    grp<- d[,p+1]
    d <- d[,1:p]
    out<-cosci_is(d,0,10^{-6})
    return(out)
  }
  sfExport('n','datalocation','cosci_is')
  outcosci.list<- sfClusterApplyLB(1:reps,wrapper1)
  sfStop()
  outcosci<- matrix(unlist(outcosci.list),ncol=p,byrow=TRUE)
    
  if(kmhc == 1){
    km.weights<-matrix(0,p,10)
    hc.weights<-km.weights
   for (idx in 1:reps){  
      temp<- paste(datalocation,toString(n),sep="/")
      name<- paste('data_',toString(idx),'.RData',sep="")
      load(paste(temp,name,sep="/"))
      grp<- d[,p+1]
      d <- d[,1:p]
      D<- scale(d,TRUE,TRUE)
      rm('d')
      
      km.perm <- KMeansSparseCluster.permute(D,K=32,nperms=5)
      km.out <- KMeansSparseCluster(D,K=32,wbounds=km.perm$bestw)
      km.weights[,idx]<- as.matrix(km.out[[1]]$ws)
      rm('km.perm','km.out')
      
      hc.perm <- HierarchicalSparseCluster.permute(D,nperms=5)
      hcdist<- hc.perm$dists
      hcbestw<- hc.perm$bestw
      rm('hc.perm')
      hc.out <- HierarchicalSparseCluster(dists=hcdist,wbound=hcbestw, method="complete")
      
      hc.weights[,idx]<- as.matrix(hc.out$ws)
      rm('hc.out')
      print(idx)
          
    }
    
    return(list("COSCI"=outcosci,"KM"=km.weights,"HC"=hc.weights))
  }
  if (kmhc == 0){
    return(list("COSCI"=outcosci))
  }
}

COSCI_select_p50<- function(scoremat,gamma){
  
  p<- ncol(scoremat)
  reps<- nrow(scoremat)
  out_select<- matrix(0,reps,p)
  results<- matrix(0,reps,2)
  sig<- c(1,2,3,4,5)
  noise<- 6:p
  for (i in 1:reps){
    
    score <- scoremat[i,]
    temp<-cosci_is_select(score,gamma)
    r<- length(temp$selected)
    print(i)
    out_select[i,1:r]<- temp$selected
    results[i,]<- cbind(sum(!(sig%in%temp$selected)),
                    sum(noise%in%temp$selected))
    
  }
  output<- cbind(mean(results[,1]),sd(results[,1])/sqrt(reps),
                 mean(results[,2]),sd(results[,2])/sqrt(reps))
  return(output)
  
}

tables<- function(output,kmhc){
  reps<-50
  p<-50
  temp3<-matrix(0,reps,1)
  temp4<-matrix(0,reps,1)
  temp5<-matrix(0,reps,1)
  temp6<-matrix(0,reps,1)
  alpha<- c(0.05,0.08,0.1,0.12,0.15,0.20)
  
  if (kmhc == 1){
    bb<-matrix(unlist(output$KM),p,reps)
    cc<-matrix(unlist(output$HC),p,reps)
    for (i in 1:reps){
      temp3[i]<- sum(bb[1:5,i]==0)
      temp4[i]<- sum(bb[6:p,i]>0)
      temp5[i]<- sum(cc[1:5,i]==0)
      temp6[i]<- sum(cc[6:p,i]>0)
    }
    FN.km<- cbind(sum(temp3)/reps,sd(temp3)/sqrt(reps))
    FP.km<- cbind(sum(temp4)/reps,sd(temp4)/sqrt(reps))
    FN.hc<- cbind(sum(temp5)/reps,sd(temp5)/sqrt(reps))
    FP.hc<- cbind(sum(temp6)/reps,sd(temp6)/sqrt(reps))
    
    FN.COSCI<- matrix(0,length(alpha),2)
    FP.COSCI<- FN.COSCI
    for (i in 1:length(alpha)){
      FN.COSCI[i,]<- cbind(sum(rowSums(output$COSCI[,1:5]<=alpha[i]))/reps,
                         sd(rowSums(output$COSCI[,1:5]<=alpha[i]))/sqrt(reps))
      FP.COSCI[i,]<- cbind(sum(rowSums(output$COSCI[,6:p]>alpha[i]))/reps,
                         sd(rowSums(output$COSCI[,6:p]>alpha[i]))/sqrt(reps))
    }
    
    tab<- as.data.frame(cbind(rbind(FN.COSCI,FN.km,FN.hc),
                              rbind(FP.COSCI,FP.km,FP.hc)))
  }
  if (kmhc == 0){
    
    FN.COSCI<- matrix(0,length(alpha),2)
    FP.COSCI<- FN.COSCI
    for (i in 1:length(alpha)){
      FN.COSCI[i,]<- cbind(sum(rowSums(output$COSCI[,1:5]<=alpha[i]))/reps,
                         sd(rowSums(output$COSCI[,1:5]<=alpha[i]))/sqrt(reps))
      FP.COSCI[i,]<- cbind(sum(rowSums(output$COSCI[,6:p]>alpha[i]))/reps,
                         sd(rowSums(output$COSCI[,6:p]>alpha[i]))/sqrt(reps))
    }
    
    tab<- as.data.frame(cbind(FN.COSCI,FP.COSCI))
     
  }
  return(tab)
  
}

#-------------------------------------------------

if (generate.data.first){
  generate.data(200,datalocation)
  generate.data(1000,datalocation)
  generate.data(2500,datalocation)
}
out_200_scores<- COSCI_p50(200,datalocation,1)
tab_200<- tables(out_200_scores,1)
out_200_select<- COSCI_select_p50(out_200_scores$COSCI,0.9)

out_1000_scores<- COSCI_p50(1000,datalocation,1)
tab_1000<- tables(out_1000_scores,1)
out_1000_select<- COSCI_select_p50(out_1000_scores$COSCI,0.9)

out_2500_scores<- COSCI_p50(2500,datalocation,1)
tab_2500<- tables(out_2500_scores,1)
out_2500_select<- COSCI_select_p50(out_2500_scores$COSCI,0.9)





