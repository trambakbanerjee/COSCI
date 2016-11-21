
#--------------------------------------------------
# Memory intensive computation. Please have atleast 8GB of RAM
#-------------------------------------------------------

wd<- 'Enter your working directory here'
ncpu<-8 #Number of workers
generate.data.first<- TRUE #FALSE if data already exists at datalocation

#-----------------------------------------------------------------------
setwd(wd)
datalocation <- wd
source("coscilibrary.R")

generate.data<- function(n,datalocation){
  
  temp<- paste(datalocation,toString(n),sep="/")
  dir.create(temp)
  p<- 5000
  sigs<- 7
  nnoise<- 1500
  tnoise<- 1500
  enoise<- p-sigs-nnoise-tnoise
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
    set.seed(10*r)
    x5 <- c(rnorm(0.3*n,-2.5,1),rnorm(0.3*n,0,1),rnorm((n-0.6*n),2.5,1))
    set.seed(10*r)
    x6 <- c(rnorm(0.5*n,-1.1,1.5),rnorm(0.5*n,1.1,1))
    set.seed(10*r)
    y = matrix(rnorm(n*nnoise,0,1),n,nnoise)
    set.seed(10*r)
    z = matrix(rt(n*tnoise,5),n,tnoise)
    set.seed(10*r)
    e = matrix(rexp(n*enoise,1),n,enoise)
    d <- cbind(x1,x2,x3,x4,x5,x6,y,z,e)
    save(d,file=paste(temp,name1,sep="/"))
    write.csv(d, file = paste(temp,name2,sep="/"))
    
  }
}

COSCI_p5k<- function(n,datalocation){
  
  library(snowfall)  
  sfInit(parallel = TRUE,cpu = ncpu)
  #sfSource('coscilibrary.R')
  
  p<- 5000
  sigs<- 7
  reps<- 50
  
  wrapper1<- function(r){
    temp<- paste(datalocation,toString(n),sep="/")
    name<- paste('data_',toString(r),'.RData',sep="")
    load(paste(temp,name,sep="/"))
    d <- d[,1:p]
    out<-cosci_is(d,0,10^{-6})
    return(out)
  }
  sfExport('n','datalocation','cosci_is')
  outcosci.list<- sfClusterApplyLB(1:reps,wrapper1)
  sfStop()
  outcosci<- matrix(unlist(outcosci.list),ncol=p,byrow=TRUE)
  return(list("COSCI"=outcosci))
  
}

COSCI_select_p5k<- function(scoremat,gamma){
  
  p<- ncol(scoremat)
  reps<- nrow(scoremat)
  out_select<- matrix(0,reps,p)
  results<- matrix(0,reps,2)
  sig<- c(1,2,3,4,5,6,7)
  noise<- 8:p
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

tables<- function(output){
  reps<-50
  p<-5000
  temp3<-c()
  temp4<-c()
  temp5<-c()
  temp6<-c()
  alpha<- c(0.05,0.08,0.1,0.12,0.15,0.20)
  
    FN.COSCI<- matrix(0,length(alpha),2)
    FP.COSCI<- FN.COSCI
    for (i in 1:length(alpha)){
      FN.COSCI[i,]<- cbind(sum(rowSums(output$COSCI[,1:7]<=alpha[i]))/reps,
                         sd(rowSums(output$COSCI[,1:7]<=alpha[i]))/sqrt(reps))
      FP.COSCI[i,]<- cbind(sum(rowSums(output$COSCI[,8:p]>alpha[i]))/reps,
                         sd(rowSums(output$COSCI[,8:p]>alpha[i]))/sqrt(reps))
    }
    
    tab<- as.data.frame(cbind(FN.COSCI,FP.COSCI))
    
  return(tab)
  
}
#--------------------------------------------------------------------------------------------

if (generate.data.first){
  generate.data(200,datalocation)
  generate.data(1000,datalocation)
  generate.data(2500,datalocation)
}
ptm<-proc.time()
out_200_scores<- COSCI_p5k(200,datalocation)
proc.time()-ptm
tab_200_scores<- tables(out_200_scores)
out_200_select<- COSCI_select_p5k(out_200_scores$COSCI,0.9)

out_1000_scores<- COSCI_p5k(1000,datalocation)
tab_1000_scores<- tables(out_1000_scores)
out_1000_select<- COSCI_select_p5k(out_1000_scores$COSCI,0.9)

out_2500_scores<- COSCI_prodstruct_p5k(2500,datalocation)
tab_2500_scores<- tables(out_2500_scores)
out_2500_select<- COSCI_select_p5k(out_2500_scores$COSCI,0.9)





