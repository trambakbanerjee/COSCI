
wd<- 'Enter your working directory here'
#--------------------------------------------
setwd(wd)
datalocation <- wd


n<-200
reps<-50

for (i in 1:reps){
  datafile<- paste(datalocation,'/',toString(n),'/data_',toString(i),'.RData',sep="")
  load(datafile)
  name<- paste(datalocation,'/',toString(n),'/data_',toString(i),'.csv',sep="")
  write.csv(d, file = name)
}

n<-1000
reps<-10

for (i in 1:reps){
  datafile<- paste(datalocation,'/',toString(n),'/data_',toString(i),'.RData',sep="")
  load(datafile)
  name<- paste(datalocation,'/',toString(n),'/data_',toString(i),'.csv',sep="")
  write.csv(d, file = name)
}

n<-2500
reps<-15

for (i in 1:reps){
  datafile<- paste(datalocation,'/',toString(n),'/data_',toString(i),'.RData',sep="")
  load(datafile)
  name<- paste(datalocation,'/',toString(n),'/data_',toString(i),'.csv',sep="")
  write.csv(d, file = name)
}