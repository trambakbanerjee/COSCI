

#' Rank the p features in an n by p design matrix
#'
#' Ranks the p features in an n by p design matrix
#' where n represents the sample size and p is the number of features.
#'
#' @param dat n by p data matrix
#' @param min.alpha the smallest threshold (typically set to 0)
#' @param small.perturbation a small positive number to remove ties. Default value is 10^(-6)
#'
#' @return  a p vector of scores
#'
#' @details Uses the univariate merging algorithm \code{\link{bmt}} and produces a score
#'     for each feature that reflects its relative importance for clustering.
#'
#' @seealso \code{\link{bmt}},\code{\link{cosci_is_select}}
#'
#' @examples
#' \donttest{
#' x<-matrix(rnorm(10000),nrow=500,ncol=20)
#' s<- cosci_is(x,0)
#' }
#'
#' @references
#' \enumerate{
#' \item  Banerjee, T., Mukherjee, G. and Radchenko P., Feature Screening in
#' Large Scale Cluster Analysis, Journal of Multivariate Analysis,
#' Volume 161, 2017, Pages 191-212
#' \item P. Radchenko, G. Mukherjee, Convex clustering via l1 fusion penalization,
#'  J. Roy. Statist, Soc. Ser. B (Statistical Methodology) (2017)
#'  doi:10.1111/rssb.12226.
#' }
#'
#' @export

cosci_is<-function(dat,min.alpha,small.perturbation=10^(-6)){

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

  }
  dim.score[is.nan(dim.score)]<-0

  return(dim.score[,1])
}

#' Use a data driven approach to select the features
#'
#' Once you have the feature scores from \code{\link{cosci_is}}, you can select the features
#' \enumerate{
#' \item based on a pre-defined threshold,
#' \item using table A.10 in the paper[1] to determine an appropriate threshold or,
#' \item using a data driven approach described in the references to select the features
#' and obtain an implicit threshold value.
#' }
#' cosci_is_select implements option 3.
#'
#' @import bbmle
#' @importFrom graphics hist
#' @importFrom stats dbeta
#' @importFrom stats glm
#' @importFrom stats pbeta
#' @importFrom stats rnorm
#'
#' @param score a p vector of scores
#' @param gamma what proportion of the p features is noise? If your sample size n
#' is smaller than 100, setting gamma = 0.85 is recommended. Otherwise set gamma = 0.9
#'
#' @return  a vector of selected features
#'
#' @details Converts the problem of screening out features with lower scores into a
#'     problem in large scale multiple testing and uses the procedure described in
#'     reference [2] to determine the signal features.
#'
#' @seealso \code{\link{cosci_is}}
#'
#' @examples
#' \dontrun{
#' x<-rbeta(1000,0.5,5)
#' s<-cosci_is_select(x,0.9)
#' }
#'
#' @references
#' \enumerate{
#' \item  Banerjee, T., Mukherjee, G. and Radchenko P., Feature Screening in
#' Large Scale Cluster Analysis, Journal of Multivariate Analysis,
#' Volume 161, 2017, Pages 191-212
#' \item T. Cai, W. Sun, W., Optimal screening and discovery of sparse signals
#' with applications to multistage high throughput studies,
#' J. Roy.Statist. Soc. Ser. B (Statistical Methodology) 79, no. 1 (2017) 197-223
#' }
#'
#' @export

cosci_is_select<-function(score,gamma){

  warnings('off')
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

  histinfo<-hist(x,breaks=min(p/2,150),plot=FALSE)
  xx<- histinfo$mids
  X<- cbind(xx,xx^{2},xx^{3},xx^{4},xx^{5})
  y<- histinfo$counts
  m<-glm(y~X,family="poisson")
  beta<-as.vector(m$coefficients)
  beta<- matrix(beta,ncol=1,byrow=TRUE)
  X<-cbind(matrix(1,p,1),x,x^{2},x^{3},x^{4},x^{5})
  v<- exp(X%*%beta)
  gap<-histinfo$breaks[2]-histinfo$breaks[1]
  fhat<- v/(length(x)*gap)

  fdr<- (pi0*f.null)/fhat
  fdr[fdr>1]=1
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

