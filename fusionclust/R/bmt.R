
#' Big Merge Tracker
#'
#' Solves an L1 relaxed univariate clustering criterion and returns a
#'     sequence of \eqn{\lambda} values where the clusters merge
#'
#' @importFrom stats rnorm
#'
#' @param x observation vector
#' @param alpha merging threshold. Default is 0.1
#' @param small.perturbation a small positive number to remove ties. Default is 10^(-6)
#'
#' @return
#' \enumerate{
#' \item path - number of clusters on the big merge path
#' \item lambda.path - sequence of lambda where clusters merge
#' \item index - cluster index at the point where clusters merge
#' \item merge - merge points
#' \item split - split points
#' \item prob - merging proportion
#' \item boundaries - cluster boundaries
#' }
#'
#' @details solves a convex relaxation of the univariate clustering criterion given by equation
#'     (2) in the referenced paper and generates a sequence of cluster merges and corresponding
#'      \eqn{\lambda} values. See algorithm 1 in the referenced paper for more details.
#'
#' @seealso \code{\link{nclust}}
#'
#' @examples
#' \donttest{
#' x<- c(rnorm(500,-2,1), rnorm(500,2,1))
#' out<- bmt(x)
#' }
#'
#' @references
#' \enumerate{
#' \item P. Radchenko, G. Mukherjee, Convex clustering via l1 fusion penalization,
#'  J. Roy. Statist, Soc. Ser. B (Statistical Methodology) (2017)
#'  doi:10.1111/rssb.12226.
#' }
#'
#' @export

bmt<-function(x,alpha=0.10,small.perturbation=10^(-6))
{
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
      print(paste('BEWARE: there is a big merge when total no. of clusters are ',C));
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
