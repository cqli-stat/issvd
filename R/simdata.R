#' Generate mupltiple-sources simulation datasets.
#'
#' Generate mupltiple-sources simulation datasets following the settings.
#'
#' @param n The number of x.
#' @param ncomp The number of factors in the simulated data.
#' @param p p1 p2 p3 represent the number of variables in three datasets, respectively. 
#' @param sig sig1 sig2 sig3 represent the noise ratio in three datasets, respectively. 
#' @param d d1 d2 d3 represent the scale parameters in three datasets, respectively.
#' @return \item{x}{multi-source data.}
#'         \item{truex}{multi-source data without noise.}

simdata=function(n,p1,p2,p3,sig1,sig2,sig3,d1,d2,d3){

  p=p1+p2+p3
  u=matrix(0,n,1)

  u = matrix(c(10,9,8,7,6,5,4,3,rep(2,17),rep(0,n-25)))
  u=apply(u,2,function(x) x/norm(x,'2'))
  v1=matrix(c(10,-10,8,-8,5,-5, rep(3,5),rep(-3,5),rep(0,p1-16)))
  v2=matrix(c(10,-10,8,-8,5,-5, rep(3,5),rep(-3,5),rep(0,p2-16)))
  v3=matrix(c(10,-10,8,-8,5,-5, rep(3,5),rep(-3,5),rep(0,p3-16)))
  v1=apply(v1,2,function(x) x/norm(x,'2'))
  v2=apply(v2,2,function(x) x/norm(x,'2'))
  v3=apply(v3,2,function(x) x/norm(x,'2'))
  
  d=d1
  truex1=u%*%d%*%t(v1)
  x1=truex1+sig1*max(truex1)*matrix(rnorm(n*p1),n,p1)
  
  d=d2
  truex2=u%*%d%*%t(v2)
  x2=truex2+sig2*max(truex2)*matrix(rnorm(n*p2),n,p2)
  
  d=d3
  truex3=u%*%d%*%t(v3)
  x3=truex3+sig3*max(truex3)*matrix(rnorm(n*p3),n,p3)
  
  x=list(x1,x2,x3)
  truex=list(truex1,truex2,truex3)
  sim_load=list(v1,v2,v3)
  index=list(seq(17,p1),seq(17,p2),seq(17,p3))
  uindex=seq(26,n)
  
  list(x=x,truex=truex)
}

