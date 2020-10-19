#' Performs cross-validation method to select tuning parameters in the ISSVD.
#' @param x list with K Data matrices of dimension $n x pk$.
#'          We are interested in finding sparse singular value decomposition of
#'          dimension $pk$.
#' @param ncomp The number of factors; default is 1.
#' @param varnumv list with optional user-supplied varnumv sequence.
#' @param varnumu list with optional user-supplied varnumu sequence.
#' @param type Which kind of penalty is chosen? Possible options are \code{"soft"}, \code{"hard"}, \code{"bic"},
#'             and \code{"SCAD"}. Default is \code{"soft"}.
#' @param nfolds The number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets.
#' Smallest value allowable is nfolds=3.
#' @return \item{varv_opt}{Optimal parameters varnumv for v.}
#'         \item{varu_opt}{Optimal parameter varnumu for u.}
#' @examples
#'    n=100;p1=p2=p3=50;sig1=sig2=sig3=0.1;
#'    d1=d2=1;d3=50;
#'    simdata=simdata(n,p1,p2,p3,sig1,
#'	                 sig2,sig3,d1,d2,d3)
#'    varnumv=list(c(15,15,15),c(16,16,16),
#'               c(17,17,17))
#'    varnumu=list(24,25,26)
#'    issvd_model1=cv.issvd(x,1,varnumv,
#'                          varnumu,type='soft',10)
#'    varv_opt=issvd_model1$varv_opt
#'    varu_opt=issvd_model1$varu_opt
#'    issvd_model=issvd(x,1,varv_opt,varu_opt,type='soft')

cv.issvd=function(x,ncomp,varnumv,varnumu,
                  type='soft',nfolds){

  len_varv=length(varnumv)
  len_varu=length(varnumu)
  K=length(x)
  p=unlist(lapply(x,ncol))
  sparse_ratio=matrix(0,len_varv,len_varu)
  varu_opt=list(); varv_opt=list();
  for(comp in 1:ncomp){
   for(m in 1:len_varv){
    for(h in 1:len_varu){
	percentRemove <- 1/nfolds
	total_error=rep(0,nfolds)
      for(i in 1L:nfolds){
         ToRemove=list()
	 for(l in 1:K){
	   randmat <- matrix(stats::runif(nrow(x[[l]]) * ncol(x[[l]])),
	   	                    ncol = ncol(x[[l]]))
	   ToRemove[[l]] <- ((i - 1) * percentRemove < randmat) & (randmat < i * percentRemove)
	   DATArm <- x
	   DATArm[[l]][ToRemove[[l]]] <- NA
	   for(c in 1:ncol(x[[l]])){
             indexc <- !is.na(DATArm[[l]][, c])
             DATArm[[l]][, c][!indexc] <- mean(DATArm[[l]][, c][indexc]) #missing values are replaced by column means
        }
      }
        varv=list(varnumv[[m]])
        varu=list(varnumu[[h]])
        model=issvd(DATArm,1,varv,varu, type='soft',maxiter=100,eps=10^(-4))

        cv_error=rep(0,K)
	for(k in 1:K){
           a1=model$U%*%model$D[[k]]%*%t(model$V[[k]])
           cv_error[k]=sum((a1[ToRemove[[k]]]-x[[k]][ToRemove[[k]]])^2)/sqrt(sum(ToRemove[[k]]))
	}
          total_error[i]=sum(cv_error)
	}
	  sparse_ratio[m,h]=mean(ratio)
    }
   }
    varnumv_index=which.min(sparse_ratio)%%len_varv
    if(varnumv_index==0){
      varnumv_index=len_varv
    }
    varnumu_index=ceiling(which.min(sparse_ratio)/len_varv)
    varv_opt[[comp]]=varnumv[[varnumv_index]]
    varu_opt[[comp]]=varnumu[[varnumu_index]]
  }
  
  list(varv_opt=varv_opt,
       varu_opt=varu_opt)
}
