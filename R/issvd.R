#' Perform intergative and sparse singular value decomposition analysis
#'
#' Performs intergative and sparse singular value decomposition analysis 
#' to multiple data matrix with penalty on the columns and the rows.
#'
#' @param x list with K Data matrices of dimension $n x pk$, which can contain NA for missing
#' values. We are interested in finding sparse singular value decomposition of
#' dimension $pk$.
#' @param ncomp The number of factors in the ISSVD to be returned; default is 1.
#' @param varnumv vectors with K elements. How sparse do you want vk to be in the K data matrices?
#' Each elements of varnumv represents how many variables is left, hence its value must be between 1 and pk.
#' The smaller it is, the sparser v will be.
#' @param varnumu How sparse do you want u to be in the K data matrices?
#' @param type Which kind of penalty is chosen? Possible options are \code{"soft"}, \code{"hard"}, \code{"bic"},
#' and \code{"SCAD"}. Default is \code{"soft"}.
#' @param maxiter How many iterations should be performed. It is best to run at
#' least 100 of so. Default is 100.
#' @param eps convergence threshhold. The default threshhold is 10^(-4).
#' @return \item{U}{The common left singular matrix U.}
#'         \item{D}{The individual scale matrices D.}
#'         \item{V}{The individual right singular matrices V.}
#' @examples
#'    n=100;p1=p2=p3=50;sig1=sig2=sig3=0.1;
#'    d1=d2=1;d3=50;
#'    simdata=simdata(n,p1,p2,p3,sig1,
#'	                  sig2,sig3,d1,d2,d3)
#'    varv_opt=list(c(16,16,16))
#'    varu_opt=25
#'    issvd_model=issvd(x,1,varv_opt,varu_opt,type='soft')

issvd=function(x,ncomp=1,varnumv,varnumu,type="soft",maxiter=100,eps=10^(-4)){

  K=length(x);n=nrow(x[[1]])
  Tx=c(); p<-rep(0,K);
  V=list();D=list();
  for(k in 1:K){
    Tx=cbind(Tx,x[[k]])
    p[k]<-ncol(x[[k]])
    V[[k]]=matrix(0,p[k],ncomp);
    D[[k]]=matrix(0,ncomp,ncomp);
  }
  U=matrix(0,n,ncomp);
  
  for(i in 1:ncomp){
    if(varnumu[[i]]>n||varnumu[[i]]<1) stop("The value of varnumu must be between 1 and n")
    if(length(varnumv[[i]])!=K) stop("The length of varnumv must be equal with the number of x ")
    for(k in 1:K){
	  if(varnumv[[i]][k]>p[k]||varnumv[[i]][k]<1) stop("Each elements of varnumv must be between 1 and pk")
	}
  }

  for(i in 1:ncomp){
    # initialization step
    U.old=svd(Tx,ncomp,ncomp)$u[,1]
    V.old=list()
    for(k in 1:K){
     V.old[[k]]=svd(x[[k]],ncomp,ncomp)$v[,1]
    }
    V.cur=V.old
    D.old=list()
    for(k in 1:K){
      D.old[[k]]=svd(x[[k]],ncomp,ncomp)$d[1]
    }
    D.cur=D.old

    # iteration step
    u.d <- v.d <- d.d<-1;iter=1;
    while((u.d>eps|v.d>eps|d.d>eps)&&iter<maxiter){
      # Updating v
      for(k in 1:K){
        lambda <- c()
        D_old_k_i=as.numeric(D.old[[k]])
        f=as.matrix(U.old)%*%D_old_k_i;
        lambda[k] <- sort(abs(t(x[[k]])%*%f))[p[k]-varnumv[[i]][k]]
         V.cur[[k]]=t(solve(D_old_k_i%*%D_old_k_i)%*%
                      thresh(t(f)%*%x[[k]],type,lambda[k], a=3.7));
         V.cur[[k]] <- V.cur[[k]]/norm(V.cur[[k]],'2')
      }

      TF=matrix(0,n,1)
      for(k in 1: K){
        TF=TF+x[[k]]%*%V.cur[[k]]%*%t(D.old[[k]])
      }

      TV=0
      for(k in 1:K){
        TV=TV+D.old[[k]]%*%t(D.old[[k]])
      }

      # Updating u
      lambda <- sort(abs(TF))[n-varnumu[[i]]]
      TF=thresh(TF,type,lambda,a=3.7)
      U.cur=TF%*%solve(TV)
      U.cur <-U.cur/norm(U.cur,'2')

      # Updating d
      for(k in 1:K){
        D.cur[[k]]=t(U.cur)%*%x[[k]]%*%V.cur[[k]]
      }

      u.d <- sqrt(sum((U.cur-U.old)^2))
      v.d <- sqrt(sum((unlist(V.cur)-unlist(V.old))^2))
      d.d <- sqrt(sum((unlist(D.cur)-unlist(D.old))^2))
      D.old=D.cur;V.old=V.cur;U.old=U.cur
      iter=iter+1
    }

    U[,i]=U.cur[,1]
    Tx=c();
    for(k in 1:K){
      V[[k]][,i]=V.cur[[k]]
      D[[k]][i,i]=D.cur[[k]]
      x[[k]]=x[[k]]-U[,i]%*%D.cur[[k]]%*%t(V.cur[[k]])
      Tx=cbind(Tx,x[[k]])
    }
  }

  return(list(U=U,V=V,D=D))
}

