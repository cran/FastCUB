#' @title Recursive computation of the inverse of a matrix
#' @description Compute the variance-covariance matrix of the incomplete score vector involved in Louis' identity for the observed information matrix
#' @aliases invmatgen
#' @usage invmatgen(G,H,listE)
#' @param G Primary matrix for the sum decomposition of $G+H$
#' @param H Secondary matrix for the sum decomposition of $G+H$
#' @param listE  Auxiliary matrices that sum up to H
#' @return The inverse of matrix G + H computed recursively thanks to matrices listed in listE
#' @import stats
#' @references  Miller K. (1981).  On the inverse of the sum of matrices, 
#' \emph{Mathematics Magazine}, \bold{54}, 67--72 \cr
#' @seealso \code{\link{fastCUB}}
#' @keywords stats 

########################################

########################################

invmatgen<-function(G,H,listE){
  
  ### G deve essere diagonale
  
  r<-qr(H)$rank
  
  if (r != nrow(G)){
    cat("Procedure not applicable","\n")
    return(try(solve(G+H),silent=TRUE))
    
  } else {
    
    Clist<-Cinv<-list()
    Clist[[1]]<-G
    Cinv[[1]]<- diag(1/(diag(G)))
    nu<-c()
    for(l in 1:r){
      
      Clist[[l+1]]<- Clist[[l]] + listE[[l]]
      
      nu[l]<-1/(1+sum(diag(Cinv[[l]]%*%listE[[l]])))
      
      Cinv[[l+1]]<- Cinv[[l]] - nu[l]*(Cinv[[l]]%*%listE[[l]]%*%Cinv[[l]])
    }
    
    return(Cinv[[r+1]])
    
    
  }
  
  
  
}
