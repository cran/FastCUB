#' @title Hadamard product of a matrix with a vector
#' @description Return the Hadamard product between the given matrix and vector: this operation corresponds 
#' to multiply every row of the matrix by the corresponding element of the vector, and it is equivalent to the 
#' standard matrix multiplication to the right with the diagonal matrix whose diagonal is the given vector. 
#' It is possible only if the length of the vector equals the number of rows of the matrix, otherwise it prints 
#'  an error message. 
#' @aliases Hadprod
#' @usage Hadprod(Amat, xvett)
#' @param  Amat  A generic matrix
#' @param xvett A generic vector
#' @return A matrix of the same dimensions as Amat
#' @keywords internal 
#' @details It is an auxiliary function needed for computing the variance-covariance matrix of the estimated model 
#' with covariates.


Hadprod<-function(Amat,xvett){
  ra=nrow(Amat);ca=ncol(Amat); 
  dprod=matrix(0,ncol=ca,nrow=ra);
  for(i in 1:ra){dprod[i,]=Amat[i,1:ca]*xvett[i]};
  return(dprod);
}