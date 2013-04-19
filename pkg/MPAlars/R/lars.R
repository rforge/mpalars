#' lars algorithm
#'
#' @title Lars algorithm
#' @author Quentin Grimonprez
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param maxSteps Maximal number of steps for lars algorithm.
#' @param eps Tolerance of the algorithm.
#' @param verbose If TRUE print details during the process.
#' @return An object of type LarsPath.
#' @examples 
#' data=MPA.simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=MPA.lars(data$data,data$response)
#' @export
#' @useDynLib MPAlars
#' 
MPA.lars <- function(X,y,maxSteps=3*min(dim(X)),eps=.Machine$double.eps^0.5,verbose=FALSE)
{
	#check arguments
	if(missing(X))
		stop("X is missing.")
	if(missing(y))
		stop("y is missing.")
	.check(X,y,maxSteps,eps,verbose)

	# call lars algorithm
	val=.Call( "lars",X,y,nrow(X),ncol(X),maxSteps,eps,verbose ,PACKAGE = "MPAlars" )
	
	#create the output object
	path=new("LarsPath",variable=val$varIdx,coefficient=val$varCoeff,lambda=val$lambda,addIndex=val$evoAddIdx,dropIndex=val$evoDropIdx,
           nbStep=val$step,mu=val$mu,ignored=val$ignored,p=ncol(X))
	return(path)
}

#' fusion algorithm
#'
#' @title Fusion algorithm
#' @author Quentin Grimonprez
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param maxSteps Maximal number of steps for lars algorithm.
#' @param eps Tolerance of the algorithm.
#' @param verbose If TRUE print details during the process.
#' @return An object of type LarsPath.
#' @examples  
#' data=MPA.simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=MPA.fusion(data$data,data$response)
#' @export 
#' 
MPA.fusion <- function(X,y,maxSteps=3*min(dim(X)),eps=.Machine$double.eps^0.5,verbose=FALSE)
{
	#check arguments
	if(missing(X))
		stop("X is missing.")
	if(missing(y))
		stop("y is missing.")
	.check(X,y,maxSteps,eps,verbose)

	# call fusion algorithm
	val=.Call( "fusion",X,y,nrow(X),ncol(X),maxSteps,eps,verbose ,PACKAGE = "MPAlars" )
	
	#create the output object
	path=new("LarsPath",nbStep=val$step,variable=val$varIdx,coefficient=val$varCoeff,lambda=val$lambda,addIndex=val$evoAddIdx,dropIndex=val$evoDropIdx,p=ncol(X),fusion=TRUE)

	return(path)
}

# check arguments from lars and fusion algorithm
.check=function(X,y,maxSteps,eps,verbose)
{
	## X: matrix of real
	if(!is.numeric(X) || !is.matrix(X))
		stop("X must be a matrix of real")
	
	## y: vector of real
	if(!is.numeric(y) || !is.vector(y))
		stop("y must be a vector of real")
	if(length(y)!=nrow(X))
		stop("The number of rows of X doesn't match with the length of y")
		
	## maxSteps
	if(!.is.wholenumber(maxSteps))
		stop("maxSteps must be a positive integer")
	if(maxSteps<=0)
		stop("maxSteps must be a positive integer")

	## eps
	if(!is.double(eps))
		stop("eps must be a positive real")	
	if(eps<=0)
		stop("eps must be a positive real")	

	##verbose
	if(!is.logical(verbose))
		stop("verbose must be a logical")

}

#check if a number is an integer
.is.wholenumber=function(x, tol = .Machine$double.eps^0.5)  
{
	abs(x - round(x)) < tol
}
