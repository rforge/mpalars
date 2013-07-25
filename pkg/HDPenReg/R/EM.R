#' EM algorithm for lasso penalty
#'
#' @title EM algorithm for lasso penalty
#' @author Quentin Grimonprez
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param lambda a sequence of l1 penalty regularization term. If no sequence is provided, the function computes his own sequence.
#' @param maxSteps Maximal number of steps for lars algorithm.
#' @param eps Tolerance of the algorithm.
#' @return An object of type LarsPath.
#' @examples 
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=EMlasso(dataset$data,dataset$response)
#' @export
#' 
EMlasso <- function(X,y,lambda,maxSteps=3*max(dim(X)),eps=.Machine$double.eps^0.5)
{
	#check arguments
	if(missing(X))
		stop("X is missing.")
	if(missing(y))
		stop("y is missing.")
	.check(X,y,maxSteps,eps,TRUE)

	if(missing(lambda))
		lambda=-1#lambda will be generated in C code
	else
	{
		lambda=sort(lambda)
		.check.lambda(lambda)
		lambda=lambda[lambda>0]
	}


	# call lars algorithm
	val=.Call("EMlasso",X,y,nrow(X),ncol(X),lambda,maxSteps,eps,PACKAGE = "HDPenReg" )
	
	return(val)
}


# check lambda for EM
.check.lambda=function(lambda)
{
	## lambda: vector of real
	if(!is.numeric(lambda) || !is.vector(lambda))
		stop("lambda must be a vector of positive real.")
	if(length(which(lambda<0))>0)
		stop("lambda must contain only positive number.")
}
