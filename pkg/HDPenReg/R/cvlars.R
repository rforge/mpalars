#' cross validation function for lars algorithm
#'
#' @title cross validation
#' @author Quentin Grimonprez
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param nbFolds the number of folds for the cross-validation.
#' @param index Values at which prediction error should be computed. This is the fraction of the saturated |beta|. The default value is seq(0,1,by=0.01).
#' @param maxSteps Maximal number of steps for lars algorithm.
#' @param eps Tolerance of the algorithm.
#' @return A list containing 
#' \describe{
#'   \item{cv}{Mean prediction error for each value of index.}
#'   \item{cvError}{Standard error of cv.}
#'   \item{minCv}{Minimal cv criterion.}
#'   \item{fraction}{Value of the l1norm fraction for which the cv criterion is minimal.}
#'   \item{index}{Values at which prediction error should be computed. This is the fraction of the saturated |beta|. The default value is seq(0,1,by=0.01).}
#'   \item{maxSteps}{Maximum number of steps of the lars algorithm.}
#' }
#' @examples 
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDcvlars(dataset$data,dataset$response,5)
#' @export
#' 
HDcvlars <- function(X,y,nbFolds=10,index=seq(0,1,by=0.01),maxSteps=3*min(dim(X)),eps=.Machine$double.eps^0.5)
{
	#check arguments
	if(missing(X))
		stop("X is missing.")
	if(missing(y))
		stop("y is missing.")
	index=unique(index)
	.checkcvlars(X,y,maxSteps,eps,nbFolds,index)

	# call lars algorithm
	val=.Call( "cvlars",X,y,nrow(X),ncol(X),maxSteps,eps,nbFolds,index,PACKAGE = "HDPenReg" )
	
	#create the output object
	cv=list(cv=val$cv,cvError=val$cvError,minCv=min(val$cv),fraction=index[which.min(val$cv)],index=index,maxSteps=maxSteps)

	class(cv)="cvlars"
	#plot.cv(cv)

	return(cv)
}


# check arguments from cvlars 
.checkcvlars=function(X,y,maxSteps,eps,nbFolds,index)
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

	## nbFolds
	if(!.is.wholenumber(nbFolds))
		stop("nbFolds must be a positive integer")
	if(nbFolds<=0 || nbFolds>length(y))
		stop("nbFolds must be a positive integer")

	## eps
	if(!is.double(eps))
		stop("eps must be a positive real")	
	if(eps<=0)
		stop("eps must be a positive real")	

	## index
	if(!is.numeric(index) || !is.vector(index))
		stop("index must be a vector of real between 0 and 1")
	if(max(index)>1 || min(index)<0)
		stop("index must be a vector of real between 0 and 1")
}

#' plot cross validation mean square error
#'
#' @title plot cross validation mean square error
#' @author Quentin Grimonprez
#' @param x Output from HDcvlars function.
#' @examples 
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDcvlars(dataset$data,dataset$response,5)
#' plot.cv(result)
#' @export
#' 
plot.cv=function(x)
{
	if(missing(x))
		stop("x is missing.")
	if(class(x)!="cvlars")
		stop("x must be an output of the HDcvlars function.")
	plot(x$index, x$cv, type = "b", ylim = range(x$cv, x$cv + x$cvError, x$cv - x$cvError),xlab="Fraction L1 Norm",ylab="Cross-Validated MSE")
	lines(x$index, x$cv+x$cvError,lty=2)
	lines(x$index, x$cv-x$cvError,lty=2)
}
