#' cross validation function for EM lasso.
#'
#' @title cross validation for EM lasso
#' @author Quentin Grimonprez, Serge Iovleff
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param nbFolds the number of folds for the cross-validation.
#' @param lambda Values at which prediction error should be computed.
#' @param maxSteps Maximal number of steps for EM algorithm.
#' @param burn Number of steps for the burn period.
#' @param intercept If TRUE, there is an intercept in the model.
#' @param model "linear" or "logistic".
#' @param threshold Zero tolerance. Coefficients under this value are set to zero.
#' @param eps Tolerance of the EM algorithm.
#' @param epsCG Epsilon for the convergence of the conjugate gradient.
#' @return A list containing 
#' \describe{
#'   \item{cv}{Mean prediction error for each value of index.}
#'   \item{cvError}{Standard error of cv.}
#'   \item{minCv}{Minimal cv criterion.}
#'   \item{lambda}{Values of lambda at which prediction error should be computed.}
#'   \item{lambda.optimal}{Value of lambda for which the cv criterion is minimal.}
#' }
#' @examples 
#' dataset=simul(50,1000,0.4,1,10,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=EMcvlasso(dataset$data,dataset$response,5,lambda=5:1,intercept=FALSE)
#' @export
#' 
EMcvlasso <- function( X , y, lambda=NULL, nbFolds=10
                     , maxSteps=1000, intercept=TRUE, model="linear"
                     , burn=30, threshold=1.e-08, eps=1e-5, epsCG=1e-8)
{
  #check arguments
  if(missing(X)) stop("X is missing.")
  if(missing(y)) stop("y is missing.")
  if(is.null(lambda)) { lambda=-1;}
  else
  {
    lambda=unique(lambda)
    lambda=sort(lambda)
  }

  ## threshold
  if(!is.double(threshold)) stop("threshold must be a positive real") 
  if(threshold<=0) stop("threshold must be a positive real") 

  ## epsCG
  if(!is.double(epsCG)) stop("epsCG must be a positive real") 
  if(epsCG<=0) stop("epsCG must be a positive real") 
  # check cv      
  mode=c("fraction","lambda")
  .checkcvlars(X,y,maxSteps,eps,nbFolds,c(0,1),intercept,mode)
  ## maxSteps
  if(!.is.wholenumber(burn)) stop("burn must be a positive integer.")
  if( (burn<=0) || (burn>maxSteps) ) stop("burn must be a positive integer lesser than maxSteps.")
  #model
  if(!(model%in%c("linear","logistic"))) stop("The model must be \"linear\" or \"logistic\"")
  # call cv for lasso
  val=list()
  if(model=="linear")
       val=.Call( "cvEMlasso",X,y,lambda,nbFolds,intercept,maxSteps,burn,threshold,eps,epsCG,PACKAGE = "HDPenReg" )
  else
       val=.Call( "cvEMlogisticLasso",X,y,lambda,nbFolds,intercept,maxSteps,burn,threshold,eps,epsCG,PACKAGE = "HDPenReg" )
  return(val)
}


#' cross validation function for EM lasso.
#'
#' @title cross validation for EM fused-lasso
#' @author Quentin Grimonprez, Serge Iovleff
#' @param X the matrix (of size n*p) of the covariates.
#' @param y a vector of length n with the response.
#' @param nbFolds the number of folds for the cross-validation.
#' @param lambda1 Values of lambda1 at which prediction error should be computed. Can be a single value.
#' @param lambda2 Values of lambda2 at which prediction error should be computed. Can be a single value.
#' @param maxSteps Maximal number of steps for EM algorithm.
#' @param burn Number of steps for the burn period.
#' @param intercept If TRUE, there is an intercept in the model.
#' @param model "linear" or "logistic".
#' @param eps0 Zero tolerance. Coefficients under this value are set to zero.
#' @param eps Tolerance of the algorithm.
#' @param epsCG Epsilon for the convergence of the conjugate gradient.
#' @return A list containing 
#' \describe{
#'   \item{cv}{Mean prediction error for each value of index.}
#'   \item{cvError}{Standard error of cv.}
#'   \item{minCv}{Minimal cv criterion.}
#'   \item{lambda1}{Values of lambda1 at which prediction error should be computed.}
#'   \item{lambda2}{Values of lambda2 at which prediction error should be computed.}
#'   \item{lambda.optimal}{Value of (lambda1,lambda2) for which the cv criterion is minimal.}
#' }
#' @examples 
#' dataset=simul(50,1000,0.4,1,10,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=EMcvfusedlasso(dataset$data,dataset$response,5,lambda1=3:1,lambda2=3:1,intercept=FALSE)
#' @export
#' 
EMcvfusedlasso <- function( X, y, lambda1, lambda2
                          , nbFolds=10
                          , maxSteps=1000, burn=50, intercept=TRUE, model="linear"
                          , eps=1e-5, eps0=1e-8, epsCG=1e-8)
{
  #check arguments
  if(missing(X)) stop("X is missing.")
  if(missing(y)) stop("y is missing.")
  if(missing(lambda1)) stop("lambda1 is missing.")
  if(missing(lambda2)) stop("lambda2 is missing.")

  ## threshold
  if(!is.double(eps0)) stop("eps0 must be a positive real") 
  if(eps0<=0) stop("eps0 must be a positive real") 
  ## epsCG
  if(!is.double(epsCG)) stop("epsCG must be a positive real") 
  if(epsCG<=0) stop("epsCG must be a positive real") 

  .checkcvlars(X,y,maxSteps,eps,nbFolds,c(0,1),intercept,"fraction")
  
  ## maxSteps
  if(!.is.wholenumber(burn)) stop("burn must be a positive integer.")
  if( (burn<=0) || (burn>maxSteps) ) stop("burn must be a positive integer lesser than maxSteps.")
  #lambda1
  .check.lambda(lambda1)
  lambda1=sort(lambda1)
   #lambda
  .check.lambda(lambda2)
  lambda2=sort(lambda2)
  
  #model
  if(!(model%in%c("linear","logistic"))) stop("The model must be \"linear\" or \"logistic\"")
  
  val=list()
  if( (length(lambda1)==1) && (length(lambda2)==1) )
  {
    val$lambda1=lambda1
    val$lambda2=lambda2
    val$lambda.optimal=c(lambda1,lambda2)
  }
  else
  {
    if(length(lambda1)==1)
    {
      optimL1=FALSE
      if(model=="linear")
  		{
        val=.Call( "cvEMfusedLasso1D",X,y,lambda1,lambda2,optimL1,nbFolds,intercept,maxSteps,burn,eps0,eps,epsCG,PACKAGE = "HDPenReg" )
  		}
      else
         val=.Call( "cvEMlogisticFusedLasso1D",X,y,lambda1,lambda2,optimL1,nbFolds,intercept,maxSteps,burn,eps0,eps,epsCG,PACKAGE = "HDPenReg" )
      names(val)[1]="lambda2"
      val$lambda1=lambda1
    }
    else
    {
      if(length(lambda2)==1)
      {
        optimL1=TRUE
        if(model=="linear")
            val=.Call( "cvEMfusedLasso1D",X,y,lambda1,lambda2,optimL1,nbFolds,intercept,maxSteps,burn,eps0,eps,epsCG,PACKAGE = "HDPenReg" ) 
        else
            val=.Call( "cvEMlogisticFusedLasso1D",X,y,lambda1,lambda2,optimL1,nbFolds,intercept,maxSteps,burn,eps0,eps,epsCG,PACKAGE = "HDPenReg" ) 
        names(val)[1]="lambda1"
        val$lambda2=lambda2
      }
      else # 2D
      {
        if(model=="linear")
        	val=.Call( "cvEMfusedLasso2D",X,y,lambda1,lambda2,nbFolds,intercept,maxSteps,burn,eps0,eps,epsCG,PACKAGE = "HDPenReg" )
        else
          val=.Call( "cvEMLogisticFusedLasso2D",X,y,lambda1,lambda2,nbFolds,intercept,maxSteps,burn,eps0,eps,epsCG,PACKAGE = "HDPenReg" )
        val$lambda1=lambda1
        val$lambda2=lambda2            
      }
    }
  }
  return(val)
}
