#' Predict response of a new sample Xnew at a given level of penalty
#'
#' @title Prediction of response
#' @author Quentin Grimonprez
#' @param x a LarsParth object
#' @param Xnew a matrix (of size n*x@@p) of covariates.
#' @param lambda If mode =\"norm\", lambda represents the l1-norm of the coefficients with which we want to predict. If mode=\"fraction\", lambda represents the ratio (l1-norm of the coefficientswith which we want to predict)/(l1-norm maximal of the LarsPath object).
#' @param mode "fraction" or "norm".
#' @return The predicted response
#' @examples 
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDlars(dataset$data[1:40,],dataset$response[1:40])
#' y=HDpredict(result,dataset$data[41:50,],0.3,"fraction")
#' @export
#' 
HDpredict=function(x,Xnew, lambda, mode="fraction")
{
	
	if(missing(x))
		stop("x is missing.")
	if(class(x)!="LarsPath")
		stop("x must be a LarsPath object.")

	if( !(mode%in%c("fraction","norm")) ) 
		stop("mode must be \"fraction\" or \"norm\".")
	if(!is.numeric(lambda))
		stop("lambda must be a positive real.")
	if(length(lambda)>1)
		stop("lambda must be a positive real.")
	if(lambda < 0)
		stop("lambda must be a positive real.")

	fraction = lambda
	if(mode == "norm")
		fraction = lambda/x@l1norm[x@nbStep+1]

 
    yPred=rep(x@mu,nrow(Xnew))

    ##fraction = 0 : all coefficients are equal to 0
    if(fraction == 0)
      return(yPred)

    ##fraction = 1 : coefficients of the last step
	if (fraction >= 1)
    {
		yPred=yPred + Xnew[,x@variable[[x@nbStep+1]]]%*%x@coefficient[[x@nbStep+1]]	
		return(yPred)
    }

    ##fraction >0 and <1
    coeff=computeCoefficients(x,fraction);

	yPred=yPred + Xnew[,coeff$variable,drop=FALSE]%*%coeff$coefficient

	return(yPred);
}


#' Compute coefficients at a given level of penalty
#'
#' @title Compute coefficients
#' @author Quentin Grimonprez
#' @param x a LarsParth object
#' @param lambda If mode ="norm", lambda represents the l1-norm of the coefficients with which we want to predict. If mode="fraction", lambda represents the ratio (l1-norm of the coefficientswith which we want to predict)/(l1-norm maximal of the LarsPath object).
#' @param mode "fraction" or "norm" or "lambda.
#' @return A list containing 
#' \describe{
#'   \item{variable}{Index of non-zeros coefficients.}
#'   \item{coefficient}{non-zeros coefficients.}
#' }
#' @examples 
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDlars(dataset$data[1:40,],dataset$response[1:40])
#' coeff=computeCoefficients(result,0.3,"fraction")
#' @export
#' 
computeCoefficients = function(x,lambda,mode="fraction")
{
	if(missing(x))
		stop("x is missing.")
	if(class(x)!="LarsPath")
		stop("x must be a LarsPath object.")
	if( !(mode%in%c("fraction","norm","lambda"))  ) 
		stop("mode must be \"fraction\" or \"norm\".")
	if(!is.numeric(lambda))
		stop("lambda must be a positive real.")
	if(length(lambda)>1)
		stop("lambda must be a positive real.")
	if(mode == "fraction")
		lambda = lambda * x@l1norm[x@nbStep+1]
	abscissa = x@l1norm
	if(mode == "lambda")
	   abscissa = x@lambda
	if(lambda < 0)
		stop("lambda must be a positive real.")
	if(lambda >= x@l1norm[x@nbStep+1])
		return(list(variable=x@variable[[x@nbStep+1]],coefficient=x@coefficient[[x@nbStep+1]]))
	if(lambda == 0)
		return(list(variable=c(),coefficient=c()))
	

	index = 1;
	if(mode=="lambda")
	{
	   while( abscissa[index] > lambda )
        index=index+1;
	}
	else
	{
	   while( abscissa[index] < lambda )
        index=index+1;
	}

	index=index-1
	if(x@dropIndex[index]==0)#no drop variable
	{
		#add variable case
		if(x@addIndex[index]!=0)
		{
			coeff=c(
			.computeOrdinate(abscissa[index], abscissa[index+1], lambda, x@coefficient[[index]], x@coefficient[[index+1]][-length(x@coefficient[[index+1]])]), 
			.computeOrdinate(abscissa[index], abscissa[index+1], lambda, 0, x@coefficient[[index+1]][length(x@coefficient[[index+1]])]))
			variable=x@variable[[index+1]]

			return(list(variable=variable,coefficient=coeff))
		}
		else
		{
			coeff=.computeOrdinate(abscissa[index], abscissa[index+1], lambda, x@coefficient[[index]], x@coefficient[[index+1]])
			variable=x@variable[[index+1]]

			return(list(variable=variable,coefficient=coeff))
		}
		
	}
	else
	{
		#delete variable case
		i = 1;
		delIndex = which(x@variable[[index]] == x@dropIndex[index])

		##drop with an add variable
		if(x@addIndex[index]!=0)
		{
			coeff=c(
			.computeOrdinate(abscissa[index], abscissa[index+1], lambda, x@coefficient[[index]][1:(delIndex-1)], x@coefficient[[index+1]][1:(delIndex-1)]),
			.computeOrdinate(abscissa[index], abscissa[index+1], lambda, x@coefficient[[index]][delIndex], 0),
			.computeOrdinate(abscissa[index], abscissa[index+1], lambda, x@coefficient[[index]][-c(1:delIndex)], x@coefficient[[index+1]][delIndex:(length(x@coefficient[[index]])-1)]),
			.computeOrdinate(abscissa[index], abscissa[index+1], lambda, 0, x@coefficient[[index+1]][length(x@coefficient[[index+1]])]))
			variable=c(x@variable[[index]], x@variable[[index+1]][length(x@variable[[index+1]])])

			return(list(variable=variable,coefficient=coeff))
		}
		else
		{	##no add
			coeff=c(
			.computeOrdinate(abscissa[index], abscissa[index+1], lambda, x@coefficient[[index]][1:(delIndex-1)], x@coefficient[[index+1]][1:(delIndex-1)]),
			.computeOrdinate(abscissa[index], abscissa[index+1], lambda, x@coefficient[[index]][delIndex], 0),
			.computeOrdinate(abscissa[index], abscissa[index+1], lambda, x@coefficient[[index]][-c(1:delIndex)], x@coefficient[[index+1]][delIndex:length(x@coefficient[[index]])]))

			variable=x@variable[[index]]
			
			return(list(variable=variable,coefficient=coeff))
		}

	}
	
}

.computeOrdinate=function(x1,x2,x3,y1,y2)
{
  return(y1 + (y2-y1) * ((x3-x1)/(x2-x1)) )
}

#beta=rep(0,x@p)
#  if(x@fusion)
#  {
#    a=0
#    index=sort(variable,index.return=TRUE)$ix
#    for(i in 1:(length(variable)-1))
#    {
#      a=a+coefficient[index[i]]
#      beta[variable[index[i]]:(variable[index[i+1]]-1)]=rep(a, variable[index[i+1]]-variable[index[i]])
#      
#    }
#    a=a+coefficient[index[length(index)]]
#    beta[variable[index[length(index)]]:x@p]=rep(a,x@p-variable[index[length(index)]]+1)
#  }
