###################################################################################
#' Constructor of LarsPath class
#'
#' This class stores the results of lars and fusion algorithms.
#'
#' \describe{
#'   \item{nbStep}{Number of steps of the algorithm.}
#'   \item{variable}{List of vector of size "step+1". The i+1-th item contains the index of non-zero coefficients at the i-th step.}
#'   \item{coefficient}{List of vector of size "step+1". The i+1-th item contains the non-zero coefficients at the i-th step.}
#'   \item{l1norm}{Vector of length "step+1", containing the L1-norm of the coefficients at each step.}
#'   \item{lambda}{Vector of length "step+1", containing the lambda at each step.}
#'   \item{dropIndex}{Vector of length "step" containing the index of the dropped variable at the i-th step, 0 means no variable has been dropped at this step.}
#'   \item{addIndex}{Vector of length "step" containing the index of the added variable at the i-th step, 0 means no variable has been added at this step.}
#'	 \item{mu}{Intercept.}
#'   \item{meanX}{Mean of columns of X.}
#'	 \item{ignored}{A vector containing index of ignored variables during the algorithm.}
#'   \item{p}{Total number of covariates.}
#'	 \item{fusion}{If TRUE,  results from HDfusion function.}
#'   \item{error}{Error message from lars.}
#' }
#'
#'
#' @name LarsPath-class
#' @rdname LarsPath-class
#' @exportClass LarsPath
#'
setClass(
  Class="LarsPath",
  representation=representation(
    variable="list",
    coefficient="list",
    l1norm="numeric",
    lambda="numeric",
    dropIndex="list",
    addIndex="list",
    nbStep="numeric",
    mu="numeric",
    meanX="numeric",
    ignored="numeric",
    fusion="logical",
    p="numeric",
    error="character"
    ),
  prototype=prototype(
    variable=list(),
    coefficient=list(),
    l1norm=numeric(0),
    lambda=numeric(0),
    dropIndex=list(),
    addIndex=list(),
    nbStep=numeric(0),
    mu=numeric(0),
    meanX=numeric(0),
    ignored=numeric(0),
    fusion=FALSE,
    p=numeric(0),
    error=character()
    )
  )


###################################################################################
#' 
#'  plot the path of the lars algorithm.
#'  
#' 
#' 
#' @title plot methods for LarsPath object
#' @param x LarsPath object
#' @param sep.line If TRUE, print vertical dashed line when a variable is added or dropped in the path
#' @param ... Other plot arguments
#' @docType methods
#' @rdname plot-methods
#' @name plot-methods 
#' @aliases plot,LarsPath-method plot-methods
#'
#' @export
setMethod(
  f="plot",
  signature="LarsPath",
  definition=function(x,sep.line=TRUE,...)
  {
    miny=0
    maxy=0
    for(i in 2:length(x@coefficient))
    {
      miny=min(miny,x@coefficient[[i]])
      maxy=max(maxy,x@coefficient[[i]])
    }
    var=unique(unlist(x@variable))
    plot(NA,xlim=c(min(x@l1norm),max(x@l1norm)),ylim=c(miny,maxy),main="Path",xlab="l1norm",ylab="coefficients") 
    abline(h=0)
    lines(x@l1norm[1:2],c(0,x@coefficient[[2]][1]),col=which(var==x@variable[[2]][1]))
    
    for(i in 2:(length(x@l1norm)-1))
    {
      if(sep.line)
        abline(v=x@l1norm[i],col="blue",lty=2)
      
      if(length(x@dropIndex[[i]])==0)##plot add case 
      {
        for(j in 1:length(x@coefficient[[i]]))
          lines(x@l1norm[i:(i+1)],c(x@coefficient[[i]][j],x@coefficient[[i+1]][j]),col=which(var==x@variable[[i]][j]))
        
        if(length(x@addIndex[[i]])!=0)
        {
          for(j in (length(x@coefficient[[i]])+1):(length(x@coefficient[[i]])+length(x@addIndex[[i]])) )
            lines(x@l1norm[i:(i+1)],c(0,x@coefficient[[i+1]][j]),col=which(var==x@variable[[i+1]][j]))
          
        }
      }
      else
      {
        j=1  
        for(k in 1:length(x@dropIndex[[i]]))
        {
          #plot line for variables before a dropped variable
          while(x@variable[[i]][j]!=x@dropIndex[[i]][k])
          {
            lines(x@l1norm[i:(i+1)],c(x@coefficient[[i]][j],x@coefficient[[i+1]][j]),col=which(var==x@variable[[i]][j]))
            j=j+1
          }
          
          #plot the line of the dropped variable
          lines(x@l1norm[i:(i+1)],c(x@coefficient[[i]][j],0),col=which(var==x@variable[[i]][j]))
          j=j+1
          
        }
        
        #plot the line of variables after every drop variable
        while(j<=length(x@coefficient[[i]]))
        {
          lines(x@l1norm[i:(i+1)],c(x@coefficient[[i]][j],x@coefficient[[i+1]][j-1]),col=which(var==x@variable[[i]][j]))
          j=j+1
        }
        
        #plot the line for added variables
        if(length(x@addIndex[[i]])!=0)
        {         
          for(j in length(x@coefficient[[i+1]]):(length(x@coefficient[[i+1]])-length(x@addIndex[[i]])+1) )
            lines(x@l1norm[i:(i+1)],c(0,x@coefficient[[i+1]][j]),col=which(var==x@variable[[i+1]][j]))        
        }
      }
    }
    
    if(sep.line)
      abline(v=x@l1norm[length(x@l1norm)],col="blue",lty=2)
    axis(4, at=x@coefficient[[length(x@coefficient)]],labels=x@variable[[length(x@variable)]])
  }
)



#' Plot of the coefficients of a step
#'
#' @title Plot of coefficients
#' @param x A LarsPath object.
#' @param step The step at which you want to plot the coefficients.
#' @param ylab Name of the y axis.
#' @param xlab Name of the x axis.
#' @param ... Other plot arguments.
#' @examples 
#' dataset=simul(50,1000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDfusion(dataset$data,dataset$response) 
#' plotCoefficient(result,result@@nbStep) #plot coefficients at the last step
#' @export 
#' 
plotCoefficient=function(x,step,ylab="coefficients",xlab="variables",...)
{
  if(missing(x))
    stop("x is missing.")
  if(missing(step))
    stop("step is missing.")
  if(class(x)!="LarsPath")
    stop("x must be a LarsPath object")
  if(!.is.wholenumber(step))
    stop("step must be a positive integer smaller than x@nbStep")
  if( (step<0) || (step>x@nbStep) )
    stop("step must be a positive integer smaller than x@nbStep")
    
  if(x@fusion)
  {
    index=sort(x@variable[[step+1]],index.return=TRUE)$ix
    if(x@variable[[step+1]][index[1]]!=1)
    { plot(c(1,x@variable[[step+1]][index[1]]-1),rep(0,2),xlim=c(1,x@p),ylim=c(min(0,cumsum(x@coefficient[[step+1]][index])),max(0,cumsum(x@coefficient[[step+1]][index]))),type="l",xlab=xlab,ylab=ylab)
    }else{
      plot(NA,xlim=c(1,x@p),ylim=c(min(cumsum(x@coefficient[[step+1]][index])),max(cumsum(x@coefficient[[step+1]][index]))),xlab=xlab,ylab=ylab)}
    
    a=0
    if(length(x@variable[[step+1]])>1)
      for(i in 1:(length(x@variable[[step+1]])-1))
      {
        a=a+x@coefficient[[step+1]][index[i]]
        lines(c(x@variable[[step+1]][index[i]],x@variable[[step+1]][index[i+1]]-1),rep(a,2))
      }
    
    a=a+x@coefficient[[step+1]][index[length(index)]]
    lines(c(x@variable[[step+1]][index[length(index)]],x@p),rep(a,2))
  }
  else
  {
    plot((1:x@p)[-x@variable[[step+1]]],rep(0,x@p-length(x@coefficient[[step+1]])),ylim=c(min(x@coefficient[[step+1]],0),max(x@coefficient[[step+1]])),xlim=c(1,x@p),xlab=xlab,ylab=ylab,...)
    points(x@variable[[step+1]],x@coefficient[[step+1]],...)    
  }
}



#' Get the vector of coefficients at a given step
#'
#' @title get coefficients at a given step.
#' @param x A LarsPath object.
#' @param step The step at which you want to get the coefficients.
#' @return a vector of size p containing the value of coefficients at the desired step.
#' @examples 
#' dataset=simul(50,1000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDfusion(dataset$data,dataset$response)
#' coefficient=coeff(result,result@@nbStep) #get the coefficients
#' @export 
#'
coeff=function(x,step)
{
  if(missing(x))
    stop("x is missing.")
  if(missing(step))
    stop("step is missing.")
  if(class(x)!="LarsPath")
    stop("x must be a LarsPath object")
  if(!.is.wholenumber(step))
    stop("step must be a positive integer smaller than x@step")
  if( (step<0) || (step>x@nbStep) )
    stop("step must be a positive integer smaller than x@step")
  
  beta=rep(0,x@p)
  
  if(x@fusion)
  {
    a=0
    index=sort(x@variable[[step+1]],index.return=TRUE)$ix
    for(i in 1:(length(x@variable[[step+1]])-1))
    {
      a=a+x@coefficient[[step+1]][index[i]]
      beta[x@variable[[step+1]][index[i]]:(x@variable[[step+1]][index[i+1]]-1)]=rep(a, x@variable[[step+1]][index[i+1]]-x@variable[[step+1]][index[i]])
      
    }
    a=a+x@coefficient[[step+1]][index[length(index)]]
    beta[x@variable[[step+1]][index[length(index)]]:x@p]=rep(a,x@p-x@variable[[step+1]][index[length(index)]]+1)
  }
  else
  {
    beta[x@variable[[step+1]]]=x@coefficient[[step+1]]
  }
  
  return(beta)
}


#' Compute coefficients at a given level of penalty
#'
#' @title Compute coefficients
#' @author Quentin Grimonprez
#' @param object a LarsParth object
#' @param index If mode ="norm", index represents the l1-norm of the coefficients with which we want to predict.
#'  If mode="fraction", index represents the ratio (l1-norm of the coefficientswith which we want to predict)/(l1-norm maximal of the LarsPath object).
#'  If mode="lambda", index represents the value of the penalty parameter. If mode="step", index represents the numer of the step at which we want coefficients.
#' @param mode "fraction" or "norm" or "lambda" or "step".
#' @param ... other arguments. Not used
#' @return A vector containing the estimated coefficient for index
#' @method coef LarsPath
#' @examples 
#' dataset=simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=HDlars(dataset$data[1:40,],dataset$response[1:40])
#' coeff=coef(result,0.3,"fraction")
#' @export
coef.LarsPath=function(object,index=NULL,mode=c("lambda","step","fraction","norm"),...)
{
  mode <- match.arg(mode)
  
  if(missing(object))
    stop("objectx is missing.")
  if(missing(index) || is.null(index))
    stop("index is missing.")
  if(class(object)!="LarsPath")
    stop("object must be a LarsPath object")
  
  beta=rep(0,object@p)
  
  if(mode=="step")
  {
    if(!.is.wholenumber(index))
      stop("index must be a positive integer smaller than object@step")
    if( (index<0) || (index>object@nbStep) )
      stop("index must be a positive integer smaller than object@step")
    
    if(object@fusion)
    {
      a=0
      ind=sort(object@variable[[index+1]],index.return=TRUE)$ix
      for(i in 1:(length(object@variable[[index+1]])-1))
      {
        a=a+object@coefficient[[index+1]][ind[i]]
        beta[object@variable[[index+1]][ind[i]]:(object@variable[[index+1]][ind[i+1]]-1)]=rep(a, object@variable[[index+1]][ind[i+1]]-object@variable[[index+1]][ind[i]])
        
      }
      a=a+object@coefficient[[index+1]][ind[length(ind)]]
      beta[object@variable[[index+1]][ind[length(ind)]]:object@p]=rep(a,object@p-object@variable[[index+1]][ind[length(ind)]]+1)
    }
    else
    {
      beta[object@variable[[index+1]]]=object@coefficient[[index+1]]
    }
  }
  else
  {
    betatemp=computeCoefficients(object,index,mode)
    beta[betatemp$variable]=betatemp$coefficient
  }

  return(beta)
}
