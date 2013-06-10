###################################################################################
#' Constructor of LarsPath class
#'
#' This class stores the results of lars and fusion algorithms.
#'
#' \describe{
#'   \item{nbStep}{Number of steps of the algorithm.}
#'   \item{variable}{List of vector of size "step+1". The i+1-th item contains the index of non-zero coefficients at the i-th step.}
#'   \item{coefficient}{List of vector of size "step+1". The i+1-th item contains the non-zero coefficients at the i-th step.}
#'   \item{lambda}{Vector of length "step+1", containing the L1-norm of the coefficients at each step.}
#'   \item{dropIndex}{Vector of length "step" containing the index of the dropped variable at the i-th step, 0 means no variable has been dropped at this step.}
#'   \item{addIndex}{Vector of length "step" containing the index of the added variable at the i-th step, 0 means no variable has been added at this step.}
#'	 \item{mu}{Intercept.}
#'	 \item{ignored}{A vector containing index of ignored variable during the algorithm.}
#'   \item{p}{Total number of covariates.}
#'	 \item{fusion}{If TRUE,  results from MPA.fusion function.}
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
    lambda="numeric",
    dropIndex="numeric",
    addIndex="numeric",
    nbStep="numeric",
    mu="numeric",
    ignored="numeric",
    fusion="logical",
    p="numeric"
    ),
  prototype=prototype(
    variable=list(),
    coefficient=list(),
    lambda=numeric(0),
    dropIndex=numeric(0),
    addIndex=numeric(0),
    nbStep=numeric(0),
    mu=numeric(0),
    ignored=numeric(0),
    fusion=FALSE,
    p=numeric(0)
    )
  )


###################################################################################
#' 
#'  plot the path of the lars algorithm.
#'  
#' \describe{
#'   \item{x}{LarsPath object.}
#'   \item{sep.line}{If TRUE, print vertical dashed line when a variable is added or dropped in the path.}
#'   \item{...}{Other plot arguments.}
#' }
#' 
#' @title plot methods for LarsPath object
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
      plot(NA,xlim=c(min(x@lambda),max(x@lambda)),ylim=c(miny,maxy),main="Path",xlab="lambda",ylab="coefficients") 
      abline(h=0)
      lines(x@lambda[1:2],c(0,x@coefficient[[2]][1]),col=which(var==x@variable[[2]][1]))
      
      for(i in 2:(length(x@lambda)-1))
      {
        if(sep.line)
          abline(v=x@lambda[i],col="blue",lty=2)
        
        if(x@dropIndex[i]==0)##plot add case 
        {
          for(j in 1:length(x@coefficient[[i]]))
            lines(x@lambda[i:(i+1)],c(x@coefficient[[i]][j],x@coefficient[[i+1]][j]),col=which(var==x@variable[[i]][j]))
          
          if(x@addIndex[i]!=0)
            lines(x@lambda[i:(i+1)],c(0,x@coefficient[[i+1]][j+1]),col=which(var==x@variable[[i+1]][j+1]))
        }
        else
        {
          j=1  

          while(x@variable[[i]][j]!=x@dropIndex[i])
          {
            lines(x@lambda[i:(i+1)],c(x@coefficient[[i]][j],x@coefficient[[i+1]][j]),col=which(var==x@variable[[i]][j]))
            j=j+1
          }
          lines(x@lambda[i:(i+1)],c(x@coefficient[[i]][j],0),col=which(var==x@variable[[i]][j]))
          j=j+1
          while(j<=length(x@coefficient[[i]]))
          {
            lines(x@lambda[i:(i+1)],c(x@coefficient[[i]][j],x@coefficient[[i+1]][j-1]),col=which(var==x@variable[[i]][j]))
            j=j+1
          }
          if(x@addIndex[i]!=0)
            lines(x@lambda[i:(i+1)],c(0,x@coefficient[[i+1]][length(x@coefficient[[i]])]),col=which(var==x@variable[[i+1]][length(x@coefficient[[i]])]))        
        }
      }
      
      if(sep.line)
        abline(v=x@lambda[length(x@lambda)],col="blue",lty=2)
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
#' dataset=MPA.simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=MPA.fusion(dataset$data,dataset$response) 
#' plot.coefficient(result,45) #plot coefficients at the step 45
#' @export 
#' 
plot.coefficient=function(x,step,ylab="coefficients",xlab="variables",...)
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
    plot(x@variable[[step+1]],x@coefficient[[step+1]],xlab=xlab,ylab=ylab,...)
    points((1:x@p)[-x@variable[[step+1]]],rep(0,x@p-length(x@coefficient[[step+1]])),...)
  }
}



#' Get the vector of coefficients at a given step
#'
#' @title get coefficients at a given step.
#' @param x A LarsPath object.
#' @param step The step at which you want to get the coefficients.
#' @return a vector of size p containing the value of coefficients at the desired step.
#' @examples 
#' dataset=MPA.simul(50,10000,0.4,10,50,matrix(c(0.1,0.8,0.02,0.02),nrow=2))
#' result=MPA.fusion(dataset$data,dataset$response)
#' MPA.coef(result,45) #get the coefficients
#' @export 
#'
MPA.coef=function(x,step)
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
