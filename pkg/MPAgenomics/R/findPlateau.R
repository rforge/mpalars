#'
#' From the results of a segmentation of a signal for different values of a segmentation parameter lambda, this function will
#' search an optimal value of lambda corresponding to the biggest plateau (stabilization in the number of breakpoints).
#'
#' @title Find the best choice of segmentation parameter.
#' @param resSeg a list, each element of the list is a vector with the breakpoints for a value of Lambda.
#' @param Lambda vector with the values of Lambda.
#' @param plot if TRUE, some graphics will be plotted.
#' @param verbose if TRUE print some informations.
#'
#' @return a list containing:
#' \describe{
#'   \item{lambda}{Optimal parameter found.}
#'   \item{maxPlateau}{A vector with the first and the last position of the biggest plateau.}
#'   \item{plateau}{A matrix of 3 columns, each row corresponds to a different plateau. The first colum is the starting value of a plateau,
#'    the second, the length of the plateau and the third, the number of values of lambda contained in the plateau.}
#' }
#' 
#' @author Quentin Grimonprez
#' 
#' @export
#' 
findPlateau=function(resSeg,Lambda,plot=TRUE,verbose=TRUE)
{
  seg=c(0,ncol=2)
    
  #number of segment per lambda
  nbSeg=sapply(resSeg,length)

  if(plot)
  {
    plot(Lambda[1:length(nbSeg)],nbSeg,ylim=c(min(nbSeg)-1,min(max(nbSeg),min(nbSeg)+50)),type="l",xlab="lambda",ylab="Number of segments")
    tab=table(unlist(resSeg))
    plot(tab/max(tab),ylab="Frequency",xlab="Probes")
  }
  
  #difference beween the number of segment between 2 consecutive values of Lambda
  chute=diff(nbSeg)
    
  #element included in a plateau
  elementOfPlateau=which(chute==0)
      
  if(length(elementOfPlateau)>1)#THERE IS AT LEAST 1 PLATEAU
  {
    #seg will contains all the plateaux
    seg=c(elementOfPlateau[1],0)
      
    #search all the plateaux
    for(i in 2:length(elementOfPlateau))
    {
      #if the difference between 2 index of elementOfPlateau is different from 1, the second index belong to an other plateau
      if( (elementOfPlateau[i]-elementOfPlateau[i-1])!=1 )
      {
        # the second index is the start of a new plateau
        seg=rbind(seg,c(elementOfPlateau[i],0))
        #the first index +1 is the end of the previous plateau
        seg[nrow(seg)-1,2]=elementOfPlateau[i-1]+1
      }
    }
      
    #complete the index of the end of the last plateau
    if(length(seg)>2)
    {
       seg[nrow(seg),2]=elementOfPlateau[length(elementOfPlateau)]+1
    }
    else
    {
      seg[2]=elementOfPlateau[length(elementOfPlateau)]+1  
      seg=matrix(seg,nrow=1)
    }    
      
    #find the biggest plateau
    A=matrix(Lambda[seg],ncol=2)  
    maxPlateau=A[which.max(A[,2]-A[,1]),]
        
  }
  else if(length(elementOfPlateau==1))#1 PLATEAU
  {
    maxPlateau=Lambda[elementOfPlateau:(elementOfPlateau+1)]
    seg=matrix(elementOfPlateau:(elementOfPlateau+1),nrow=1)
  }
  else #NO PLATEAU
  {
    maxPlateau=c(NA,NA)
    seg=matrix(maxPlateau,nrow=1)
  } 
  
  #print the optimal lambda : the first lambda of the biggest plateau
  if(!is.na(maxPlateau[1]))
  {
    #formatting the results
    plateau=cbind(Lambda[seg[,1]],Lambda[seg[,2]]-Lambda[seg[,1]],seg[,2]-seg[,1]+1)  
    colnames(plateau)=c("start lambda","length","number of lambda")
    rownames(plateau)=rep(NULL,nrow(seg))
    
    #sort the plateaux by their size
    ind=sort(plateau[,2],decreasing=TRUE,index.return=TRUE)$ix
    plateau=plateau[ind,]
    
    lambda=maxPlateau[1]
    if(verbose)
      cat("optimal parameter: ",maxPlateau[1],"\n")
  }
  else
  {
    plateau=matrix(rep(NA,3),nrow=1)
    lambda=Lambda[length(resSeg)]
    if(verbose)
      cat("optimal parameter: ",lambda,"\n")
  }
  
  res=list(lambda=lambda,maxPlateau=maxPlateau,plateau=plateau)
  invisible(res)
}
