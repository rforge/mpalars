#'
#' This function selects, for each chromosome, the most relevant markers according to a response.
#'
#' @title markers selection
#' 
#' @param dataSetName The name of the data-set folder.
#' @param dataResponse A csv files or a data.frame with 2 columns : "files" and "response". The column "files" contains the filename to extract and the second column the response associated to the file.
#' @param chromosome A vector containing the number of the chromosomes for the SNPs selection.
#' @param signal either \"CN\" or \"fracB\". corresponding to which signal will be analyzed (default=\"CN\").
#' @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files.
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param onlySNP (only if signal=\"CN\"). If TRUE, only the SNPs probes are used (default=FALSE).
#' @param nbFolds number of folds in the cross validation (default=10).
#' @param loss either \"logistic\" (binary response) or \"linear\" (quantitative response), default is ”logistic”
#' @param plot If TRUE, cross-validation mean squared error is plotted (default=TRUE).
#' @param ... Other parameters for HDlars function or glmnet function.
#' 
#' @return a list containing length(chromosme) elements. Each element is a list containing
#' \describe{
#'   \item{chr}{chromosome corresponding to the signal.}
#'   \item{markers.index}{A vector containing the index of all selected markers.}
#'   \item{markers.position}{A vector containing the position of all selected markers.}
#'   \item{markers.names}{A vector containing the names of all selected markers.}
#'   \item{coefficient}{A vector containing the coefficients of all selected markers.}
#'   \item{intercept}{Intercept of the model.}
#   \item{fraction}{fraction of the l1 norm of coefficient selected by cross validation.}
#' }
#'
#' @details This function requires to use the aroma folder architecture. In your working directory, there must have the rawData folder and totalAndFracBData folder.
#' This function launches the lars algorithm on the CN or fracB data and uses a cross-validation to select the most appropriate solution.
#' 
#' @examples
#' #DO NOT EXECUTE
#' # res=markerSelection("DataB",rnorm(27),chromosome=22,signal="CN",normalTumorArray,onlySNP=TRUE)
#' 
#' @seealso HDPenReg
#'
#' @author Quentin Grimonprez
#'
#' @export
#'
markerSelection=function(dataSetName,dataResponse,chromosome=1:22,signal="CN",normalTumorArray,onlySNP=FALSE,nbFolds=min(length(dataResponse),10),loss="logistic",plot=TRUE,...)
{
  allpkg=TRUE
  if(!suppressPackageStartupMessages(require("aroma.affymetrix", quietly=TRUE) ) )
  {
    cat("Package not found: aroma.affymetrix. For download it:\n")
    cat("source(\"http://www.braju.com/R/hbLite.R\")\n")
    cat(" hbLite(\"sfit\")\n")
    cat("source(\"http://bioconductor.org/biocLite.R\")\n")
    cat("biocLite(\"affxparser\")\n")
    cat("source(\"http://aroma-project.org/hbLite.R\")\n")
    cat("hbInstall(\"aroma.affymetrix\")\n")
    allpkg=FALSE
  }
#   else
#     cat("Package aroma.affymetrix loaded.\n")
  
  if(!suppressPackageStartupMessages(require("aroma.cn", quietly=TRUE) ) )
  {
    cat("source(\"http://aroma-project.org/hbLite.R\")\n")
    cat("hbInstall(\"aroma.cn\")\n") 
    allpkg=FALSE
  }
#   else
#     cat("Package aroma.cn loaded.\n")
  
  if(!suppressPackageStartupMessages(require(HDPenReg,quietly=TRUE) ) )
  {
    allpkg=FALSE
    cat("The package HDPenReg is missing. You can install it with the following command:\n","install.packages(HDPenReg, repos=\"http://R-Forge.R-project.org\") \n")
  }
#   else
#     cat("Package HDPenReg loaded.\n")
 
  if(!suppressPackageStartupMessages(require(glmnet,quietly=TRUE) ) )
  {
    allpkg=FALSE
    cat("The package glmnet is missing. You can install it with the following command:\n","install.packages(glmnet) \n")
  }
#   else
#     cat("Package glmnet loaded.\n")
  
  if(!allpkg)
    stop("You have to install some packages : Follow the printed informations.")
  
  ################## check parameters
  #check if the user is in the right place
  if(!("totalAndFracBData"%in%list.files()))
    stop("There is no \"totalAndFracBData\", check if you are in the good working directory or if you have run the signalPreProcess function before.")
  
  #signal
  if(!(signal%in%c("CN","fracB","both")))
    stop("signal must be either \"CN\", \"fracB\" or \"both\".")
  #dataSetName
  if(!is.character(dataSetName))
    stop("dataSetName must be the name of a folder in rawData.")
  if(length(grep(dataSetName,list.files("totalAndFracBData")))==0) #check if data are available for this dataSetName
    stop("The dataSetName you specify doesn't exist in the totalAndFracBData, run the ... function at first.")
  
  #dataResponse
  if(is.character(dataResponse))
    dataResponse=read.csv(dataResponse)
  else
  {
    if(!is.data.frame(dataResponse))
      stop("dataResponse must be either the path to the normalTumorArray csv file or a data.frame containing the response.") 
  }
  if(!("files"%in%names(dataResponse)))
    stop("dataResponse does not contain the column \"files\".")
  if(!("response"%in%names(dataResponse)))
    stop("dataResponse does not contain the column \"names\".")
    
  #loss
  if(!(loss%in%c("logistic","linear")))
    stop("loss must be either \"logistic\" or \"linear\".")
  
  #chromosome : vector of integer between 1 and 25
  if( !is.numeric(chromosome) || !is.vector(chromosome) )
    stop("chromosome must be a vector containing the different chromosomes to analyze.")
  chromsome=unique(chromosome)  #keep only unique values
  if( sum(sapply(chromosome,FUN=function(chr){!is.wholenumber(chr)}))>0 )
    stop("chromosome must be a vector containing integers corresponding to the different chromosomes to analyze.") 
  if(min(chromosome)<1)
    stop("chromosome must contain integer between 1 and 25.")
  if(max(chromosome)>25)
    stop("chromosome must contain integer between 1 and 25.")
  
  ########################
  
  #launch the right function depending of the signal
  res=switch(signal,
         "CN"=SNPselectionCNsignal(dataSetName,dataResponse,chromosome,normalTumorArray,onlySNP,nbFolds,loss),
         "fracB"=SNPselectionFracBsignal(dataSetName,dataResponse,chromosome,normalTumorArray,nbFolds,loss),
         "both"=stop("Not yet implemented."))
  
  return(res)
}


#
# This function selects, for each chromosome, the most relevant SNPs.
#
# @title SNPs selection from CN signal
# 
# @param dataSetName Name of the dataset folder in rawData 
# @param dataResponse response associated to the data
# @param chromosome chromosome used in the study
# @param normalTumorArray only if you have normal and tumor profile in your data folder. A csv file or a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
# @param onlySNP (only if signal=\"CN\"). If TRUE, only the SNP markers will be used.
# @param nbFolds number of folds in the cross validation
# @param loss either \"logistic\" (binary response) or \"linear\" (quantitative response).
# @param plot if TRUE, plot the cross validation graphic
#
# @return a list containing length(chromosme) elements. Each element is a list containing
# \describe{
#   \item{chr}{chromosome corresponding to the signal.}
#   \item{variable}{A vector containing the index of all selected SNPs.}
#   \item{coefficient}{A vector containing the coefficeints of all selected SNPs.}
# }
#
# @details This function requires to use the aroma folder architecture. In your working directory, there lust have the rawData folder and totalAndFracBData folder.
#
# @author Quentin Grimonprez
#
SNPselectionCNsignal=function(dataSetName,dataResponse,chromosome,normalTumorArray,onlySNP,nbFolds=10,loss="logistic",plot=TRUE,...)
{
  res=list()
  for(chr in chromosome)
  {
    #get the matrix with the copy number signal for the chromosome chr
    C=getCopyNumberSignal(dataSetName,chr,normalTumorArray,onlySNP,listOfFiles=as.character(dataResponse$files),verbose=FALSE)
    C=C[[paste0("chr",chr)]]
    gc()
    #extract the response in the right order
    ind=sapply(names(C)[3:(ncol(C)-1)],FUN=function(x){match(x,dataResponse$files)})    
    if(sum(is.na(ind))!=0)
      stop(paste0("A response is missing for the following files : ",paste(C$sampleNames[is.na(ind)],collapse=", ")))      

    response=dataResponse$response[ind]
        
    if(loss=="linear")
    {
      #cross validation to choose the best l1 norm ratio
      rescv=HDcvlars(t(as.matrix(C[3:(ncol(C)-1)])), response, nbFolds,index = seq(0, 1, by = 0.01),...)
      
      if(plot)
      {
        plotCv(rescv)
        title(paste0("chr",chr))
      }
      
      #lars algorithm for obtaining all the path
      reslars=HDlars(t(as.matrix(C[3:(ncol(C)-1)])), response,...) #ajouter critere d'arret sur la norme, le nb de variable ??
      
      #we compute the coefficients for the value given by the HDcvlars function
      coeff=computeCoefficients(reslars,rescv$fraction,mode="fraction")
      
      
      intercept=reslars@mu
      if(length(coeff$variable)!=0)
      {
        index=order(coeff$variable)
        var=coeff$variable[index]
        pos=C$position[coeff$variable[index]]
        name=as.character(C$featureNames)[coeff$variable[index]]
        coef=coeff$coefficient[index]
      }
      else
      {
        pos=c()
        var=c()
        name=c()
        coef=c()
      }
      rm(reslars,coeff)
      gc()
    }
    else
    {
      if(loss=="logistic")
      {
        rescv=cv.glmnet(t(as.matrix(C[3:(ncol(C)-1)])), response, nfolds=nbFolds,family="binomial",...)
        if(plot)
        {
          plot(rescv)
          title(paste0("chr",chr))
        }
        coef=coef(rescv,s=rescv$lambda.min)
        
        ind=which(coeff!=0)
        
        var=ind[-1]-1
        pos=C$position[var]
        name=as.character(C$featureNames)[var]
        intercept=coef[1]
        coef=coef[var]
      }
    }
        
    res[[paste0("chr",chr)]]=list(chr=chr,markers.index=var,markers.position=pos,
                                  markers.names=name,coefficient=coef,
                                  intercept=intercept)
    
    #delete objects created during this loop
    rm(C,rescv,coef,pos,intercept,name,var)
    gc()
  }

  return(res)
}

#
# This function selects, for each chromosome, the most relevant SNPs.
#
# @title SNPs selection from fracB signal
# 
# @param dataSetName Name of the dataset folder in rawData 
# @param dataResponse response associated to the data
# @param chromosome chromosome used in the study
# @param normalTumorArray only if you have normal and tumor profile in your data folder. A csv file or a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
# @param onlySNP (only if signal=\"CN\"). If TRUE, only the SNP markers will be used.
# @param nFolds number of folds in the cross validation
# @param loss either \"logistic\" (binary response) or \"linear\" (quantitative response).
# @param if TRUE, plot the cross validation graphic
# 
# @return a list containing length(chromosme) elements. Each element is a list containing
# \describe{
#   \item{chr}{chromosome corresponding to the signal.}
#   \item{variable}{A vector containing the index of all selected SNPs.}
#   \item{coefficient}{A vector containing the coefficeints of all selected SNPs.}
# }
#
# @details This function requires to use the aroma folder architecture. In your working directory, there lust have the rawData folder and totalAndFracBData folder.
#
# @author Quentin Grimonprez
#
SNPselectionFracBsignal=function(dataSetName,dataResponse,chromosome,normalTumorArray,nbFolds=10,loss="logistic",plot=TRUE,...)
{
  res=list()
  for(chr in chromosome)
  {
    #besoin que de la tumeur, rajouter onlyTumor comme parametre dans getFracB
    #ajouter listOfFile quand on veux ou on ne peux pas travailler sur ttes les donnees?
    
    #get the fracB signal for the chromosome chr
    fracB=getFracBSignal(dataSetName,chr,normalTumorArray,listOfFiles=as.character(dataResponse$files),verbose=FALSE)
    

    fracB=fracB[[paste0("chr",chr)]]$tumor
    gc()
    
    #extract the response in the right order
    ind=sapply(names(fracB)[3:(ncol(fracB)-1)],FUN=function(x){match(x,dataResponse$files)})    
    if(sum(is.na(ind))!=0)
      stop(paste0("A response is missing for the following files : ",paste(C$sampleNames[is.na(ind)],collapse=", ")))      
    
    response=dataResponse$response[ind]
      
    if(loss=="linear")
    {
      #cross validation to choose the best l1 norm ratio
      rescv=HDcvlars(t(as.matrix(fracB[,3:(ncol(fracB)-1)])), response, nbFolds,index = seq(0, 1, by = 0.01),...)
      
      if(plot)
      {
        plotCv(rescv)
        title(paste0("chr",chr))
      }
      
      #lars algorithm for obtaining all the path
      reslars=HDlars(t(as.matrix(fracB[,3:(ncol(fracB)-1)])), response,...) #ajouter critere d'arret sur la norme, le nb de variable ??
      
      #we compute the coefficients for the value given by the HDcvlars function
      coeff=computeCoefficients(reslars,rescv$fraction,mode="fraction")      
      

      intercept=reslars@mu
      if(length(coeff$variable)!=0)
      {
        index=sort(coeff$variable,index.return=TRUE)$ix
        var=coeff$variable[index]
        pos=fracB$position[coeff$variable[index]]
        name=as.character(fracB$featureNames)[coeff$variable[index]]
        coef=coeff$coefficient[index]
      }
      else
      {
        pos=c()
        var=c()
        name=c()
        coef=c()
      }
      
      rm(reslars,coeff)
      gc()
    }
    else
    {
      if(loss=="logistic")
      {
        rescv=cv.glmnet(t(as.matrix(fracB[,3:(ncol(fracB)-1)])), response, nfolds=nbFolds,family="binomial",...)
        if(plot)
        {
          plot(rescv)
          title(paste0("chr",chr))
        }
        coef=coef(rescv,s=rescv$lambda.min)
        
        ind=which(coeff!=0)
        
        var=ind[-1]-1
        pos=fracB$position[var]
        name=as.character(fracB$featureNames)[var]
        intercept=coef[1]
        coef=coef[var]
      }
    }
    
    res[[paste0("chr",chr)]]=list(chr=chr,markers.index=var,markers.position=pos,
                                  markers.names=name,coefficient=coef,
                                  intercept=intercept)
    
    rm(fracB,rescv,coef,pos,intercept,name,var)
    gc()
  }
  
  return(res)
}


#'
#' This function selects the most relevant variables according to a response.
#'
#' @title SNPs selection
#' 
#' @param dataMatrix Matrix containing the data, each row is a different sample.
#' @param dataResponse response associated to the data.
#' @param nbFolds number of folds in the cross validation.
#' @param loss either \"logistic\" (binary response) or \"linear\" (quantitative response).
#' @param plot If TRUE plot cross-validation mean squared error (default=TRUE).
#' @param ... spplementary arguments for cv.glmnet function in case of logistic loss or for HDlars function for linear loss.
#' 
#' @return a list containing 
#' \describe{
#'   \item{variable}{A vector containing the index of all selected variables.}
#'   \item{coefficient}{A vector containing the coefficients of all selected variables.}
#'   \item{intercept}{Intercept of the model.}
#'}
#'
#' @author Quentin Grimonprez
#' 
#' @export
variableSelection=function(dataMatrix,dataResponse,nbFolds=min(length(dataResponse),10),loss="logistic",plot=TRUE,...)
{
  allpkg=TRUE
  if(!require(HDPenReg,quietly=TRUE))
  {
    allpkg=FALSE
    cat("The package changepoint is missing. You can install it with the following command:\n","install.packages(HDPenReg, repos=\"http://R-Forge.R-project.org\") \n")
  }
  
  if(!allpkg)
    stop("You have to install some packages : Follow the printed informations.")
  
  #check plot (other parameters will be checcked in HDcvlars function)
  if(!is.logical(plot))
    stop("plot must be a logical.")
  
  if(loss=="linear")
  {
    #cross validation to choose the best l1 norm ratio
    rescv=HDcvlars(dataMatrix, dataResponse, nbFolds,index = seq(0, 1, by = 0.01),...)
    
    if(plot)
    {
      plotCv(rescv)
    }
    
    #lars algorithm for obtaining all the path
    reslars=HDlars(dataMatrix, dataResponse,...)
    
    #we compute the coefficients for the value given by the HDcvlars function
    coeff=computeCoefficients(reslars,rescv$fraction,mode="fraction")
        
    var=coeff$variable
    coef=coeff$coefficient
    intercept=reslars@mu
    if(length(coeff$variable)!=0)
    {
      index=order(coeff$variable)
      var=coeff$variable[index]
      coef=coeff$coefficient[index]
    }

    rm(reslars,coeff)
    gc()
  }
  else
  {
    if(loss=="logistic")
    {
      rescv=cv.glmnet(dataMatrix, dataResponse, nfolds=nbFolds,family="binomial",...)
      if(plot)
      {
        plot(rescv)
      }
      coef=coef(rescv,s=rescv$lambda.min)
      
      ind=which(coeff!=0)
      
      var=ind[-1]-1
      pos=C$position[var]
      name=C$featureNames[var]
      intercept=coef[1]
      coef=coef[var]
    }
  }
  

  res=list(markers.index=var,coefficient=coef,intercept=intercept)  
}