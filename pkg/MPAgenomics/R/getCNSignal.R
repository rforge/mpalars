#' 
#' Extract copy-number signals from aroma files. It requires to have executed the normalization process of aroma 
#' package.
#'
#' @title Extract copy-number signal from aroma files
#' @param dataSetName name of the dataset (it must correpond to a folder name in rawData folder.).
#' @param chromosome A vector containing the chromosomes for which the CN signal must be extract. 
#' @param normalTumorArray only if you have normal and tumor profile in your data folder. A csv file or a data.frame with 2 columns: "normal" and "tumor".
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param onlySNP if TRUE only the copy-number at SNPs probes will be returned.
#' @param listOfFiles vector containing the names of the files in dataSetName to extract.
#' @param verbose if TRUE print some informations.
#' 
#' @return a list of length the number of chromosome containing a data.frame with columns:
#' \describe{
#'   \item{chromosome}{chromosome corresponding to the signal.}
#'   \item{position}{Positions associated to the allele B fraction.}
#'   \item{copynumber}{Many columns (the number of selected files), the name of each column is the name of the associated data file name. It contains the CN signal for a different profile.}
#'   \item{featureNames}{Names of the probes.}
#' }
#'  
#' @details You have to respect the aroma architecture. Your working directory must contain rawData folder and totalAndFracBData folder.
#' To access easily the names of the files available in your dataset, you can use the \link{getListOfFiles} function.
#' 
#' @examples 
#' #DO NOT EXECUTE
#' #C=getCopyNumberSignal("data1",5,normalTumorArray,TRUE)
#' #C=getCopyNumberSignal("data2",5,onlySNP=TRUE)
#'
#' @export 
getCopyNumberSignal=function(dataSetName,chromosome,normalTumorArray,onlySNP=FALSE,listOfFiles=NULL,verbose=TRUE)
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
  
  if(!suppressPackageStartupMessages(require("aroma.cn", quietly=TRUE) ) )
  {
    cat("source(\"http://aroma-project.org/hbLite.R\")\n")
    cat("hbInstall(\"aroma.cn\")\n") 
    allpkg=FALSE
  }
    
  if(!allpkg)
    stop("You have to install some packages : Follow the printed informations.")

  if(!("totalAndFracBData"%in%list.files()))
    stop("There is no \"totalAndFracBData\", check if you are in the good working directory or if you have run the signalPreProcess function before.")

  ################### check 
  #onlySNP
  if(!is.logical(onlySNP))
    stop("onlySNP must be either TRUE or FALSE.")
  
  #chromosome
  if(missing(chromosome))
    stop("chromosome is missing.")
  if(!is.numeric(chromosome))
    stop("chromosome must contain integers between 1 and 25.")
  if(sum(unlist(lapply(chromosome,is.wholenumber)))!=length(chromosome))
    stop("chromosome must contain integers between 1 and 25.")
  chromosome=unique(chromosome)
  if(length(chromosome)>25 || length(chromosome)==0)
    stop("chromosome contains too many elements.")
  chromosome=sort(chromosome)
  for(i in 1:length(chromosome))
    if( (chromosome[i]<1) || (chromosome[i]>25) )
      stop("chromosome must contain integers between 1 and 25.")
  
  #dataSetName
  if(missing(dataSetName))
    stop("dataSetName is missing.")
  if(!is.character(dataSetName))
    stop("dataSetName must be the name of an existing folder in rawData.")
  
 
  #check if we are in a normal-tumor study or in a single array study
  singleStudy=TRUE
  if(missing(normalTumorArray))
  {
    if(verbose)
      cat("No normalTumorArray specified.\n The copy-number signal will be extracted for all the specified data and the median of the dataset will be used as reference.\n")
  }
  else
  {
    if(verbose)
      cat("The copy-number signal will be extracted for all the specified data and the copy-number signal of normal DNA will be used as reference.\n")
    singleStudy=FALSE
  }
  
  #################### import the dataSet to have the name of all the files
  
  #path where find the CN data
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
  dataSet <- paste0(dataSetName,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
  
  #load CN
  dsC <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
  #print(dsC);
  
  ###### check listOfFiles
  pos=c()
  if(is.null(listOfFiles) || missing(listOfFiles))
  {
    listOfFiles=getNames(dsC)
    pos=1:length(dsC)
  }
  else
  {
    #check the format
    if(!is.vector(listOfFiles) || !is.character(listOfFiles))
      stop("listOfFiles must be a vector of string.")
    listOfFiles=unique(listOfFiles)
    
    #check if all the files of listOfFiles are in the folder
    pos=sapply(listOfFiles,match,getNames(dsC))#position of the files of listOfFiles in the folder
    if(length(which(pos>0))!=length(pos))
      stop("Wrong name of files in listOfFiles")
  } 
  
  ################### check normalTumorArray
  if(!singleStudy)
  {
    #normalTumorArray
    if(is.character(normalTumorArray))
      normalTumorArray=read.csv(normalTumorArray)
    else
    {
      if(!is.data.frame(normalTumorArray))
        stop("normalTumorArray must be either the path to the normalTumorArray csv file or a data.frame containing the data.\n")
    }
    
    #check normalTumorArray
    if(!("normal"%in%colnames(normalTumorArray)) || !("tumor"%in%colnames(normalTumorArray)))
      stop("normalTumorArray doesn't contain a column \"normal\" or \"tumor\".\n")
    
    #check is the normalTumorArray file contains all the files
#     isArrayComplete=sapply(getNames(dsC),FUN=function(name,nameOfFiles){name%in%nameOfFiles},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
#     if(sum(isArrayComplete)!=length(isArrayComplete))
#       stop("normalTumorArray doesn't contain all the filenames of dataSetName.")
    
    if(!missing(normalTumorArray))
    {
      isArrayComplete=sapply(listOfFiles,FUN=function(name,listOfNames){name%in%listOfNames},c(as.character(normalTumorArray$normal),as.character(normalTumorArray$tumor)))
      if(sum(isArrayComplete)!=length(isArrayComplete))
        stop("normalTumorArray doesn't contain all the filenames you specified in listOfFiles parameter.")
    }
  }
  
  #if paired study, we keep the name of normal files
  normalFiles=NULL
  if(!singleStudy)
  {
    #if normal-tumor study, we need the tumor and normal files
      
    #we obtain the complementary files
    compFiles=getComplementaryFile(listOfFiles,normalTumorArray)
    allFiles=unique(c(listOfFiles,compFiles))
    
    #get the status ("normal" or "tumor") of each files
    status=getStatus(allFiles,normalTumorArray)
    
    #keep only the normal files
    normalFiles=allFiles[which(status=="normal")]
        
    rm(compFiles,allFiles)
  }  
  
  ########### END CHECK ARGUMENT
  
  
  #get names and position of the probes
  ugp <- getAromaUgpFile(dsC);
  unf <- getUnitNamesFile(ugp);
  
  #get the prefix of SNP probes
  platform <- getPlatform(ugp);
  if (platform == "Affymetrix") 
  {
    require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix");
    snpPattern <- "^SNP";
  } 
  else if (platform == "Illumina") 
  {
    snpPattern <- "^rs[0-9]";
  }
  else 
    throw("Unknown platform: ", platform);
  
  allCN=list()
  for(chr in chromosome)
  {
    units <- getUnitsOnChromosome(ugp, chromosome=chr);
    unitNames <- getUnitNames(unf,units=units);##names of the probes
    if(onlySNP)
    {
      #keep the SNP units
      units=units[grep(snpPattern,unitNames)]
      unitNames=unitNames[grep(snpPattern,unitNames)]
    }
    posChr <- getPositions(ugp, units=units);#positions of the probes
    #sort signal by position
    indSort=order(posChr)
    
    units=units[indSort]
    unitNames=unitNames[indSort]
    posChr=posChr[indSort]
    
    
    #get the normalized CN signal
    if(singleStudy)
      C=getCopyNumberSignalSingleStudy(dsC,units,chromosome,pos,ugp)
    else
      C=getCopyNumberSignalPairedStudy(dsC,units,normalTumorArray,chromosome,normalFiles,ugp)
    

    allCN[[paste0("chr",chr)]]=data.frame(chromosome=rep(chr,length(posChr)),position=posChr,C$CN[,,drop=FALSE],featureNames=unitNames)
  }

  return(allCN)
}

##########################################################################################################################################

#
# @title Get CN signal in case of single study
# @param dsC object return by AromaUnitTotalCnBinarySet$byName function
# @param units position to keep
# @param chromosome chromosome to extract
# @param indexOfFiles position (in dsC) of the files to extract
# 
# @return A list of 2 elements :
# \describe{
#   \item{CN}{A matrix. Each column contains the CN signal for a different profile.}
#   \item{sampleNames}{Names of the extracted signal.}
# }
# 
getCopyNumberSignalSingleStudy=function(dsC,units,chromosome,indexOfFiles,ugp)
{    
  #compute the median
  ceR <- getAverageFile(dsC, verbose=0)
  
  #reduce to the files from indexOfFiles
  dsC=extract(dsC,indexOfFiles)
  
  # Extract total CNs
  C <- extractMatrix(dsC,units=units);
  
  #CN median for units
  thetaR <- extractMatrix(ceR,units=units)
  
  #the ploidy depends of the gender for the cromosome 23 and 24
  if(chromosome==23 || chromosome ==24)
  {
    gender=findGender(getName(dsC),1:length(dsC),ugp)
    ploidy=rep(0,length(gender))
    ploidy[gender=="XY"]=1
    if(chromosome==23)
      ploidy[gender=="XX"]=2
    else
      ploidy[gender=="XX"]=NA
    
    C <- sapply(1:ncol(C),FUN=function(x){ploidy[x]*C[,x]/thetaR});
  }
  else
  {
    #normalization
    C <- apply(C,2,FUN=function(x){2*x/thetaR});
  }

  sampleNames=getNames(dsC)
  return(list(CN=C,sampleNames=sampleNames))
}

##########################################################################################################################################

#
# @title Get CN signal in case of normal-tumor study
# @param dsC object return by AromaUnitTotalCnBinarySet$byName function
# @param units position to keep
# @param normalTumorArray a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
# @param chromosome chromosome to extract
# @param normalFiles names of the normal files
# 
# @return A list of 2 elements :
# \describe{
#   \item{CN}{A matrix. Each column contains the CN signal for a different profile.}
#   \item{sampleNames}{Names of the extracted signal.}
# }
# 
getCopyNumberSignalPairedStudy=function(dsC,units,normalTumorArray,chromosome,normalFiles,ugp)
{  
  #id of normal files
  normalId=sapply(normalFiles,FUN=function(x,names){which(names==x)},getNames(dsC))
  
  #names of tumor files
  tumorFiles=getComplementaryFile(normalFiles,normalTumorArray)
  
  #id of tumor files
  tumorId=sapply(tumorFiles,FUN=function(x,names){which(names==x)},getNames(dsC))

  # Extract normal
  dsPairCNormal <- extract(dsC, normalId);  
  
  # Extract tumor
  dsPairCTumor <- extract(dsC, tumorId);  
  
  #tumor CN
  C <- extractMatrix(dsPairCTumor, units=units);
  
  #normalCN
  Cnormal <- extractMatrix(dsPairCNormal, units=units);
  
  #the ploidy depends of the gender for the cromosome 23 and 24
  if(chromosome==23 || chromosome ==24)
  {
    gender=findGender(getName(dsC),normalId,ugp)
    ploidy=rep(0,length(gender))
    ploidy[gender=="XY"]=1#male:chr23=chrX, chr24=chrY, the ploidy is 1
    if(chromosome==23)
      ploidy[gender=="XX"]=2#female:chr23=chrXX, the ploidy is 2
    else
      ploidy[gender=="XX"]=NA#female:chr24=nothing, the ploidy is NA
    
    #normalize by the normal sample
    C <- C/Cnormal
    #apply the ploidy
    C <- sapply(1:ncol(C),FUN=function(x){ploidy[x]*C[,x]});
  }
  else
  {
    #normalization of CN
    C <- 2*C/Cnormal;
  }

  sampleNames=getNames(dsC)[tumorId]
  
  return(list(CN=C,sampleNames=sampleNames))
}


##########################################################################################################################################

#find gender of a profile CEL file
#
# @param name of the dataset (it must correpond to a folder name in rawData folder.).
# @param indexFileToExtract index of the files to test
#
# @return a vector containing the gender ("XX" or "XY") of the specified files
findGender=function(dataSetName,indexFileToExtract,ugp)
{
  #path where find the data
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
  dataSet <- paste0(dataSetName,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
  
  #get the allele B fraction
  dsF <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
        
  # Identify units on ChrX and ChrY 
  units23 <- getUnitsOnChromosome(ugp, 23);
  units24 <- getUnitsOnChromosome(ugp, 24);
    
  gender=lapply(indexFileToExtract,FUN=function(index)
  {
    df=getFile(dsF,index)
  
    beta23=df[units23,1,drop=TRUE]
    beta24=df[units24,1,drop=TRUE]
    genderTemp <- callXXorXY(beta23, beta24, adjust=1.5, from=0, to=1);
  })
  
  return(unlist(gender))
}

##########################################################################################################################################

# check if a number is an integer
# @param x number
# @param tol tolerance
#
# @ return TRUE if the number is an integer, FALSE else
#
is.wholenumber=function(x, tol = .Machine$double.eps^0.5)  
{
  #if(!is.double(x))
  #  return(FALSE)
  
  abs(x - round(x)) < tol
}

##########################################################################################################################################

# get the normal files associated to a tumor files (or vice-versa)
#
# @param listOfFiles vector containing names of iles
# @param normalTumorArray a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
#
# @return a vector of the same size as listOfFiles containing the complementary files
#
getComplementaryFile=function(listOfFiles,normalTumorArray)
{
  sapply(listOfFiles,FUN=function(names,data)
  {
    ind=match(names,data$normal)
    if(!is.na(ind))
    {
      return(as.character(data$tumor[ind]))
    }
    else
    {
      ind=match(names,data$tumor)
      if(!is.na(ind))
        return(as.character(data$normal[ind]))
      else
        return(NA)
    }
  },normalTumorArray);
}

##########################################################################################################################################

# get the status of a file ("normal" or "tumor")
#
# @param listOfFiles vector containing names of iles
# @param normalTumorArray a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
#
# @return a vector containing the status of each files of listOfFiles
#
getStatus=function(listOfFiles,normalTumorArray)
{
  sapply(listOfFiles,FUN=function(names,data)
  {
    if(names%in%data$normal)
    {
      return("normal")
    }
    else
    {
      if(names%in%data$tumor)
        return("tumor")
      else
        return(NA)
    }
  },normalTumorArray);
}

##########################################################################################################################################

# get the normal files associated to a tumor files (or vice-versa)
#
# @param listOfFiles vector containing names of iles
# @param normalTumorArray a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
#
# @return a vector of the same size as listOfFiles containing the complementary files
#
getId=function(listOfFiles,normalTumorArray)
{
  sapply(listOfFiles,FUN=function(names,data)
  {
    ind=match(names,data$normal)
    if(!is.na(ind))
    {
      return(ind)
    }
    else
    {
      ind=match(names,data$tumor)
      if(!is.na(ind))
        return(ind)
      else
        return(NA)
    }
  },normalTumorArray);
}