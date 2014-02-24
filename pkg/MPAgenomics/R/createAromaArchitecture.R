#' Create the architecture required by aroma packages
#' 
#' @title Create aroma architecture
#' @param dataSetName name of the data set
#' @param chipType type of the chip used for obtaining the data 
#' @param path path where create folders
#' @param verbose if TRUE, print details of the process
#'
#' @author Quentin Grimonprez
#'
#' @details This function creates the following architecture:
#'Architecture to create:
#'   <path>
#'     +- annotationData/
#'     |  +- chipTypes/
#'     |     +- <chipType>/ <-- must match exactly the name of the CDF file (fullname minus tags)
#'     |        +- CDF file(s) and other annotation (possibly subdirectories)
#'     |
#'     +- rawData/
#'        +- <dataSetName>/
#'           +- <chipType>/ <-- must match exactly a chip type folder under annotationData/
#'              +- CEL files
#'
#' @examples
#' createEmptyArchitecture("data1","GenomeWideSNP_6",".",TRUE)
#'
#' @export
createEmptyArchitecture=function(dataSetName,chipType,path=".",verbose=TRUE)
{
  #check arguments
  if(missing(dataSetName))
    stop("dataSetName is missing.")
  if(missing(chipType))
    stop("chipType is missing.")
  if(!is.character(dataSetName))
    stop("dataSetName must be a string.")
  if(!is.character(chipType))
    stop("chipType must be a string.")
  
  #create rawData
  if(!file.exists(paste0(path,"/rawData/",dataSetName,"/",chipType)))
  {
    dir.create(paste0(path,"/rawData/",dataSetName,"/",chipType),recursive=TRUE)
    ok1=file.access(paste0(path,"/rawData/",dataSetName,"/",chipType), mode = 0)
    if(ok1 != 0)
      stop(paste0("Can not create \"",path,"/rawData/",dataSetName,"/",chipType,"\" folder."))
    else
    {
      if(verbose)
        cat(gsub("//","/",paste0("Folder \"",path,"/rawData/",dataSetName,"/",chipType,"\" created.")),"\n")
    }
  }
  else
    cat("The folder",paste0(path,"/rawData/",dataSetName,"/",chipType),"already exits.")
  
  #create annotationData
  if(!file.exists(paste0(path,"/annotationData/chipTypes/",chipType)))
  {
    dir.create(paste0(path,"/annotationData/chipTypes/",chipType),recursive=TRUE)
    
    #test existence of the created folders
    ok2=file.access(paste0(path,"/annotationData/chipTypes/",chipType), mode = 0)

    if(ok2 != 0)
      stop(paste0("Can not create \"",path,"/annotationData/chipTypes/",chipType,"\" folder."))  
    else
    {
      if(verbose)
        cat(gsub("//","/",paste0("Folder \"",path,"/annotationData/chipTypes/",chipType,"\" created.")),"\n")
    }
  }
  else
    cat("The folder",paste0(path,"/annotationData/chipTypes/",chipType),"already exits.")
  
  return(invisible(TRUE))
}

# @title Copy data files in the appropriate folder
# @param dataSetName name of the data set
# @param dataSetPath path where to find the data files
# @param chipName type of the chip used for the data 
# @param path path where to find rawData folder
# @param verbose if TRUE, print details of the process
#
# @author Quentin Grimonprez

#@note WARNING : Files are copied, not moved. 
copyDataFiles=function(dataSetName,dataSetPath,chipName,path,verbose)
{
  if(verbose)
    cat("\nCopy data files:\n")
  
  #list all the files in the dataSet
  fileNames=list.files(dataSetPath)
  
  lapply(fileNames,FUN=function(fileName)
    {
      #check the extension
      if(toupper(substr(fileName,nchar(fileName)-3,nchar(fileName)))==".CEL")
      {
        #copy the file in rawData/<dataSetName>/<chipName>/
        if(file.exists(paste0(path,"/rawData/",dataSetName,"/",chipName,"/",fileName)))
        {
          if(verbose)
            cat("File",fileName," already exists.\n")
        }
        else
        {
          if(verbose)
            cat(paste0("Copying ",fileName,"..."))
          file.copy(paste0(dataSetPath,"/",fileName),paste0(path,"/rawData/",dataSetName,"/",chipName,"/",fileName),overwrite=FALSE)
          if(verbose)
            cat("\t DONE.\n")
        }

      }
    }
  )#end lapply
  
  return(invisible(TRUE)) 
}

# @title Copy chip files in the appropriate folder
# @param pathToChipFiles path where to find the chip files
# @param chipName type of the chip used for the data 
# @param path path where to find rawData folder
# @param verbose if TRUE, print details of the process
#
#
#@note WARNING : Files are copied, not moved.  Overwrite existing files with the same name.
# All the chip files have the form : chipName,Tags.extension.
# Extension can be ugp, ufl, cdf or acs
#
# @author Quentin Grimonprez
#
copyChipFiles=function(pathToChipFiles,chipName,path,verbose)
{
  if(verbose)
    cat("\nCopy chip files:\n")
  
  #valid extensions for chip files
  validExtension=c(".ugp",".ufl",".cdf",".acs")
  
  #list all the files in the dataSet
  fileNames=list.files(pathToChipFiles)
    
  lengthChipName=nchar(chipName)
  
  lapply(fileNames,FUN=function(fileName)
  {
    #check the extension
    if(tolower(substr(fileName,nchar(fileName)-3,nchar(fileName)))%in%validExtension)
    { 
      #check if it is a file for the good chip
      if(chipName==substr(fileName,1,nchar(chipName)))
      {
        #copy the file in annotationData/<chipType>/
        if(verbose)
          cat(paste0("Copying ",fileName,"..."))
        file.copy(paste0(pathToChipFiles,"/",fileName),paste0(path,"/annotationData/chipTypes/",chipName,"/",fileName),overwrite=TRUE)
        
        #change ".Full" in ",Full" in cdf files
        if(substr(fileName,nchar(fileName)-8,nchar(fileName)-4)==".Full")
        {
          newfileName=paste0(substr(fileName,1,nchar(fileName)-9),",Full",substr(fileName,nchar(fileName)-3,nchar(fileName)))
          renameFile(paste0(path,"/annotationData/chipTypes/",chipName,"/",fileName),paste0(path,"/annotationData/chipTypes/",chipName,"/",newfileName))          
        }
        
        if(verbose)
          cat("\t DONE.\n")
      }
      else
      {
        if(verbose)
          cat("File name: ",fileName,"doesn't correspond with specified chip name :",chipName ,"\n")
      }
    }

  }
  )#end lapply
  
  return(invisible(TRUE))   
}

#check if the chipfiles are good
#@author Quentin Grimonprez
.checkChipType=function(chipType,tag,path)
{
  if(!suppressPackageStartupMessages(require("aroma.affymetrix", quietly=TRUE) ) )
  {
    cat("Package not found: aroma.affymetrix. For download it:\n")
    cat("source(\"http://www.braju.com/R/hbLite.R\")\n")
    cat(" hbLite(\"sfit\")\n")
    cat("source(\"http://bioconductor.org/biocLite.R\")\n")
    cat("biocLite(\"affxparser\")\n")
    cat("source(\"http://aroma-project.org/hbLite.R\")\n")
    cat("hbInstall(\"aroma.affymetrix\")\n")
  }
  else
    cat("Package aroma.affymetrix loaded.\n")
  actualPath=getwd()
  setwd(path)
  result <- try(cdf <- AffymetrixCdfFile$byChipType(chipType, tags=tag),silent=TRUE)
  if(class(result)[1]=="try-error")
  {
    setwd(actualPath)
    stop(geterrmessage())
    #stop(paste0("Problem with cdf files for the chip ",chipType," with tag ",tag))
  }

  result <- try(gi <- getGenomeInformation(cdf),silent=TRUE)
  if(class(result)[1]=="try-error")
  {
    setwd(actualPath)
    stop(geterrmessage())
    #stop(paste0("No ugp files for the chip ",chipType," with tag ",tag))
  }

  result <- try(si <- getSnpInformation(cdf),silent=TRUE)
  if(class(result)[1]=="try-error")
  {
    setwd(actualPath)
    stop(geterrmessage())
    #stop(paste0("No ufl files for the chip ",chipType," with tag ",tag))
  }

  result <- try(acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE)),silent=TRUE)
  if(class(result)[1]=="try-error")
  {
    setwd(actualPath)
    stop(geterrmessage())
    #stop(paste0("No acs files for the chip ",chipType," with tag ",tag))
  }

  
  setwd(actualPath)
  
  return(invisible(TRUE))
}


#' Create the architecture required by aroma.* packages and copy files into created folders.
#' 
#' @title Create aroma architecture and copy files
#' 
#' @param dataSetName name of the data set
#' @param chipType type of the chip used for the data 
#' @param dataSetPath Path to the folder containing the data CEL files
#' @param chipFilesPath Path to the folder containing the chip files
#' @param path path where create folders
#' @param verbose if TRUE, print details of the process
#'
#' @seealso copyChipFiles, copyDataFiles, createAromaArchitecture
#' 
#' @examples
#' #createArchitecture("data1","GenomeWideSNP_6",dataSetPath="./data/",chipFilesPath="./chip/",path="./aroma",verbose=TRUE)
#' 
#' @author Quentin Grimonprez
#' 
#' @export
createArchitecture=function(dataSetName,chipType,dataSetPath,chipFilesPath,path=".",verbose=FALSE)
{
  if(!suppressPackageStartupMessages(require("aroma.affymetrix", quietly=TRUE) ) )
  {
    cat("Package not found: aroma.affymetrix. For download it:\n")
    cat("source(\"http://www.braju.com/R/hbLite.R\")\n")
    cat(" hbLite(\"sfit\")\n")
    cat("source(\"http://bioconductor.org/biocLite.R\")\n")
    cat("biocLite(\"affxparser\")\n")
    cat("source(\"http://aroma-project.org/hbLite.R\")\n")
    cat("hbInstall(\"aroma.affymetrix\")\n")
  }
  else
    cat("Package aroma.affymetrix loaded.\n")
  
   createEmptyArchitecture(dataSetName,chipType,path,verbose)
   copyChipFiles(chipFilesPath,chipType,path,verbose)
  tag="Full"
  .checkChipType(chipType,tag,path)
   copyDataFiles(dataSetName,dataSetPath,chipType,path,verbose)

  
  return(invisible(TRUE))  
}

