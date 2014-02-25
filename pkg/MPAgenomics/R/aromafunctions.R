#'
#' @title Get the contents of a data folder
#' 
#' @description Get the cel files of the specified dataSetName
#' 
#' @param dataSetName Name of a folder in rawData folder
#' @param chipType Name of the chip type used for the data
#' 
#' @return The filenames of all the files in rawData/dataSetName/chipType
#' 
#' @details if chipType is not provided, the function returns the files for the first chip (in the alphabetic order)
#' 
#' @author Quentin Grimonprez
#' 
#' @export
getListOfFiles=function(dataSetName,chipType)
{
  if(!("annotationData"%in%list.files()))
    stop("There is no annotationData folder in the current architecture.")
  
  if(!("rawData"%in%list.files()))
    stop("There is no rawData folder in the current architecture.")
  
  #dataSetName
  if(missing(dataSetName))
    stop("dataSetName is missing")
  if(!is.character(dataSetName))
    stop("dataSetName must be a string.")
  
  #chipType
  if(missing(chipType))
    chipType=list.files(paste0("rawData/",dataSetName,"/"))[1]
  if(!is.character(chipType))
    stop("chipType must be a string.")
  
  files=list.files(paste0("rawData/",dataSetName,"/",chipType,"/"))
  withoutextension=gsub(".CEL$","",files)
  withoutextension=gsub(".cel$","",withoutextension)
  return(withoutextension)
}

#'
#' Create a folder in "annotationData/chipTypes" and copy the specified files in this folder.
#'
#' @title add a new chip type to the existing aroma architecture
#' 
#' @param chipType Name of the new chiptype to add.
#' @param chipPath Path to the files to add.
#' 
#' @author Quentin Grimonprez
#' 
#' @export
addChipType=function(chipType,chipPath)
{
  if(!("annotationData"%in%list.files()))
    stop("There is no annotationData folder in the current architecture.")
    
  existingChip=list.files("annotationData/chipTypes")
  
  if(chipType%in%existingChip)
    stop(paste0("A ",chipType," folder already exits."))
  
  copyChipFiles(chipPath,chipType,".",TRUE)
  
}


#'
#' Create a folder in "rawData" and copy the specified files in this folder.
#'
#' @title add a new data-set to the existing aroma architecture
#' 
#' @param dataSetName Name of the new data-set to add.
#' @param dataPath Path to the files to add.
#' @param chipType Name of the chip used for the data.
#' 
#' @author Quentin Grimonprez
#' 
#' @export
addData=function(dataSetName,dataPath,chipType)
{
  if(!("annotationData"%in%list.files()))
    stop("There is no annotationData folder in the current architecture.")
  
  if(!("rawData"%in%list.files()))
    stop("There is no rawData folder in the current architecture.")
  
  existingChip=list.files("annotationData/chipTypes")
  
  if(!(chipType%in%existingChip))
    stop(paste0("The ",chipType," folder does not exit."))
  
  copyDataFiles(dataSetName,dataPath,chipType,".",TRUE)
  
}
  