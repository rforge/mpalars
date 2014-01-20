# 
# @title genotypage calls
# @param dataSetName name of the data set.
# @param normalTumorArray  A data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
# @param singleArray if TRUE, it is a single Array analysis, no need of normalTumorArray
# @param plot, if TRUE, graphics will be created in the case of a singleArray study
#
# @export
naiveGenotypeCalls<-function(dataSetName,normalTumorArray,singleArray,plot)
{  
  require(aroma.cn)
  ##Setup
  log <- verbose <- Arguments$getVerbose(-1, timestamp=TRUE);
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
 
  #name of the folder containing the fracB of dataSetName
  dataSet <- paste0(dataSetName,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load the raw data set 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  ds <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
  dsN=ds;
  #print(dsN);

  #in case of normal-tumor study, we extract the normal file and launch the genotype call only fot them
  if(!singleArray)
  {
    #### Extract the normals
    
    #get the type of each file ("tumor" or "normal")
    types <- getStatus(getNames(ds),normalTumorArray)
    
    #index of the normal file
    normals <- grep("normal", types);
    
    #we keep normal files
    dsN <- extract(ds, normals);
  }  
  
  #path for the output  
  fullname <- paste(c(getFullName(dsN), "NGC"), collapse=",");
  chipType <- getChipType(dsN, fullname=FALSE);
  outPath <- file.path("callData", fullname, chipType);
  
  #
  df <- getFile(dsN, 1);
  units <- seq(length=nbrOfUnits(df));
  rm(df);
  
  #adjust parameter for callNaiveGenotypes function
  adjust <- 1.5;
  
  # Identify units on ChrX and ChrY (will require a specific treatment)
  ugp <- getAromaUgpFile(dsN);
  units23 <- getUnitsOnChromosome(ugp, 23);
  is23 <- is.element(units, units23);
  units24 <- getUnitsOnChromosome(ugp, 24);
  is24 <- is.element(units, units24);
  
  #launch genotype call for each file in the data set
  for(numberOfFile in 1:(length(dsN)))
  {   
    #get the data of the numberOfFile-th file
    dfN <- getFile(dsN, numberOfFile);
    
    #create the output filename
    tags <- getTags(dfN);
    tags <- setdiff(tags, "fracB");
    
    #filename for genotype call
    tags <- c(tags, "genotypes");
    fullname <- paste(c(getName(dfN), tags), collapse=","); 
    filename <- sprintf("%s.acf", fullname);
    gcPathname <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE);
    
    #filename for confidence score
    csTags <- c(tags, "confidenceScores");
    fullname <- paste(c(getName(dfN), csTags), collapse=",");
    filename <- sprintf("%s.acf", fullname);
    csPathname <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE);
    
    
    if (isFile(gcPathname) && isFile(csPathname)) 
    {
      next;
    }
    
    #get the fracB for the specified units
    betaN <- dfN[units,1,drop=TRUE];
    
    # Call gender
    gender <- callXXorXY(betaN[is23], betaN[is24], adjust=adjust, from=0, to=1);
    
    # initialization of mu (genotype) and cs (confident score)
    naValue <- as.double(NA);
    fit <- NULL;
    mu <- rep(naValue, times=length(units));
    cs <- rep(naValue, times=length(units));
    
    #male case
    if (gender == "XY")
    {
      # All but ChrX & ChrY in male
      isDiploid <- (!(is23 | is24));
      
      #look for diploid position
      use <- which(isDiploid);
      #launch genotype calls with copy number=2 for diploid position
      muT <- callNaiveGenotypes(betaN[use], cn=2, adjust=adjust, from=0, to=1,verbose=less(verbose,10));
      fit <- attr(muT, 'modelFit');
      mu[use] <- muT;
      
      #look for non diploid position
      use <- which(!isDiploid);
      #launch genotype calls with copy number=1 for non diploid position
      muT <- callNaiveGenotypes(betaN[use], cn=1, adjust=adjust, from=0, to=1,verbose=less(verbose,10));
      mu[use] <- muT;
    } 
    else ##female case
    {
      # All but ChrY in female
      isDiploid <- (!is24);
      
      #look for diploid position
      use <- which(isDiploid);
      #launch genotype calls with copy number=2 for diploid position
      muT <- callNaiveGenotypes(betaN[use], cn=2, adjust=adjust, from=0, to=1,verbose=less(verbose,10));
      fit <- attr(muT, 'modelFit');
      mu[use] <- muT;
    }
    
    # Translate genotype calls in fracB space to (AA,AB,BB)
    calls <- rep(as.character(NA), times=length(mu));
    calls[mu ==   0] <- "AA";
    calls[mu == 1/2] <- "AB";
    calls[mu ==   1] <- "BB";
    print(table(calls, exclude=NULL));
    
    # Calculate confidence scores
    a <- fit[[1]]$fitValleys$x[1];
    b <- fit[[1]]$fitValleys$x[2];
    cs[isDiploid] <- rowMins(abs(cbind(betaN[isDiploid]-a, betaN[isDiploid]-b)));
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Writing genotype calls (via temporary file)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pathname <- gcPathname;
    pathnameT <- sprintf("%s.tmp", pathname);
    nbrOfUnits <- nbrOfUnits(dfN);
    gfN <- AromaUnitGenotypeCallFile$allocate(pathnameT, platform=getPlatform(dfN), chipType=getChipType(dfN), nbrOfRows=nbrOfUnits);
    footer <- readFooter(gfN);
    footer$method <- "NaiveGenotypeCaller";
    writeFooter(gfN, footer);
    rm(footer);

    updateGenotypes(gfN, units=units, calls=calls);
    rm(calls);
    
    res <- file.rename(pathnameT, pathname);
    if (!isFile(pathname))
    {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    if (isFile(pathnameT)) 
    {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    rm(pathnameT);
    
    gfN <- AromaUnitGenotypeCallFile(pathname);
    
    #print(gfN);
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Writing confidence scores (via temporary file)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pathname <- csPathname;
    pathnameT <- sprintf("%s.tmp", pathname);
    nbrOfUnits <- nbrOfUnits(dfN);
    csfN <- AromaUnitSignalBinaryFile$allocate(pathnameT, platform=getPlatform(dfN), chipType=getChipType(dfN), nbrOfRows=nbrOfUnits, types="double", size=4, signed=TRUE);
    footer <- readFooter(csfN);
    footer$method <- "NaiveGenotypeConfidenceScoreEstimator";
    writeFooter(csfN, footer);
    rm(footer);
    
    csfN[units, 1] <- cs
    rm(cs);
    
    res <- file.rename(pathnameT, pathname);
    if (!isFile(pathname)) 
    {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    if (isFile(pathnameT)) 
    {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    rm(pathnameT);
     
    cfN <- AromaUnitSignalBinaryFile(pathname);
    #print(cfN);
    
    
    #load the data
#     gcN <- AromaUnitGenotypeCallSet$byName(dataSet, tags="NGC", chipType="*");
#     print(gcN);
#     
#     csN <- AromaUnitSignalBinarySet$byName(dataSet, tags="NGC", chipType="*", pattern="confidenceScores", paths="callData");
#     print(csN)
  }

  if(singleArray && plot)
    SingleStudyPlot(dataSetName)
    
}

################################################


# plot function in case of single study
#
# @param dataFolder name of the folder containing the data
#
SingleStudyPlot=function(dataFolder)
{ 
  log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
  
  #path where download cn and fracB 
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
  dataSet <- paste0(dataFolder,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
  
  #load CN and fracB data
  dsC <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
  dsFracB <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
  
  # Identify available genotype calls
  rootPath <- "callData";
  rootPath <- Arguments$getReadablePath(rootPath); 
  genotypeTag <- "NGC";
  gsN <- AromaUnitGenotypeCallSet$byName(dataSet, tags=genotypeTag, chipType="*");
  
  #print(gsN);
  
  #we compute the median of CN
  ceR <- getAverageFile(dsC, verbose=verbose)
  
  #names of the probes
  ugp <- getAromaUgpFile(dsC);
  unf <- getUnitNamesFile(ugp);
  sampleNames=getNames(dsFracB)
  #human genome reference
  hg=getTags(ugp)[grep("hg",getTags(ugp))]
  
  # prefix of SNP
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
  
  #loop on all the files
  for(i in 1:length(sampleNames))
  {
    #id of the i-th file
    sampleName=sampleNames[i]
    
    #CN and fracB for the i-th file
    dsCtemp=extract(dsC,i)
    dsFracBtemp=extract(dsFracB,i)
    gsNtemp=extract(gsN,i)
    gsNtemp=getFile(gsNtemp,1)
        
    cat("Saving graphics for sample",sampleName,"\n")
    
    for(chromosome in 1:25)
    {
      chrTag <- sprintf("Chr%02d", chromosome);
  
      testGraph=TRUE
      figName <- sprintf("%s,%s", sampleName, chrTag);
      figPath <- Arguments$getWritablePath(paste0("figures/",dataFolder,"/signal/"));
      pathname <- filePath(figPath, sprintf("%s.png", figName));
      
      #check existence of the file, if the file already exists, we don't create it
      testGraph=!isFile(pathname)

      if(testGraph)
      {      
        
        units <- getUnitsOnChromosome(ugp, chromosome=chromosome);
        unitNames <- getUnitNames(unf,units=units);##names of the probes
        posChr <- getPositions(ugp, units=units);#positions of the probes      
        
        
        ##########################################################
        #####################  COPY NUMBER  ######################
        ##########################################################
        
        # Extract total CNs
        C <- extractMatrix(dsCtemp,units=units);##units
        
        #CN median for units
        thetaR <- extractMatrix(ceR,units=units)
        
        #normalization
        C <- 2*C[,1]/thetaR;
        
        ##########################################################
        ################  FRACTION ALLELE B  #####################
        ##########################################################
        
        # Identify SNP units (fracB is only defined on SNP probes)
        keep <- (regexpr(snpPattern, unitNames) != -1);
        SNPunits <- units[keep];
        posSNP <- getPositions(ugp, units=SNPunits);#SNP position
        
        
        # Extract Allele B fractions
        df <- getFile(dsFracBtemp,1);
        beta <- df[SNPunits,1,drop=TRUE];
        beta <- as.data.frame(beta);
        beta <- as.matrix(beta);
        names <- colnames(beta);
        names <- "tumor";  
        
        ####################### SORTIE GRAPHIQUE
    
        # Plot dimensions
        x <- posChr/1e6;
        xlim <- range(x, na.rm=TRUE);
        xlab <- "Position (Mb)";
        
        width <- 1280;
        aspect <- 0.6*1/3;
        
        # Plot total CNs
        ylim <- c(0,6);
        ylab <- "Copy number";
        
        fig <- devNew("png", pathname, label=figName, width=width, height=2*aspect*width);
        par(mfrow=c(2,1))
        
        plot(NA,ylim=c(0,6),xlim=xlim, xlab=xlab, ylab=ylab, axes=FALSE,pch=".")
        axis(side=1);
        axis(side=2, at=c(0,2,4,6));
        points(x, C, pch=".");
        label <- sprintf("%s", sampleName);
        stext(side=3, pos=0, label);
        stext(side=3, pos=1, chrTag);
        
        # Plot Allele B fractions
        x <- posSNP/1e6;
        ylim <- c(-0.05,1.05);
        ylim <- c(-0.1,1.1);
        ylab <- "Allele B Fraction"; 

        calls=gsNtemp[SNPunits,1,drop=TRUE]
        cols <- as.integer(calls != 1) + 1L;
        
        name <- names;
        plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE);
        axis(side=1);
        axis(side=2, at=c(0,1/2,1));
        points(x, beta, pch=".",col=cols);
        label <- sprintf("%s (%s)", sampleName, name);
        stext(side=3, pos=0, label);
        stext(side=3, pos=1, chrTag); 
        devDone();
        cat("*")
      } 
    }##end loop chromosome 
    cat("\n")
  }##end loop sample
}##end function