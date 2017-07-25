########################
# FUNCTIONS!!
########################

######################
# HARVEST DATA
######################

HarvestData <- function(masterDF,dirs){
  newdf <- masterDF
  newdf$ISRES <- NA
  newdf$ISEXT <- NA
  newdf$MAXFREQ <- NA
  
  
  setwd(dirs$RESULTS_DIR)
  nreps <- nrow(masterDF)
  i=1
  for(i in 1:nreps){
    filename <- sprintf("Rep%s_results.RData",i)
    if(file.exists(filename)){
      load(filename)
      attach(ResultsList)
        newdf$ISRES[i] <- ifelse((finalabund>1000)&(mean(finalfreq)>0.5),1,0)
        newdf$ISEXT[i] <- ifelse(finalabund==0,1,0)
        newdf$MAXFREQ[i] <- max(byYear$resfreq)
      detach(ResultsList)
    }
  }
  return(newdf)
}


######################
# MAKE WORKER: FOR PARALLELIZATION
######################

MakeWorker <- function(NYEARS, masterDF, dirs, fake){
  
  force(NYEARS)  # not sure why this is necessary but whatever
  
  force(masterDF)
  
  force(dirs)
  
  force(fake)
  
  #library(raster)   # for now, load key packages (too much work to reference explicitly)
  

         # note: could me more efficient not to load packages, but just to reference the packages explicitly
             # actually that's not right- don't need that
  
  setwd(dirs$GIT_DIR)
  source("PlagueResistanceEvolution_FUNCTIONS.R")     # load necessary functions
  
  #LoadPackages(env=environment())   # maybe try to avoid loading packages, for better parallelization
  
 # thisEnvironment <- environment()  # need to be able to reference the worker environment
 
  
  Worker <- function(i){
    DoSimulateResistancePar(rep=i,fake=fake)
  }
  
  return(Worker) 
  
   
}  # end closure "MakeWorker"



#######################
# DoSimulateResistance
#######################


DoSimulateResistancePar <- function(rep=1,fake=F){
  #plot(PlagueRaster)
  
  #env <- environment()    # surrogate for global environment
  
  ############
  ## USER-DEFINED VARIABLES
  ############
  
  dmat <- list()
  dmat[[1]] <- matrix(c(1,0,0, 1,0,0), nrow=2,ncol=3,byrow = T)  # recessive
  dmat[[2]] <- matrix(c(1,1,0, 1,0,0), nrow=2,ncol=3,byrow = T)  # gene 1 dominant
  dmat[[3]] <- matrix(c(1,0,0, 1,1,0), nrow=2,ncol=3,byrow = T)  # gene 2 dominant
  dmat[[4]] <- matrix(c(1,1,0, 1,1,0), nrow=2,ncol=3,byrow = T)  # dominant
  
  UserParams <- DefineUserParams(PER_SUITABLE=masterDF$PER_SUITABLE[rep],SNUGGLE=masterDF$SNUGGLE[rep],NFOCI=1,MAXDISPERSAL=500,BASELINE_DISPERSAL=0.05,
                                 MAXDISPERSAL_PLAGUE = 1000,PLAGUE_DISPERSAL=masterDF$PLAGUE_DISPERSAL[rep], MAXDENS = masterDF$MAXDENS[rep],
                                 MINDENS = 15, BASELINE_MEANSURV = 0.6, BASELINE_PLAGUESURV=masterDF$BASELINE_PLAGUESURV[rep],
                                 BASELINE_PLAGUESURV_RESIST=0.5,BASELINE_MEANFEC=masterDF$BASELINE_MEANFEC[rep],
                                 FITNESS_COST=rep(masterDF$FITNESS_COST[rep],2),INITFREQ=rep(masterDF$INITFREQ[rep],2),DOMINANCE=dmat[[masterDF$DOMINANCE[rep]]])
  
  #assign(x="UserParams",value=UserParams, envir = env)
  
  InitList <- DoInitialization(UserParams,fake=fake)  #BaseLandscape
  PlagueRaster <- InitList$PlagueRaster
  YearsSincePlague <- InitList$YearsSincePlague
  #PopArray <- InitList$PopArray
  InitFreqList <- InitList$InitFreqList
  InitDensRaster <- InitList$InitDensRaster
  PlagueModel <- InitList$PlagueModel
  DispList <- InitList$DispList
  UserParams <- InitList$UserParams
  BaseLandscape <- InitList$BaseLandscape
  EnvCovs <- InitList$EnvCovs
  # UserParams <- get("UserParams",envir=env)
  # BaseLandscape <- get("BaseLandscape",envir=env)

  ResultsList <- SetUpResults(NYEARS,UserParams)
  
  ####################
  # START LOOP THROUGH YEARS
  ####################
  
  
  # names of important raster maps to save to file etc...
  #rasterNames  <- c("PopArray","NextPlagueSurvRaster","NextNormalSurvRaster","PlagueResistancePotentialRaster")   # Deprecate?
  
  
  # t=which(plagueyear)[1]
  t<-0
  t<-t+1
  for(t in 1:(NYEARS)){
    deviate <- rnorm(1)   #determine if this is a good year or a bad year (for now, survival and fecundity are perfectly correlated)
    #assign(x="deviate",value=deviate, envir = env)
    
    cv=UserParams$Popbio$CV_SURVIVAL   # set up for using the getYearVariate function
    #assign(x="cv",value=cv, envir = env)
    
    if(t==1){ 
      FreqList<-InitFreqList; 
      DensRaster<-InitDensRaster; 
      newFociRaster <- raster::reclassify(DensRaster,rcl=c(-Inf,Inf,0))     # initial conditions
      # assign(x="FreqList",value=FreqList, envir = env)
      # assign(x="DensRaster",value=DensRaster, envir = env)
      # assign(x="newFociRaster",value=newFociRaster, envir = env)
    }
    
    ##################
    # DENSITY INDEPENDENT SURVIVAL (including plague survival)
    ##################
    temp <- doSurvival(DensRaster,PlagueRaster,FreqList,UserParams,BaseLandscape,deviate,cv)  #env
    DensRaster <- temp$NewDensRaster  
    FreqList <- temp$FreqList
    
    #assign(x="DensRaster",value=DensRaster, envir = env)
    #FreqList <- get("FreqList",envir=env) # get("UserParams",envir=env)
    # plot(DensRaster)
    # plot(FreqList[["gene1"]])
    
    ################
    # REPRODUCTION
    ################
    
    DensRaster <- doReproduce(UserParams,DensRaster,PlagueRaster,deviate)    # TODO: make specific to each resistance type...?
    #assign(x="DensRaster",value=DensRaster, envir = env)
    # plot(DensRaster)
    # plot(FreqList[["gene1"]])
    
    ###############
    # DISPERSAL: Move individuals around the landscape (this takes a while!)
    ###############
    temp <- doDispersal(UserParams,DensRaster,PlagueRaster,newFociRaster,FreqList,BaseLandscape,DispList)
    DensRaster <- temp$newPop
    newFociRaster <- temp$newFociRaster
    FreqList <- temp$newFreqList
    
    #assign(x="DensRaster",value=DensRaster, envir = env)
    
    # raster::plot(DensRaster)    # good in t=1, not so much in t=2
    # raster::plot(FreqList[["gene2"]])
    
    ###############
    # ALLEE EFFECT: REMOVE POPULATIONS BELOW A MINIMUM ABUNDANCE THRESHOLD
    ###############
    #PopArray <- doAllee()
    #if(MINABUND>0) PopArray <- doAllee()   # don't need this! all low-dens individuals move out anyway...
    # $
    
    ###############
    # CLEAR excess individuals from cells (DD)
    
    DensRaster <- doDDSurvival(DensRaster,UserParams)
    #assign(x="DensRaster",value=DensRaster, envir = env)
    
    # raster::plot(DensRaster) 
    # plot(FreqList[["gene1"]])
    
    ##################
    # PLAGUE PRESSURE
    ##################
    
    YearsSincePlague <- YearsSincePlague + 1
    
    PlagueRaster <- doPlague(PlagueRaster=PlagueRaster,YearsSincePlague=YearsSincePlague,DensRaster=DensRaster,UserParams,PlagueModel,EnvCovs,fake,timestep=t)
    
    # raster::plot(PlagueRaster)
    
    ##################
    # update the years since plague variable...
    
    YearsSincePlague[PlagueRaster==1] <- 0 
    
    # raster::plot(YearsSincePlague)
    
    #assign(x="DensRaster",value=DensRaster, envir = env)
    
    # plot(PlagueRaster)
    
  
    
    
    # # plot(PopArray)
    # 
    #   plot(PopArray)    # okay
    #   plot(raster::reclassify(NextNormalSurvRaster,rcl=c(NA,NA,0)))   # okay in t=1,2
    #   plot(NextNormalSurvRaster)
    #   plot(raster::reclassify(NextPlagueSurvRaster,rcl=c(NA,NA,0)))    # okay t=1,2
    #   
    #   plot(newFociRaster)   # okay t=1
    #   
    #   plot(PlagueRaster)
    
    ###############
    # MAKE PLOTS
    ###############
    
    if(rep%%2==0) MakePlots(rep,t,BaseLandscape,DensRaster,PlagueRaster,FreqList)
    
    ###############
    # STORE RESULTS
    ###############
    
    ResultsList <- ComputeYearResults(ResultsList,DensRaster,FreqList,t)
    
  }   # end loop through time
  
  ###############
  # MAKE MOVIES
  ###############
  
  #if(rep%%2==0) MakeMovie(rep) 
  
  ################
  # FINAL RESULTS
  ################
  
  ResultsList <- ComputeFinalResults(ResultsList,DensRaster,FreqList)
  
  setwd(dirs$RESULTS_DIR)
  filename <- sprintf("Rep%s_results.RData",rep)
  save(ResultsList,file = filename)
  
}  # end main function DoSimulateResistance


#######################
# SET UP RESULTS
#######################

SetUpResults <- function(NYEARS,UserParams){
  results <- list()
  
  results$byYear <- list()
  
  results$byYear$totabund <- numeric(NYEARS)
  results$byYear$resfreq <- array(0,dim=c(NYEARS,UserParams[["Genetics"]]$NGENES))
  
  results$finalfreq<-numeric(UserParams[["Genetics"]]$NGENES)
  results$finalabund <- 0
  
  return(results)
}

######################
# COMPUTE YEARLY RESULTS
######################

ComputeYearResults <- function(ResultsList,DensRaster,FreqList,year){
  results <- ResultsList
  
  totabund <- raster::cellStats(DensRaster,stat='sum')
  results$byYear$totabund[year] <- totabund
  
  ngenes <- ncol(results$byYear$resfreq)
  
  for(i in 1:ngenes){
    freq <- raster::cellStats(FreqList[[i]],stat='mean')
    results$byYear$resfreq[year,i] <- freq
  }
  
  return(results)
}

######################
# COMPUTE FINAL RESULTS
######################

ComputeFinalResults <- function(ResultsList,DensRaster,FreqList){
  results <- ResultsList
  
  totabund <- raster::cellStats(DensRaster,stat='sum')
  results$finalabund <- totabund
  
  ngenes <- ncol(results$byYear$resfreq)
  
  for(i in 1:ngenes){
    freq <- raster::cellStats(FreqList[[i]],stat='mean')
    results$finalfreq[i] <- freq
  }
  
  return(results)
}

#######################
# MAKE PLOTS
#######################

MakePlots <- function(rep,t,BaseLandscape,DensRaster,PlagueRaster,FreqList){
  ################
  # MAKE PLOTS
  
  width = 500
  height= 500
  
  thisFIGS_DIR <- sprintf("%s\\rep%04d",dirs$FIGS_DIR2,rep)
  if(is.na(file.info(thisFIGS_DIR)[1,"isdir"])) dir.create(thisFIGS_DIR)
  
  # abundance figure
  setwd(thisFIGS_DIR)
  file = sprintf("AbundanceFig_year%03d.tif",t)
  tiff(file, width=width,height=height)
  raster::plot(BaseLandscape$patchRaster,col=gray(0.7),legend=F)
  #col = colorRampPalette(c("red","red"))(1)
  col = rgb(0,seq(0,1,length=10),0)
  raster::plot(raster::reclassify(DensRaster,rcl=c(-Inf,5,NA)),add=T,legend=T)
  #raster::plot(raster::reclassify(PlagueRaster,rcl=c(-Inf,0.01,NA)),col=rgb(1,0,0),add=T,alpha=0.5,legend=F)
  #raster::plot(raster::reclassify(NextPlagueSurvRaster,rcl=c(-Inf,0.001,NA)),col=heat.colors(10),add=T,legend=T)
  dev.off()
  
  # evolution figure 1
  setwd(thisFIGS_DIR)
  file = sprintf("AllGenesFreqFig_year%03d.tif",t)
  tiff(file, width=width*1.5,height=height)
  par(mfrow=c(1,2))
  raster::plot(BaseLandscape$patchRaster,col=gray(0.7),legend=F,main="Gene 1")
  #col = colorRampPalette(c("red","red"))(1)
  col2 = max(0,min(1,raster::maxValue(FreqList[["gene1"]])))
  col = rgb(0,seq(0,col2,length=10),0)
  raster::plot(raster::reclassify(FreqList[["gene1"]],rcl=c(-Inf,0.001,NA)),add=T,col=col,legend=T)
  
  raster::plot(BaseLandscape$patchRaster,col=gray(0.7),legend=F,main="Gene 2")
  #col = colorRampPalette(c("red","red"))(1)
  col2 = max(0,min(1,raster::maxValue(FreqList[["gene2"]])))
  col = rgb(0,seq(0,col2,length=10),0)
  raster::plot(raster::reclassify(FreqList[["gene2"]],rcl=c(-Inf,0.001,NA)),add=T,col=col,legend=T)
  #raster::plot(raster::reclassify(PlagueRaster,rcl=c(-Inf,0.01,NA)),col=rgb(1,0,0),add=T,alpha=0.5,legend=F)
  #raster::plot(raster::reclassify(NextPlagueSurvRaster,rcl=c(-Inf,0.001,NA)),col=heat.colors(10),add=T,legend=T)
  dev.off() 

  
  # plague figure
  setwd(thisFIGS_DIR)
  file = sprintf("PlagueFig_year%03d.tif",t)
  tiff(file, width=width,height=height)
  raster::plot(BaseLandscape$patchRaster,col=gray(0.7),legend=F)
  #col = colorRampPalette(c("red","red"))(1)
  col2 = max(0,min(1,raster::maxValue(PlagueRaster)))
  col = rgb(seq(0,col2,length=10),0,0)
  raster::plot(raster::reclassify(PlagueRaster,rcl=c(-Inf,0.01,NA)),add=T,col=col,legend=T)
  #raster::plot(raster::reclassify(PlagueRaster,rcl=c(-Inf,0.01,NA)),col=rgb(1,0,0),add=T,alpha=0.5,legend=F)
  #raster::plot(raster::reclassify(NextPlagueSurvRaster,rcl=c(-Inf,0.001,NA)),col=heat.colors(10),add=T,legend=T)
  dev.off()
  
  
  #gray.colors(10)
  
}   # end function MakePlots

###############
# MAKE MOVIE   (need ffmpeg and imagemagick installed)
###############


MakeMovie <- function(rep){
  thisMOVIE_DIR <- sprintf("%s\\rep%04d",dirs$MOVIE_DIR2,rep)
  if(is.na(file.info(thisMOVIE_DIR)[1,"isdir"])) dir.create(thisMOVIE_DIR)
  
  thisFIGS_DIR <- sprintf("%s\\rep%04d",dirs$FIGS_DIR2,rep)
  
  setwd(thisFIGS_DIR)
  
  #### remove figures that are not part of this simulation?
  
  fileyrs <- as.numeric(unlist(regmatches(list.files(), gregexpr("[[:digit:]]+", list.files()))))
  notthissim <- fileyrs > NYEARS
  toremove <- list.files()[notthissim]
  
  if(any(notthissim)) file.remove(toremove)
  
  ## NOTE: need command line like this: ffmpeg -f image2 -framerate 2 -i AbundanceFig_year%03d.tif -s 500x500 test.avi -y
  
  # MAKING THE REAL MOVIE HERE! USE IMAGE MAGICK AND FFMPEG SOFTWARE  (https://blogazonia.wordpress.com/2016/01/19/making-a-movie-with-r/)
  
  
  # create the movie
  cmd_abundmov <- paste0("ffmpeg -f image2 -framerate 2 -i AbundanceFig_year%03d.tif -s 500x500 ", 
                         sprintf("%s\\AbundanceMovie.avi",thisMOVIE_DIR)," -y")
  
  cmd_evolutionmov <- paste0("ffmpeg -f image2 -framerate 2 -i AllGenesFreqFig_year%03d.tif -s 800x500 ", 
                             sprintf("%s\\EvolutionMovie.avi",thisMOVIE_DIR)," -y")
  
  cmd_plaguemov <- paste0("ffmpeg -f image2 -framerate 2 -i PlagueFig_year%03d.tif -s 500x500 ", 
                          sprintf("%s\\PlagueMovie.avi",thisMOVIE_DIR)," -y")
  
  #sink(tempfile())
  system(cmd_abundmov,ignore.stdout = T,ignore.stderr = T)
  
  system(cmd_evolutionmov,ignore.stdout = T,ignore.stderr = T)
  
  system(cmd_plaguemov,ignore.stdout = T,ignore.stderr = T)
  
  #sink()
}  ### end function "MakeMovie"



#######################
# DO INITIALIZATION
#######################

DoInitialization <- function(UserParams,fake=FALSE){  # BaseLandscape
  #####################
  # INITIALIZE DISPERSAL   (for both plague and no plague... )   
  #####################
  
  temp <- InitializeDispersal(UserParams)   # env
  UserParams <- temp$UserParams    #get("UserParams",envir=env)
  DispList <- temp$DispList
  
  ########################
  # INITIALIZE LANDSCAPE
  ########################
  
  temp <- InitializeLandscape(solid=F,fake=fake,UserParams)   # generate patchmaps etc.  # env
  UserParams <- temp$UserParams   #get("UserParams",envir=env)
  BaseLandscape <- temp$BaseLandscape
  EnvCovs <- temp$EnvCovs
  
  ########################
  # GET PLAGUE MODEL
  ########################
  
  suppressWarnings(  PlagueModel <- GetPlagueModel(fake=fake)  )  # for now, use fake plague model- will be a statistical model
  
  #assign(x="PlagueModel",value=PlagueModel, envir = env)
  
  ########################
  # INITIALIZE POPULATION
  ########################
  
  InitDensRaster <- BaseLandscape$KRaster    # initialize abundance at carrying capacity
  # raster::plot(InitDensRaster)
  
  
  #######################
  # INITIALIZE ALLELE FREQUENCIES / RESISTANCE FACTORS [keep for now- will be multiple genes in the model somehow]
  #######################
  # NOTE: Some regions are more likely to evolve faster because they have greater percentages of those genes that can confer resistance. 
  
  InitFreqList <- GetInitFreqs(UserParams,BaseLandscape)
  
  #assign(x="InitFreqList",value=InitFreqList, envir = env)
  
  #####################
  # INITIALIZE POPULATION
  #####################
  
  ### Code block for pop starting from small loci
  
  # InitDensRaster2 <- raster::reclassify(BaseLandscape$patchIDRaster,rcl=c(-Inf,Inf,0))    # for testing
  # ndx <- sample(which(!is.na(InitDensRaster2@data@values)),size=3)
  # InitDensRaster2[ndx] <- 1000   # initialize population in random locations
  # InitDensRaster <- InitDensRaster2
  
  #PopArray2 <- InitDensRaster   # copy, for dispersal algorithm... 
  #assign(x="InitDensRaster",value=InitDensRaster, envir = env)
  
  # PopArray <- GetStructuredPop(InitDensRaster,InitFreqList,UserParams,env)
  # assign(x="PopArray",value=PopArray, envir = env)
  #raster::plot(PopArray)
  
  
  ######################
  # INITIALIZE PLAGUE PROCESS   [KTS: moving away from this and towards a purely statistical model]
  ######################
  #  for now, assume that plague hits at the patch level, and is a random process.
  
  # PROB_PLAGUE_YEAR <- 0.5 # probability that a plague event hits (landscape level?)   # for now, plague only hits one patch in a plague year
  # 
  # plagueyear = as.logical(rbinom(NYEARS,1,PROB_PLAGUE_YEAR))
  # 
  # plagueNow = floor(runif(NYEARS,1,BaseLandscape$nPatches+1)) * as.numeric(plagueyear)    ## which patch plagues out?
  
  PlagueRaster_template <- raster::reclassify(BaseLandscape$patchIDRaster,rcl=c(-Inf,Inf,0))   
  
  YearsSincePlague <- raster::reclassify(BaseLandscape$patchIDRaster,rcl=c(-Inf,Inf,10))
  
    #PlagueModel <- get("PlagueModel",envir=env)
  PlagueRaster <- doPlague(PlagueRaster=PlagueRaster_template, YearsSincePlague = YearsSincePlague,
                           DensRaster=InitDensRaster,UserParams,PlagueModel,EnvCovs,fake,timestep=0)
  
  #raster::plot(PlagueRaster)
  
  YearsSincePlague[PlagueRaster==1] <- 0

  #assign(x="PlagueRaster",value=PlagueRaster, envir = env)
  
  init.list <- list()
  
  init.list$PlagueRaster <- PlagueRaster
  init.list$YearsSincePlague <- YearsSincePlague
  #init.list$PopArray <- PopArray
  init.list$InitFreqList <- InitFreqList
  init.list$InitDensRaster <- InitDensRaster
  init.list$PlagueModel <- PlagueModel
  
  init.list$DispList <- DispList
  
  init.list$UserParams <- UserParams
  
  init.list$BaseLandscape <- BaseLandscape
  
  init.list$EnvCovs <- EnvCovs
  
  return(init.list)
  
  ###########
  # SAVE VARS TO GLOBAL ENVIRONMENT
  ###########
  

  
  # raster::plot(InitFreqList)
  
  # raster::plot(BaseLandscape$patchRaster)
  # raster::plot(BaseLandscape$patchIDRaster)
}




###########
## FUNCTION "loadPackage"
##
## GENERIC FUNCTION FOR INSTALLING/LOADING PACKAGES FROM CRAN
##
###########

loadPackage <- function(pkg, env=.GlobalEnv){
  
  if(pkg %in% rownames(installed.packages()) == FALSE) {suppressMessages(suppressWarnings(install.packages(pkg)))}
  eval(parse(text=sprintf("suppressMessages(suppressWarnings(require(%s)))",pkg)), envir= env)
  
}

###########
## FUNCTION "loadPackage"
##
## GENERIC FUNCTION FOR INSTALLING/LOADING SOURCE CODE FROM GITHUB
##
###########

source_github <- function(baseurl,scriptname) {
  # load package
  suppressMessages(suppressWarnings(require(RCurl)))
  
  # read script lines from website
  url <- sprintf("%s%s",baseurl,scriptname)
  script <- getURL(url, ssl.verifypeer = FALSE)
  
  script <- gsub("\r\n", "\n", script)     # get rid of carriage returns (not sure why this is necessary...)
  
  # parse lines and evaluate in the global environement
  eval(parse(text = script), envir= .GlobalEnv)
}



##########
## SET UP WORKSPACE
##########



SetUpDirectories <- function(){
  
  dirs <- list()
  
  if(KEVIN_LAPTOP) dirs$BASE_DIR <- "C:\\Users\\Kevin\\Dropbox\\PlagueModeling\\ResistanceEvolution"
  if(KEVIN_OFFICEPC) dirs$BASE_DIR <- "E:\\Dropbox\\PlagueModeling\\ResistanceEvolution"
  if(KEVIN_LAPTOP2) dirs$BASE_DIR <- "C:\\Users\\KevinT_Kevin\\Dropbox\\PlagueModeling\\ResistanceEvolution"
  if(MIRANDA) dirs$BASE_DIR <- "C:\\PlagueResistance"

  if(KEVIN_OFFICEPC) dirs$BASE_DIR2 <- "E:\\ResistanceEvolution"
  if(KEVIN_LAPTOP) dirs$BASE_DIR2 <- "C:\\ResistanceEvolution" 
  if(KEVIN_LAPTOP2) dirs$BASE_DIR2 <- "C:\\ResistanceEvolution"
  if(MIRANDA) dirs$BASE_DIR2 <- "C:\\PlagueResistance"
  
  dirs$plaguemod <- list()
  
  if(KEVIN_LAPTOP2) dirs$plaguemod$rootDir <- "C:\\Users\\KevinT_Kevin\\Dropbox\\PlagueModeling"
  if(KEVIN_OFFICEPC) dirs$plaguemod$rootDir <- "E:\\Dropbox\\PlagueModeling"
  if(MIRANDA) dirs$plaguemod$rootDir <- "C:\\PlagueResistance"
  
  
  dirs$plaguemod$ScriptDir <- paste(dirs$plaguemod$rootDir,"\\Rscript",sep="")
  dirs$plaguemod$DataDir <- paste(dirs$plaguemod$rootDir,"\\PrairieDogData",sep="")
  dirs$plaguemod$CovDir <- paste(dirs$plaguemod$rootDir,"\\Covariates",sep="")
  dirs$plaguemod$MODISDir <- paste(dirs$plaguemod$CovDir, "\\MODIS", sep="\\")
  
  # assign(x="BASE_DIR",value=BASE_DIR, envir = env)
  # assign(x="BASE_DIR2",value=BASE_DIR2, envir = env)
  
  
  ############
  ## SET UP WORKSPACE (define global variables)
  ############
  
  # RSCRIPT_DIR <- sprintf("%s\\Rscripts",GIT_DIR)
  dirs$DATA_DIR <- sprintf("%s\\Data",dirs$BASE_DIR)
  if(is.na(file.info(dirs$DATA_DIR)[1,"isdir"])) dir.create(dirs$DATA_DIR)
  #assign(x="DATA_DIR",value=DATA_DIR, envir = env)
  
  dirs$FIGS_DIR <- sprintf("%s\\RawFigs",dirs$BASE_DIR)
  if(is.na(file.info(dirs$FIGS_DIR)[1,"isdir"])) dir.create(dirs$FIGS_DIR)
  #assign(x="FIGS_DIR",value=FIGS_DIR, envir = env)
  
  dirs$FIGS_DIR2 <- sprintf("%s\\RawFigs",dirs$BASE_DIR2)
  if(is.na(file.info(dirs$FIGS_DIR2)[1,"isdir"])) dir.create(dirs$FIGS_DIR2)
  #assign(x="FIGS_DIR2",value=FIGS_DIR2, envir = env)
  
  dirs$MOVIE_DIR <- sprintf("%s\\Movies",dirs$BASE_DIR)
  if(is.na(file.info(dirs$MOVIE_DIR)[1,"isdir"])) dir.create(dirs$MOVIE_DIR)
#  assign(x="MOVIE_DIR",value=MOVIE_DIR, envir = env)
  
  dirs$MOVIE_DIR2 <- sprintf("%s\\Movies",dirs$BASE_DIR2)
  if(is.na(file.info(dirs$MOVIE_DIR2)[1,"isdir"])) dir.create(dirs$MOVIE_DIR2)
  #assign(x="MOVIE_DIR2",value=MOVIE_DIR2, envir = env)

  dirs$RESULTS_DIR <- sprintf("%s\\Results",dirs$BASE_DIR2)
  if(is.na(file.info(dirs$RESULTS_DIR)[1,"isdir"])) dir.create(dirs$RESULTS_DIR)
  #assign(x="MOVIE_DIR2",value=MOVIE_DIR2, envir = env)
    
  dirs$GIT_DIR <- GIT_DIR
  
  #setwd(dirs$DATA_DIR)
  
  #assign(x="dirs",value=dirs, envir = env)
  #LoadPackages(env)  # load all packages
  
  return(dirs)
}


#################
# LOAD PACKAGES
#################


LoadPackages <- function(env=.GlobalEnv){
  #############################
  #  LOAD PACKAGES
  #############################
  loadPackage("raster",env=env)           # note: loadPackage should install the package from CRAN automatically if it is not already installed
  loadPackage("secr",env=env)
  loadPackage("igraph",env=env)
  loadPackage("parallel",env=env)
  loadPackage("doParallel",env=env)
  loadPackage("RCurl",env=env)
  loadPackage("lhs",env=env)
  
}



###########
## FUNCTION "specifyLHSParam"
##
## Information necessary to translate standard uniform LHS sample into parameters of interest for paleo project 
###########

specifyLHSParam <- function(paramslist,name,type,lb,ub){
  newlist <- paramslist
  eval(parse(text=sprintf("newlist$%s <- list()",name)))
  eval(parse(text=sprintf("newlist$%s$type <- \"%s\"",name,type)))
  eval(parse(text=sprintf("newlist$%s$lb <- %s",name,lb)))
  eval(parse(text=sprintf("newlist$%s$ub <- %s",name,ub))) 	
  return(newlist)
}


###########
## FUNCTION "MakeLHSSamples"
##
## Samples from the uniform LHS and translates into desired parameter space
##  Returns a master data frame that will specify all MP files to be run
###########


###############################
#        SPECIFY PARAMETER RANGES
################################

# PER_SUITABLE=0.4,SNUGGLE=0.75,MAXDISPERSAL=500,BASELINE_DISPERSAL=0.05,
# MAXDISPERSAL_PLAGUE = 1000,PLAGUE_DISPERSAL=0.95, MAXDENS = 100,
# MINDENS = 15, BASELINE_MEANSURV = 0.6, BASELINE_PLAGUESURV=0.05,
# BASELINE_PLAGUESURV_RESIST=0.5,BASELINE_MEANFEC=3.2,
# FITNESS_COST=c(0.1,0.05),INITFREQ=c(0.09,0.1),DOMINANCE=dmat

MakeLHSSamples <- function(add){
  
  LHSParms <- list()    # initialize the container for parameter bounds
  
  ####  PER_SUITABLE
  LHSParms <- specifyLHSParam(LHSParms,"PER_SUITABLE",type="CONT",lb=0.1,ub=1)
  
  #### SNUGGLE
  LHSParms <- specifyLHSParam(LHSParms,"SNUGGLE",type="CONT",lb=0.01,ub=0.99)
  
  ### PLAGUE_DISPERSAL    0 to 100 individuals...  
  LHSParms <- specifyLHSParam(LHSParms,"PLAGUE_DISPERSAL",type="CONT",lb=0.05,ub=0.99)       
  
  ### MAXDENS
  LHSParms <- specifyLHSParam(LHSParms,"MAXDENS",type="INT",lb=25,ub=150)       
  
  ### BASELINE_PLAGUESURV
  LHSParms <- specifyLHSParam(LHSParms,"BASELINE_PLAGUESURV",type="CONT",lb=0.01,ub=0.1)
  
  ### BASELINE_MEANFEC
  LHSParms <- specifyLHSParam(LHSParms,"BASELINE_MEANFEC",type="CONT",lb=2.5,ub=6) 
  
  ### FITNESS_COST 
  LHSParms <- specifyLHSParam(LHSParms,"FITNESS_COST",type="CONT",lb=0.01,ub=0.25)
  
  ### INITFREQ
  LHSParms <- specifyLHSParam(LHSParms,"INITFREQ",type="CONT",lb=0.01,ub=0.25)
  
  #### DOMINANCE
  LHSParms <- specifyLHSParam(LHSParms,"DOMINANCE",type="CAT",lb=0,ub=4)    
  
  ##################
  ##### GENERATE LATIN HYPERCUBE SAMPLE
  
  nVars <- length(names(LHSParms))  
  
  LHS <- lhs::randomLHS(N_LHS_SAMPLES, nVars)   # generate multiple samples from parameter space according to a LHS sampling scheme
  
  temp <- as.data.frame(LHS)
  
  ### translate raw lhs samples into desired parameter space
  colnames(temp) <- names(LHSParms)
  parm=1
  for(parm in 1:nVars){
    if(LHSParms[[parm]]$type=="CONT"){
      temp[,parm] <- LHSParms[[parm]]$lb + LHS[,parm]*(LHSParms[[parm]]$ub-LHSParms[[parm]]$lb)
    }
    if(LHSParms[[parm]]$type=="CAT"){
      temp[,parm] <- ceiling(LHSParms[[parm]]$lb + LHS[,parm]*(LHSParms[[parm]]$ub-LHSParms[[parm]]$lb))
    }
    if(LHSParms[[parm]]$type=="INT"){
      temp[,parm] <- round(LHSParms[[parm]]$lb + LHS[,parm]*(LHSParms[[parm]]$ub-LHSParms[[parm]]$lb))
    }
  }
  
  if(add==FALSE) masterDF <- temp    #  storage container (data frame) to record relevant details for each MP file. Rows:MP file/LHS samples. Cols: relevant variables
  
  if(add==TRUE) masterDF <- rbind(masterDF,temp)
  
  setwd(dirs$DATA_DIR)
  ## name file for LHS parameters 
  write.csv(masterDF,sprintf("masterDF_prelim%s.csv",Sys.Date()),row.names=F)
  
  return(masterDF)
}


###########
## DEFINE USER PARAMETERS
###########

# This function defines a UserParams variable which is passed around and stores the user-defined parameters for the simulation. 
#   Note that many of these parameters can also be subject to global sensitivity analysis.

dmat <- matrix(c(1,1,0, 1,0,0), nrow=2,ncol=3,byrow = T)

DefineUserParams <- function(PER_SUITABLE=0.4,SNUGGLE=0.75,NFOCI=1,MAXDISPERSAL=500,BASELINE_DISPERSAL=0.05,
                             MAXDISPERSAL_PLAGUE = 1000,PLAGUE_DISPERSAL=0.95, MAXDENS = 100,
                             MINDENS = 15, BASELINE_MEANSURV = 0.6, BASELINE_PLAGUESURV=0.05,
                             BASELINE_PLAGUESURV_RESIST=0.5,BASELINE_MEANFEC=3.2,
                             FITNESS_COST=c(0.1,0.05),INITFREQ=c(0.09,0.1),DOMINANCE=dmat[[1]]){
    
  ############
  ## DEFINE RANGE OF POSSIBLE SCENARIOS  [not implemented yet]
  ############
  
  # Bg : "background"  - for looking at intermittent selection...
  
  # BgDiseaseSpatialExtents <- c("gridcell","subpatch","patch","regional","landscape")   
  # BgDiseaseFrequencies <- c(1,2,3,5,10,20,30,40,50)
  
  # ProbTransmissions <- seq(0,1,by=0.1)   # for now, no transmission
  
  
  UserParams <- list()
  
  #########
  # defining the landscape
  
  
  UserParams[["Landscape"]] <- list()
  
  # define the gridded landscape
  UserParams[["Landscape"]]$NROWS <- 50
  UserParams[["Landscape"]]$NCOLS <- 50
  
  UserParams[["Landscape"]]$CELLAREA_HA <- 1     # in hectares
  UserParams[["Landscape"]]$CELLAREA_M2 <- UserParams[["Landscape"]]$CELLAREA_HA*10000   # in square meters
  
  UserParams[["Landscape"]]$CELLWIDTH_M <- sqrt(UserParams[["Landscape"]]$CELLAREA_M2)    # in meters
  UserParams[["Landscape"]]$HALFCELLWIDTH_M <- UserParams[["Landscape"]]$CELLWIDTH_M/2
  
  
  # define the percent of the landscape that is suitable
  UserParams[["Landscape"]]$PER_SUITABLE <- PER_SUITABLE
  
  #########
  # defining the species biology:
  
  UserParams[["Dispersal"]] <- list()
  
  # dispersal params
  
  # define the clustering (degree to which species prefers to establish residence near members of its own kind)
  # 1 is complete tendency to cluster in space. 0 is agnostic to members of its own kind. -1 is tendency to avoid members of its own kind
  UserParams[["Dispersal"]]$SNUGGLE <- SNUGGLE
  
  # maximum annual dispersal distance (m)  
  UserParams[["Dispersal"]]$MAXDISPERSAL_M <- MAXDISPERSAL
  
  # number of colony foci to establish (within the dispersal range) after colonies get below the low-density threshold
  UserParams[["Dispersal"]]$NFOCI <- NFOCI
  
  # rate of transmission (per-disperser probability of initiating an outbreak in recipient population) (under optimal plague conditions?)
  #UserParams[["Dispersal"]]$PROB_TRANSMISSION <- 0  # [not implemented yet...] [maybe don't need...] [but critical if looking at tradeoffs]
  
  # dispersal rate (independent of density)
  UserParams[["Dispersal"]]$BASELINE_DISPERSAL <- BASELINE_DISPERSAL
  
  # dispersal distance (m) for plagued-out populations
  UserParams[["Dispersal"]]$MAXDISPERSAL_PLAGUE <-  MAXDISPERSAL_PLAGUE   # individuals from plagued-out populations might move farther than normal population
  
  # dispersal rate for plagued-out populations
  UserParams[["Dispersal"]]$PLAGUE_DISPERSAL <- PLAGUE_DISPERSAL     # individuals from plagued-out populations might have a higher tendency to move than normal populations- this can affect the spread of plague and the spread of plague resistance genes... 
  
  
  
  
  UserParams[["Popbio"]] <- list()
  
  # define the maximum per-cell number of individuals
  UserParams[["Popbio"]]$MAXDENS_HA <- MAXDENS
  UserParams[["Popbio"]]$MAXABUND <- UserParams[["Popbio"]]$MAXDENS_HA *UserParams[["Landscape"]]$CELLAREA_HA
  
  # define the minimum per-cell number of individuals (provide a simple hard Allee effect)
  UserParams[["Popbio"]]$MINDENS_HA <- MINDENS
  UserParams[["Popbio"]]$MINABUND <- UserParams[["Popbio"]]$MINDENS_HA*UserParams[["Landscape"]]$CELLAREA_HA
  
  # maximum survival under plague (limit to resistance)
  #   RESISTANCE_LIMIT <- 0.75   # deprecated  
  
  # baseline survival in a non_plague year (baseline survival for a naive population under no plague)
  UserParams[["Popbio"]]$BASELINE_MEANSURV <- BASELINE_MEANSURV
  
  # Variation in survival among years, expressed as CV
  UserParams[["Popbio"]]$CV_SURVIVAL <- 0.2 
  
  # minimum survival under plague (completely naive population)
  UserParams[["Popbio"]]$BASELINE_PLAGUESURV <- BASELINE_PLAGUESURV
  
  # survival under plague for resistant individuals
  UserParams[["Popbio"]]$BASELINE_PLAGUESURV_RESIST <- BASELINE_PLAGUESURV_RESIST
  
  # minimum survival for non-plague populations
  UserParams[["Popbio"]]$SURVMIN_NOPLAGUE <- 0.1
  
  # maximum survival for non-plague populations
  UserParams[["Popbio"]]$SURVMAX_NOPLAGUE <- 0.9
  
  # minimum survival for plague populations
  UserParams[["Popbio"]]$SURVMIN_PLAGUE <- 0.01
  
  # maximum survival for plague populations (circumscribing env stochasticity at the population)
  UserParams[["Popbio"]]$SURVMAX_PLAGUE <- 0.75
  
  
  # fecundity in a non-plague year (baseline, number of offspring per adult, not sex structured)
  UserParams[["Popbio"]]$BASELINE_MEANFEC <- BASELINE_MEANFEC
  
  # temporal variation in fecundity, expressed as a CV
  UserParams[["Popbio"]]$CV_FECUNDITY <- 0.5 
  
  
  
  # mean change in the ability to withstand plague among the survivor population. (accounts for limited heritability)
  #SURVCHANGE_NEXTPLAGUE <- 0.05  # deprecated
  
  # variation in the ability to survival plague among survivors, expressed as a SD [variation in standing genetic propensity to evolve resistance]
  #SD_SURVCHANGE_NEXTPLAGUE <- 0.05    # deprecated
  
  # plague survival rate as a function of the resistance factors?? [TODO] 
  
  # change in the ability to survive in a non-plague year, as a percentage of the fitness benefit in a plague year (fitness costs to resistance)
  #FITNESS_COST <- 0.1    # now means the degree to which survival is reduced for resistant individuals in nonplague years
  
  #########
  # defining the species genetics:
  
  UserParams[["Genetics"]] <- list()
  
  UserParams[["Genetics"]]$NGENES <- 2
  
  UserParams[["Genetics"]]$NWAYS_OF_GETTING <- 1   # ways of getting resistance
  UserParams[["Genetics"]]$RESISTANCE_SCENARIOS <- list()
  UserParams[["Genetics"]]$RESISTANCE_SCENARIOS[[1]] <- c("factor1","factor2")   # for now, needs all factors! 
  names(UserParams[["Genetics"]]$RESISTANCE_SCENARIOS[[1]]) <-c("AND","AND")   # names follow boolean conventions. In this case, both factors are required for resistance
  
  UserParams[["Genetics"]]$FITNESS_COST <- numeric(UserParams[["Genetics"]]$NGENES)
  
  UserParams[["Genetics"]]$FITNESS_COST[1] <- FITNESS_COST[1]     # fitness cost of the first gene   0.1
  UserParams[["Genetics"]]$FITNESS_COST[1] <- FITNESS_COST[2]     # fitness cost of the second gene   0.05
  
  UserParams[["Genetics"]]$INITFREQ <- numeric(UserParams[["Genetics"]]$NGENES)
  UserParams[["Genetics"]]$INITFREQ[1] <- INITFREQ[1]  #0.09
  UserParams[["Genetics"]]$INITFREQ[2] <- INITFREQ[2] #  0.1
  
  UserParams[["Genetics"]]$INITFREQ_SD <- 0.03     # degree of variation in initial frequency of resistance.
  
  temp <- matrix(0,nrow=UserParams[["Genetics"]]$NGENES,ncol=3)
  colnames(temp) <- c("2x(rr)","1x(rs)","0x(ss)")
  UserParams[["Genetics"]]$DOMINANCE <- temp
  UserParams[["Genetics"]]$DOMINANCE[1,] <- DOMINANCE[1,]   # dominant
  UserParams[["Genetics"]]$DOMINANCE[2,] <- DOMINANCE[2,]  # recessive
  
  return(UserParams)
}


############
## DO PLAGUE 
############
#
# This function determines which populations currently have plague given statistical model)

doPlague <- function(PlagueRaster=PlagueRaster,YearsSincePlague=YearsSincePlague,DensRaster=DensRaster,UserParams,PlagueModel,EnvCovs,fake,timestep){ 
  
  if(fake){
    nPlagueNeighbors <- raster::focal(PlagueRaster, w=matrix(1, nc=UserParams$Dispersal$MAXDISPERSAL_CELLS, nr=UserParams$Dispersal$MAXDISPERSAL_CELLS),na.rm=T)
    
    # raster::plot(PlagueRaster)
    # raster::plot(nPlagueNeighbors)
    
    newdf <- data.frame(
      dens = DensRaster@data@values,
      plaguepops = nPlagueNeighbors@data@values
    )
    
    # newdf <- data.frame(
    #   dens=0,
    #   plaguepops=0
    # )
    
    prob <- plogis(as.numeric(predict(PlagueModel,newdata=newdf)))
    ProbRaster <- raster::setValues(PlagueRaster,values=prob)
    ndx <- !is.na(prob)
    prediction <- prob
    prediction[ndx] <- rbinom(length(which(ndx)),1,prediction[ndx])
    
  }else{     ### IF we are using a real plague model!!
    
    #### PATCH SIZE
    
    #DensRaster2 <- raster::crop(DensRaster,extent(EnvCovs$lat.c))   ## cropped to study landscape
    
    # plot(DensRaster)
    # plot(EnvCovs$lat.c)
    
    temp.patch <- raster::clump(DensRaster, directions=8, gaps=FALSE)   # find occupied patches
    cells.patch <- raster::freq(temp.patch)     # number of patches in clump
    cells.patch <- as.data.frame(cells.patch[1:(nrow(cells.patch)-1),,drop=F])  # if no patches, this will have zero rows
    
    if(nrow(cells.patch)>0){
      PATCHSIZE <- raster::subs(temp.patch, cells.patch, by=1, which=2)  # number of cells in patch, for each pixel
    } else{
      PATCHSIZE <- raster::reclassify(DensRaster,rcl=c(-Inf,Inf,0))
    }
    
   # plot(PATCHSIZE)
    
    #### YEARS SINCE LAST PLAGUE (now computed in the main loop)
    
    
    #### SLOPE COST TO PLAGUE
    
    slope.trans <- gdistance::transition(EnvCovs$slope.c, mean, 8)   # cost to move to neighboring cells (move to initialization?)
    
    dist.noplague <- DensRaster   # set up a baseline "dist to plague" raster
    raster::values(dist.noplague) <- 50000    # 500000? #Need a raster for distance if there's no plague in the area (put distance at 500km)
    
    temp <- PlagueRaster  
    raster::values(temp) <- ifelse(raster::values(temp)==0, NA, raster::values(temp))
    if(any(!is.na(raster::values(temp)))){  # if any plagued cells?
      plague.pts <- raster::rasterToPoints(PlagueRaster, fun=function(x){x==1})[,1:2]  ##Need plagued areas as points
      slope.cost <- gdistance::accCost(slope.trans, plague.pts)
      slope.cost <- raster::resample(slope.cost, DensRaster)
    } else {
      slope.cost <- dist.noplague    # if all no plague, then dist to plague is maximized for all cells 
    }
    
    ################
    # CONSTRUCT DATA FRAME OF PREDICTOR VARIABLES
    ################
    
    
    newdf <- data.frame(
      psize = raster::values(PATCHSIZE),
      slope.cost.plague = raster::values(slope.cost),
      winspr.prcp.2prev = raster::values(EnvCovs$prcp.WinSpr[[max(1,(timestep-2))]]),
      winspr.prcp.prev = raster::values(EnvCovs$prcp.WinSpr[[max(1,(timestep-1))]]),
      winspr.prcp = raster::values(EnvCovs$prcp.WinSpr[[max(1,timestep)]]),
      sumfal.prcp.2prev = raster::values(EnvCovs$prcp.SumFal[[max(1,(timestep-2))]]),
      sumfal.prcp.prev =  raster::values(EnvCovs$prcp.SumFal[[max(1,(timestep-1))]]), 
      sumfal.prcp = raster::values(EnvCovs$prcp.SumFal[[max(1,timestep)]]), 
      year.prcp = raster::values(EnvCovs$prcp.year[[max(1,timestep)]]), 
      sand0 = raster::values(EnvCovs$sand0.c), 
      tmax.prev = raster::values(EnvCovs$tmax[[max(1,(timestep-1))]]), 
      tmax = raster::values(EnvCovs$tmax[[max(1,timestep)]]),  
      elev = raster::values(EnvCovs$NED.c), 
      lat = raster::values(EnvCovs$lat.c), 
      long = raster::values(EnvCovs$long.c), 
      years.plague = raster::values(YearsSincePlague)
    )
    
    keep <- !(raster::values(DensRaster)%in%c(NA,0))    # in an occupied patch
    #rmv <- is.na(raster::values(DensRaster))   # out of a patch    # Q: was the model fitted to only data within patches?
    
    #newdf[rmv,] <- NA   # set  
    
    newdf2 <- newdf[keep,]
    
    head(newdf2)
    
    ################
    # MAKE THE PREDICTION!
    ################
    
    
    
    temp <- predict(PlagueModel,newdata=newdf2,type="prob")     # this can take a while unfortunately, but only has to be done once per time step
    
    prob <- numeric(nrow(newdf2))
    i=1
    for(i in 1:nrow(newdf2)){
      prob[i] <- temp[[i]][,"plague.1"][1]
    }
    
    
    #ProbRaster <- raster::setValues(PlagueRaster,values=prob)
    #ndx <- !is.na(prob)
    
    #prediction <- prob
    prediction2 <- rbinom(length(prob),1,prob)    # maybe turn this into patch-level phenomenon?
    
    prediction <- rep(0,times=nrow(newdf))
    prediction[keep] <- prediction2
    
  }
  
  PlagueRaster <- raster::setValues(PlagueRaster,values=prediction)   # set the plagueraster according to the statistical model
  return(PlagueRaster)
}

############
## GET PLAGUE MODEL
############
#
#  This function (for now) builds a model of where and when plague occurs on the landscape as a function of the density of 
#   colonies etc. 

GetPlagueModel <- function(fake=TRUE){
  if(fake){
    BETADENS <- 0.02  # 0.005
    BETAPLAGUE <- 0.01  # 0.9
    INTERACTION <- 0.05  # 0.03
    INTERCEPT <- -6  # -6
    
    faken <- 1000
    fakedens <- seq(0,200,length=faken)   # fakedens <- rep(seq(0,200,by=50),each=10)
    fakeplaguepops <- sample(c(0:5),faken,replace=T)   # fakeplaguepops <- rep(c(0:9),times = 5)
    fakeplagueprob <- plogis(INTERCEPT + BETADENS*fakedens + BETAPLAGUE*fakeplaguepops + INTERACTION*fakeplaguepops*fakedens)
    # matrix(fakeplagueprob,ncol=10,byrow=T)
    
    # raster::plot(fakeplagueprob~fakeplaguepops)
    #raster::plot(fakeplagueprob~fakedens)
    fakeplague <- rbinom(faken,1,fakeplagueprob)
    
    dataFrame <- data.frame(plague=fakeplague,dens=fakedens,plaguepops=fakeplaguepops)
    
    plaguemodel <- glm(plague~1+plaguepops*dens,data=dataFrame,family="binomial")
    
    #plogis(predict(plaguemodel,newdata=data.frame(dens=c(10:10),plaguepops=c(1:10))))
    
    #summary(plaguemodel)
    
  }else{   # otherwise use real plague model!
    setwd(dirs$plaguemod$DataDir)
    load("randomForestModel_2017-05-12.RData")    # load pre-constructed plague model (random forest!)
    plaguemodel <- rfp1
  }
  
  return(plaguemodel)
}


############
## INITIALIZE DISPERSAL STRUCTURES
############
#
#  This function creates several structures useful for dispersal modeling

InitializeDispersal <- function(UserParams){   # env
  
  #UserParams <- get("UserParams",envir=env)
  
  UserParams$Dispersal$MAXDISPERSAL_CELLS <- floor(UserParams$Dispersal$MAXDISPERSAL_M/UserParams$Landscape$CELLWIDTH_M)
  UserParams$Dispersal$MAXDISPERSAL_CELLS_PLAGUE <- floor(UserParams$Dispersal$MAXDISPERSAL_PLAGUE/UserParams$Landscape$CELLWIDTH_M)
  UserParams$Dispersal$MAXDISPERSAL_M <- UserParams$Dispersal$MAXDISPERSAL_CELLS*UserParams$Landscape$CELLWIDTH_M
  UserParams$Dispersal$MAXDISPERSAL_PLAGUE <- UserParams$Dispersal$MAXDISPERSAL_CELLS_PLAGUE*UserParams$Landscape$CELLWIDTH_M
  
  DispMask <- list()
  DispKernel <- list()
  DispRates <- list()
  i="noPlague"
  for(i in c("noPlague","plague")){
    if(i=="noPlague"){ DispRates$maxdispcells = UserParams$Dispersal$MAXDISPERSAL_CELLS ; DispRates$maxdisp = UserParams$Dispersal$MAXDISPERSAL_M} 
    if(i=="plague") {DispRates$maxdispcells = UserParams$Dispersal$MAXDISPERSAL_CELLS_PLAGUE ; DispRates$maxdisp = UserParams$Dispersal$MAXDISPERSAL_PLAGUE}
    tempDispMask <- matrix(1,nrow=((DispRates$maxdispcells*2)+1),ncol=((DispRates$maxdispcells*2)+1))     # dispersal mask defines where it is possible to go...
    #DispMask_Plague <- matrix(1,nrow=((MAXDISPERSAL_CELLS_PLAGUE*2)+1),ncol=((MAXDISPERSAL_CELLS_PLAGUE*2)+1))     # dispersal mask defines where it is possible to go...
    
    # exclude values greater than the max dispersal distance from the dispersal kernel 
    xs <- seq(UserParams$Landscape$HALFCELLWIDTH_M,DispRates$maxdisp*2+UserParams$Landscape$HALFCELLWIDTH_M,by=UserParams$Landscape$CELLWIDTH_M)
    ys <- seq(UserParams$Landscape$HALFCELLWIDTH_M,DispRates$maxdisp*2+UserParams$Landscape$HALFCELLWIDTH_M,by=UserParams$Landscape$CELLWIDTH_M)
    
    focalcell <- ((length(xs)-1)/2)+1
    focalx <- xs[focalcell]
    focaly <- ys[focalcell]
    col=1;row=1
    for(col in 1:length(xs)){
      for(row in 1:length(ys)){
        dist <- sqrt((xs[col]-focalx)^2 + (ys[row]-focaly)^2)
        if(dist>DispRates$maxdisp) tempDispMask[row,col] = NA 
      }
    }
    tempDispKernel <- tempDispMask   # for now, dispersal kernel is uniform    ; dispersal kernel defines where pdogs go if they lose a home colony
    #DispKernel_Plague <- DispMask_Plague
    DispMask[[i]] <- tempDispMask
    DispKernel[[i]] <- tempDispKernel
    
  }
  
  # set up structures for colony expansion- individuals move to neighboring cells first..   [note: for now, donuts obey only the non-plague max dispersal]
  
  donuts <- list()
  j="noPlague"
  for(j in c("noPlague","plague")){
    if(j=="noPlague"){ DispRates$maxdispcells = UserParams$Dispersal$MAXDISPERSAL_CELLS ; DispRates$maxdisp = UserParams$Dispersal$MAXDISPERSAL_M} 
    if(j=="plague") {DispRates$maxdispcells = UserParams$Dispersal$MAXDISPERSAL_CELLS_PLAGUE ; DispRates$maxdisp = UserParams$Dispersal$MAXDISPERSAL_PLAGUE}
    donuts[[j]] <- list()
    donut_template <- matrix(0,nrow=((DispRates$maxdispcells*2)+1),ncol=((DispRates$maxdispcells*2)+1),byrow=T)
    newseq <- c(-DispRates$maxdispcells:DispRates$maxdispcells)
    
    ndx <- which(abs(newseq)<1)
    donutsZero <- donut_template
    donutsZero[ndx,ndx] <- 1 
    i=1
    for(i in 1:DispRates$maxdispcells){
      ndx <- which(abs(newseq)<(i+1))
      donuts[[j]][[i]] <- donut_template
      donuts[[j]][[i]][ndx,ndx] <- 1
      donuts[[j]][[i]] <- donuts[[j]][[i]] - donutsZero
      #donuts[[i]] <- apply(donuts[[i]],c(1,2),function(t) ifelse(t==0,NA,t)) 
      donutsZero <- donuts[[j]][[i]] + donutsZero
    }
  }
  
  donuts[[j]] <- lapply(donuts[[j]],function(k) apply(k,c(1,2),function(t) ifelse(t==0,NA,t)))
  
  #UserParams <- UserParams   # save to global environment
  
  templist <- list()
  templist$DispList <- list()
  templist$DispList$donuts <- donuts
  templist$DispList$DispMask <- DispMask
  templist$DispList$DispKernel <- DispKernel
  templist$UserParams <- UserParams
  return(templist)
  
  # assign(x="donuts",value=donuts, envir = env)
  # assign(x="DispMask",value=DispMask, envir = env)
  # assign(x="DispKernel",value=DispKernel, envir = env)
  # assign(x="UserParams",value=UserParams, envir = env)
  
  # DispList$donuts[['noPlague']][[4]]
}



############
## INITIALIZE LANDSCAPE (patches of habitat)
############
#
#  This function generates a fake landscape in which to model plague dynamics

InitializeLandscape <- function(solid=F,fake=F,UserParams){   #env
  
  maxdisp <- max(UserParams$Dispersal$MAXDISPERSAL_CELLS,UserParams$Dispersal$MAXDISPERSAL_CELLS_PLAGUE)  # for extending the landscape to accommodate dispersal
  
  if(fake){
    templateRaster <- raster::raster(nrows=UserParams$Landscape$NROWS, ncols=UserParams$Landscape$NCOLS, xmn=0, 
                              xmx=UserParams$Landscape$CELLWIDTH_M*UserParams$Landscape$NROWS,ymn=0, 
                              ymx=UserParams$Landscape$CELLWIDTH_M*UserParams$Landscape$NCOLS,vals=NA)    # template raster
    # raster::plot(templateRaster)
    
    if(solid){
      patchRaster <- raster::setValues(templateRaster,1) 
    }else{
    
    # use utility function from secr package to initialize landscape... 
      tempgrid <- secr::make.grid(nx = UserParams$Landscape$NCOLS, ny = UserParams$Landscape$NROWS, spacing = UserParams$Landscape$CELLWIDTH_M,
                            detector = "single", originxy = c(0,0), hollow = F,
                            ID = "alphay")
      
      # secr::plot(tempgrid)  
      
      tempmask <- secr::make.mask(traps=tempgrid, buffer = UserParams$Landscape$HALFCELLWIDTH_M, spacing = UserParams$Landscape$CELLWIDTH_M, 
                            nx = UserParams$Landscape$NCOLS, ny = UserParams$Landscape$NROWS, type =
                              c("traprect"))
      
      # plot(tempmask)
      
      temppatches <- secr::randomHabitat(mask=tempmask, p = 0.4, A = UserParams$Landscape$PER_SUITABLE, directions = 4, minpatch = 20,
                                   drop = FALSE, covname = "habitat", plt = FALSE)
      
      #BaseLandscape$patchRaster <- templateRaster
      #patchvals <- as.vector(t(as.matrix(covariates(temppatches)$habitat)))
      patchRaster <- raster::setValues(templateRaster,values=secr::covariates(temppatches)$habitat)
      
      patchRaster <- raster::reclassify(patchRaster,rcl=c(-Inf,0.5,NA, 0.6,Inf,1))   # raster of habitat patches
      # plot(patchRaster) 
    }
    ENV_COVARS <- NULL
    
  }else{   # read in real landscape: patches and covariates.
    
    ##########
    # LOAD ENVIRONMENTAL COVARIATES  (these do not change!)
    ##########
    
    setwd(dirs$plaguemod$CovDir)
    load("EnvCov_smallPawnee.RData")
    
      ##### DEFINE ENV COVARIATES
       ## note that env covariates have the right projection and resolution. 
      
      plot=F
      
      if(plot){
        raster::plot(ENV_COVARS$lat.c)
        raster::plot(ENV_COVARS$long.c)
        raster::plot(ENV_COVARS$NED.c)     ## elevation
        raster::plot(ENV_COVARS$prcp.SumFal)   # time series of maps
        raster::plot(ENV_COVARS$prcp.WinSpr)   # time series
        raster::plot(ENV_COVARS$prcp.year)     # ts
        raster::plot(ENV_COVARS$sand0.c)     #??
        raster::plot(ENV_COVARS$sand2.c)
        raster::plot(ENV_COVARS$slope.c)
        raster::plot(ENV_COVARS$tmax)   # ts  
      } 
      
      ##################
      # IF NEEDED ASSEMBLE THE COVARIATES
      ##################
      
      convert=F 
      if(convert){
        setwd(dirs$plaguemod$CovDir)
        
        load("EnvCov_stacks_smallPawnee_FUTURE.RData")
        load("EnvCov_stacks_smallPawnee_PAST.RData")
        load("EnvCov_stacks_smallPawnee.RData")
        
        ENV_COVARS <- list()
        
        ENV_COVARS$lat.c <- lat.c
        ENV_COVARS$long.c <- long.c
        ENV_COVARS$NED.c <- NED.c      # elevation 
        ENV_COVARS$prcp.SumFal <- raster::brick(c(prcp.SumFal.past,prcp.SumFal,prcp.SumFal.f.r))
        allyears <- as.numeric(unlist(regmatches(names(ENV_COVARS$prcp.SumFal), gregexpr("[[:digit:]]+", names(ENV_COVARS$prcp.SumFal)))))
        names(ENV_COVARS$prcp.SumFal) <- allyears
        
        ENV_COVARS$prcp.WinSpr <- raster::brick(c(prcp.WinSpr.past,prcp.WinSpr,prcp.WinSpr.f.r))
        names(ENV_COVARS$prcp.WinSpr) <- allyears
        
        ENV_COVARS$prcp.year <- raster::brick(c(prcp.year.past,prcp.year,prcp.year.f.r))
        names(ENV_COVARS$prcp.year) <- allyears
        
        ENV_COVARS$sand0.c <- sand0.c
        ENV_COVARS$sand2.c <- sand2.c
        ENV_COVARS$slope.c <- slope.c
        
        ENV_COVARS$tmax <- raster::brick(c(tmax.past,tmax,tmax.f.r))
        names(ENV_COVARS$tmax) <- allyears
        
        setwd(dirs$plaguemod$CovDir)
        save(ENV_COVARS, file="EnvCov_smallPawnee.RData")   # was EnvCov_stacks_smallPawnee.RData
        
        rm(lat.c,long.c,NED.c,
           prcp.SumFal.past,prcp.SumFal,prcp.SumFal.f.r,
           prcp.WinSpr.past,prcp.WinSpr,prcp.WinSpr.f.r,
           prcp.year.past,prcp.year,prcp.year.f.r,
           sand0.c,sand2.c,slope.c,
           tmax.past,tmax,tmax.f.r)
      }
      
      
      ###############
      # SET UP TEMPLATE RASTER
      ###############
      
      ## create template
      templateRaster <- raster::raster()   # template to create all other rasters for landscape of interest
      raster::extent(templateRaster) <- raster::extent(ENV_COVARS$lat.c) #Using smaller extent from covariates (snippet of Pawnee)
      raster::res(templateRaster) <- raster::res(ENV_COVARS$lat.c)   # set the correct resolution
      
      ###################
      # EXTEND THE ENVIRONMENTAL COVARIATES
      
      ENV_COVARS$lat.c <- raster::extend(ENV_COVARS$lat.c,maxdisp,value=NA)     #  raster::extend(patchRaster,maxdisp,value=NA)
      ENV_COVARS$long.c <- raster::extend(ENV_COVARS$long.c,maxdisp,value=NA)
      ENV_COVARS$NED.c <- raster::extend(ENV_COVARS$NED.c,maxdisp,value=NA)
      ENV_COVARS$prcp.SumFal <- raster::extend(ENV_COVARS$prcp.SumFal,maxdisp,value=NA)
      ENV_COVARS$prcp.WinSpr <- raster::extend(ENV_COVARS$prcp.WinSpr,maxdisp,value=NA)
      ENV_COVARS$prcp.year <- raster::extend(ENV_COVARS$prcp.year,maxdisp,value=NA)
      ENV_COVARS$sand0.c <- raster::extend(ENV_COVARS$sand0.c,maxdisp,value=NA)
      ENV_COVARS$sand2.c <- raster::extend(ENV_COVARS$sand2.c,maxdisp,value=NA)
      ENV_COVARS$slope.c <- raster::extend(ENV_COVARS$slope.c,maxdisp,value=NA)
      ENV_COVARS$tmax <- raster::extend(ENV_COVARS$tmax,maxdisp,value=NA)
      
      ##############
      # SET UP PATCH RASTER
      ##############
      
      setwd(dirs$plaguemod$DataDir)
      
      #Read in the "master colony" shapefile
      
      load("PawneeColonies.RData")
      
      #load("MasterColonies_2017-05-09.RData")    # elizabeth's version had way too much stuff in it!
      #save(master.colony.Pawne,file="PawneeColonies.RData")
      
      
      plot=F
      
      if(plot){
        plot(ENV_COVARS$NED.c)
        plot(master.colony.Pawne,add=T)
      }
      
      #Crop master colony polygon to working extent
      mcP <- raster::crop(master.colony.Pawne, raster::extent(templateRaster))   # do I need this?
      
      patchRaster <- raster::rasterize(mcP, templateRaster, background=NA, field=1)
      
      
      #############
      # UPDATE USER PARAMS
      #############
      
      UserParams$Landscape$NROWS <- raster::nrow(templateRaster)
      UserParams$Landscape$NCOLS <- raster::ncol(templateRaster)
      UserParams$Landscape$CELLAREA_HA <- prod(raster::res(templateRaster))/10000
      UserParams$Landscape$CELLAREA_M2 <- prod(raster::res(templateRaster))
      UserParams$Landscape$CELLWIDTH_M <- raster::res(templateRaster)[1]
      UserParams$Landscape$HALFCELLWIDTH_M <- UserParams$Landscape$CELLWIDTH_M / 2
      UserParams$Landscape$PER_SUITABLE <- length(which(raster::values(patchRaster)==1))/raster::ncell(templateRaster)    # overwrite the percent suitable if a real landscape is used (this is constant!)
    
  }  

    # extend patch raster to go outside the landscape bounds to the max dispersal distance...
  maxdisp <- max(UserParams$Dispersal$MAXDISPERSAL_CELLS,UserParams$Dispersal$MAXDISPERSAL_CELLS_PLAGUE)
  patchRaster <- raster::extend(patchRaster,maxdisp,value=NA)
  
  #patchMatrix <- as.matrix(newraster)    # matrix of habitat patches
  
  KRaster <- patchRaster * UserParams$Popbio$MAXDENS_HA     # matrix of carrying capacity
  # plot(KRaster)
  
  patchIDRaster <- raster::clump(patchRaster,directions=4,gaps=F)   # determine unique ID for each patch... 
  # plot(patchIDRaster)
  
  nPatches <- raster::cellStats(patchIDRaster,"max")   # number of patches in the landscape
  
  nCells <- raster::ncell(patchRaster)
  
  # Data frame of coordinates for all non-na cells
  # focalCells <- which(!is.na(patchRaster@data@values))
  # xy_df <- data.frame(
  #   raster::xyFromCell(patchRaster,focalCells)
  # )
  # 
  # head(xy_df)
  
  # extent- in terms of cell centroids... (used for dispersal)
  UserParams$Landscape$MINX <- patchRaster@extent@xmin+UserParams$Landscape$HALFCELLWIDTH_M
  UserParams$Landscape$MAXX <- patchRaster@extent@xmax-UserParams$Landscape$HALFCELLWIDTH_M
  UserParams$Landscape$MINY <- patchRaster@extent@ymin+UserParams$Landscape$HALFCELLWIDTH_M
  UserParams$Landscape$MAXY <- patchRaster@extent@ymax-UserParams$Landscape$HALFCELLWIDTH_M
  
  UserParams$Landscape$FULLEXTENT <- raster::extent(UserParams$Landscape$MINX,UserParams$Landscape$MAXX,UserParams$Landscape$MINY,UserParams$Landscape$MAXY)
  #UserParams <- UserParams  # save to global env
  
  BaseLandscape <- list(patchRaster = patchRaster,
                        KRaster = KRaster,
                        patchIDRaster = patchIDRaster,
                        nPatches = nPatches,
                        nCells = nCells
                     )
  
  returnlist <- list()
  returnlist$BaseLandscape <- BaseLandscape
  returnlist$UserParams <- UserParams
  returnlist$EnvCovs <- ENV_COVARS
  
  # assign(x="patchRaster",value=patchRaster, envir = env)
  # assign(x="KRaster",value=KRaster, envir = env)
  # assign(x="patchIDRaster",value=patchIDRaster, envir = env)
  # assign(x="nPatches",value=nPatches, envir = env)
  # assign(x="nCells",value=nCells, envir = env)
  # assign("BaseLandscape",value=BaseLandscape, envir=env)
  # assign(x="UserParams",value=UserParams, envir = env)
  return(returnlist)
}



##################
# DISPERSAL FUNCTIONS
##################
#
#  NOTE: there is a baseline dispersal, plus any individuals that are crowded out... 
#    individuals that are crowded out are dictated by the snuggle function
#    individuals that are baseline dispersers are governed by the kernel


##################
# UpdateAbund: UPDATE ABUNDANCE TO ACCOUNT FOR DISPERSAL
##################
#
#  This function takes a central point and distributes individuals to other cells as dictated by a dataframe, "toAdd_df"
#    This function also updates the allele frequency layers appropriately
#

# NOTE: need to deal with boundary conditions... 

###################
# Add (or subtract) individuals from cells during dispersal phase

UpdateAbund <- function(focalxy=xy, stack, df=toAdd_df,DensRaster,FreqList){
  newstack <- stack
  focalndx <- raster::cellFromXY(DensRaster,focalxy)
  ndx <- raster::cellFromXY(DensRaster,df[,c("x","y")])
  indiv <- df[,3]   # individuals to add
  names <- names(stack)
  i=names[2]
  for(i in names){
    if(i=="DensRaster"){
      Add <- indiv
    }else{
      Add <- eval(parse(text=sprintf("FreqList[[\"%s\"]]",i)))[focalndx]*indiv   # otherwise, add to the total in each cell for averaging (evolving!)
      # add the focal value of evolving layers to the surrounding cells   [change this!!!]
    }
    newstack[[i]][ndx] <- newstack[[i]][ndx]+Add
  }
  return(newstack)
}

##################
# UpdateStack: UPDATE ALLELE FREQUENCIES 
##################
#
#  This function computes the mean allele frequencies and saves to the global environment 

UpdateStack <- function(stack){  #env
  names <- names(stack)
  # i=names[2]
  keep <- c()
  for(i in names){
    if(i!="DensRaster"){
      ndx <- which(stack[["DensRaster"]]@data@values>0.00001)  #only consider cells with actual individuals
      stack[[i]][ndx] <- stack[[i]][ndx]/stack[['DensRaster']][ndx]    # get the population mean for key evolving variables
      stack[[i]][stack[["DensRaster"]]@data@values<=0.00001] <- 0
      keep <- c(keep,i)
    }
  }
  return <- raster::subset(stack,subset=keep,drop=FALSE)
  #assign(x="FreqList",value=return, envir = env)   # assign the variable to the global environment
  return(return)
}


##################
# GetDispersalRates: helper function for getting dispersal rates 
##################
#
#  This function determines the dispersal rates from this cell (depends on plague status)
##########


GetDispersalRates <- function(plagueStatus="noPlague",UserParams){
  if(plagueStatus=="noPlague") {maxdisp<-UserParams$Dispersal$MAXDISPERSAL_M;disprate<-UserParams$Dispersal$BASELINE_DISPERSAL;maxdispcells<-UserParams$Dispersal$MAXDISPERSAL_CELLS}
  if(plagueStatus=="plague") {maxdisp<-UserParams$Dispersal$MAXDISPERSAL_PLAGUE;disprate<-UserParams$Dispersal$PLAGUE_DISPERSAL;maxdispcells<-UserParams$Dispersal$MAXDISPERSAL_CELLS_PLAGUE}
  ret <- list()
  ret$maxdisp <- maxdisp
  ret$disprate <- disprate
  ret$maxdispcells <- maxdispcells
  # assign(x="maxdisp",value=maxdisp, envir = env) 
  # assign(x="disprate",value=disprate, envir = env) 
  # assign(x="maxdispcells",value=maxdispcells, envir = env) 
  return(ret)
}



#############
# DEFINE THE NEIGHBORHOOD TO WHICH PDOGS COULD MOVE
#############

makeNeighborhoodRasters <- function(newStack,plagueStatus,xy,UserParams,DispList){
  mask=DispList$DispMask[[plagueStatus]]
  kernel=DispList$DispKernel[[plagueStatus]]
  DispRates <- GetDispersalRates(plagueStatus,UserParams) # already done... 
  neighborhood.extent <- raster::extent(xy$x-DispRates$maxdisp-UserParams$Landscape$HALFCELLWIDTH_M,
                                xy$x+DispRates$maxdisp+UserParams$Landscape$HALFCELLWIDTH_M,
                                xy$y-DispRates$maxdisp-UserParams$Landscape$HALFCELLWIDTH_M,
                                xy$y+DispRates$maxdisp+UserParams$Landscape$HALFCELLWIDTH_M)

  neighborhood_raster <- raster::crop(newStack[["DensRaster"]],neighborhood.extent)    # possible cells to move to and from [TODO: change this to current status, not former]
  #neighborhood_raster <- neighborhood_raster - (neighborhood_raster*disprate) # make sure that there is space for incoming dispersers [revisit this]
  
  ################
  # Make a neighborhood mask
  
  vals <- as.vector(t(mask))  # convert disp mask to a vector format
  #vals[vals==0] <- NA
  #neighborhood_kernel_raster <- neighborhood_raster 
  length(vals)
  raster::ncell(neighborhood_raster)
  neighborhood_mask_raster <- raster::setValues(neighborhood_raster,values=vals)  # NOTE: this could be a dispersal kernel
  
  ###############
  # Make a neighborhood kernel
  vals <- as.vector(t(kernel))
  neighborhood_kernel_raster <- raster::setValues(neighborhood_raster,values=vals)
  
  #     plot(neighborhood_raster)
  #     plot(neighborhood_mask_raster)
  neighborhood_raster <- neighborhood_raster * neighborhood_mask_raster  # abundance in real neighborhood (assuming edge is hard barrier)  
  
  neighborhood_mask_raster <- raster::reclassify(neighborhood_raster,c(-Inf,Inf,1))
  
  neighborhood_kernel_raster <- neighborhood_kernel_raster*neighborhood_mask_raster
  
  # names(neighborhoodRasters)
  freeSpace_raster <- raster::reclassify(UserParams$Popbio$MAXDENS_HA-neighborhood_raster,rcl=c(-Inf,0,0))  # indicator of how much "free space" is available for cells in the neighborhood
  # plot(freeSpace_raster)
  
  neighborhood_stack <- raster::stack(list(
    maskRaster = neighborhood_mask_raster,
    kernelRaster = neighborhood_kernel_raster,
    abundRaster = neighborhood_raster,
    freeSpaceRaster = freeSpace_raster
  ))
  
  return(neighborhood_stack)
}


####################
# FUNCTION FOR SPREADING OUT THE POPULATION
####################
#account for the snuggle effect...

#############################
# DISPERSAL MODE #1: CLUSTER WITH EXISTING POPULATION
#############################

#names(neighborhoodRasters)

SpreadOut <- function(newStack,DensRaster,FreqList,leave_snuggle,
                      plagueStatus,
                      xy,
                      neighborhood_raster,
                      neighborhood_mask_raster,
                      freeSpace_raster,
                      UserParams,
                      DispList){  # env
  DispRates <- GetDispersalRates(plagueStatus,UserParams) 
  updatedStack <- newStack
  i=2
  #toAdd <- data.frame(x=numeric(0),y=numeric(0))
  for(i in 1:DispRates$maxdispcells){
    donutraster <- raster::setValues(neighborhood_raster,values=DispList$donuts[[plagueStatus]][[i]])
    # plot(donutraster)
    #tempmask <- raster::reclassify(neighborhood_raster,rcl=c(-Inf,Inf,1))  # NA,NA,0, 
    donutraster <- donutraster*neighborhood_mask_raster  # all non-NA cells within the donut
    #         plot(donutraster+tempmask)
    #         plot(donutraster)
    #         plot(tempmask)
    #         plot(neighborhood_mask_raster)
    #         plot(donutraster==tempmask)
    #         plot(freeSpace_raster)
    temp <- freeSpace_raster*donutraster
    #plot(raster::reclassify(neighborhood_raster,rcl=c(-Inf,Inf,1)))
    #temp <- raster::reclassify(MAXABUND-temp,rcl=c(-Inf,0,0))   # REVISIT!
    # plot(temp)
    freeSpace <- temp@data@values[which(temp@data@values>0)]
    # if there is free space... then put individuals in that space... 
    if(length(freeSpace)>0){
      toAdd_df <- as.data.frame(raster::xyFromCell(temp,which(temp@data@values>0))) # cells with empty space to go    
      # allocate those individuals 
      #toAdd_df$indiv <- 0
      toAdd <- min(sum(freeSpace),leave_snuggle)
      toAdd_df$indiv <- rmultinom(1,toAdd,prob=freeSpace)[,1]    # disperser individuals to add to each cell
      dispndx <- which(toAdd_df$indiv>0)
      updatedStack <- UpdateAbund(focalxy=xy,stack=updatedStack,df=toAdd_df[dispndx,],DensRaster,FreqList)  # add these disperser individuals and update the evolving factors
      leave_snuggle <- leave_snuggle-toAdd   # remove these individuals from the "leave" pool
      #ndx <- raster::cellFromXY(newPop,toAdd_df[,c("x","y")])
      #newPop[ndx] <- newPop[ndx]+toAdd_df$indiv
      # plot(newPop)
      # plot(updatedStack[["PopArray"]])
      
      
    } # end if free space available 
    
    # if all individuals have been dispersed, then stop the loop
    if(leave_snuggle<=0) break
    #if(i==MAXDISPERSAL_CELLS) # don't need to kill off the existing individuals explicitly... they automatically don't get propagated... 
  } #end loop through donuts
  return(updatedStack)
}


####################
# FUNCTION FOR PERFORMING LONG-DISTANCE DISPERSAL
####################
#account for the snuggle effect...

#############################
# DISPERSAL MODE #2: LONG-DISTANCE (not density dependent...)  [mostly redundant with above, move into function!]  [these dispersers, unlike the snuggle dispersers, are forced to integrate with the recipient populations]
############################# 

#names(neighborhoodRasters)
LongDistanceDispersal <- function(newStack,DensRaster,FreqList,
                                  leave_kernel=leave_kernel,
                                  plagueStatus=plagueStatus,
                                  xy=xy,
                                  neighborhood_raster=neighborhoodRasters[['abundRaster']],
                                  freeSpace_raster=neighborhoodRasters[['freeSpaceRaster']],
                                  neighborhood_kernel_raster=neighborhoodRasters[['kernelRaster']],
                                  UserParams){
  
  # disperse via the kernel...
  #leave_kernel
  #plot(neighborhood_raster)
  #plot(raster::reclassify(neighborhood_raster,rcl=c(-Inf,Inf,1)))
  
  DispRates <- GetDispersalRates(plagueStatus,UserParams)
  
  freeSpace <-  freeSpace_raster@data@values[which(freeSpace_raster@data@values>0)]
  if(length(freeSpace)>0){
    freeSpace <- (freeSpace+1)/max(freeSpace)  # convert to weighting factor in line with kernel
    toAdd_df <- as.data.frame(raster::xyFromCell(freeSpace_raster,which(freeSpace_raster@data@values>0))) # cells with empty space and habitat to go   
    kernel.weights <- neighborhood_kernel_raster[raster::cellFromXY(neighborhood_kernel_raster,toAdd_df[,c(1,2)])] 
    # allocate those individuals 
    #toAdd_df$indiv <- 0
    #toAdd <- min(sum(freeSpace),leave_kernel)
    toAdd_df$indiv <- rmultinom(1,leave_kernel,prob=freeSpace*kernel.weights)[,1]    # disperser individuals to add to each cell (favor cells that have more free space...)
  }else{  # if there is no free space, then force dispersal anyway
    toAdd_df <- as.data.frame(raster::xyFromCell(neighborhood_raster,which(!is.na(neighborhood_raster@data@values)))) # cells with habitat to go   
    kernel.weights <- neighborhood_kernel_raster[raster::cellFromXY(neighborhood_kernel_raster,toAdd_df[,c(1,2)])] 
    # allocate those individuals 
    #toAdd_df$indiv <- 0
    #toAdd <- min(sum(freeSpace),leave_kernel)
    toAdd_df$indiv <- rmultinom(1,leave_kernel,prob=kernel.weights)[,1]    # disperser individuals to add to each cell (favor cells that have more free space...)
  }
  dispndx <- which(toAdd_df$indiv>0)
  updatedStack <- UpdateAbund(focalxy=xy,stack=newStack,df=toAdd_df[dispndx,],DensRaster,FreqList)  # add these disperser individuals and update the accumulating raster stack
  #leave_kernel <- leave_kernel-toAdd   # remove these individuals from the "leave" pool
  #ndx <- raster::cellFromXY(newPop,toAdd_df[,c("x","y")])
  #newPop[ndx] <- newPop[ndx]+toAdd_df$indiv
  # plot(newPop)
  return(updatedStack)
}



####################
# FUNCTION FOR PERFORMING CONSOLIDATION OF LOW-DENSITY POPULATIONS
####################

###################
# DISPERSAL MODE #3: CONSOLIDATE AROUND NEW FOCI (establish new colonies)
###################

#names(neighborhoodRasters)
ConsolidatePops <- function(newStack,DensRaster,FreqList,leave,
                            #newFociRaster=newFociRaster,
                            plagueStatus,
                            xy,
                            newFociRaster,
                            neighborhood_raster,
                            neighborhood_mask_raster,
                            UserParams){
  
  
  DispRates <- GetDispersalRates(plagueStatus,UserParams)
  # plot(newFociRaster)
  
  tempFociRaster <- raster::crop(newFociRaster,raster::extent(neighborhood_mask_raster)) 
  
  # plot(tempFociRaster)
  
  localfoci <- as.data.frame(raster::xyFromCell(tempFociRaster,which(tempFociRaster@data@values==1)))   # data frame of local foci
  nlocalfoci <- nrow(localfoci)              # number of existing new foci in the area
  
  if(nlocalfoci>UserParams$Dispersal$NFOCI) localfoci <- localfoci[sample(c(1:nlocalfoci),UserParams$Dispersal$NFOCI),]   # if multiple local foci, reduce the number of foci to NFOCI
  
  newfocineeded <- max(0,UserParams$Dispersal$NFOCI-nlocalfoci)  # new foci needed
  
  if(newfocineeded>0){
    
    candidates <- neighborhood_raster@data@values   # candidates as new foci
    candidates[candidates>UserParams$Popbio$MINDENS_HA] <- 0   # discount cells that are already fairly full  [note, if no cells are below minimum abundace, then another cell will be selected!]
    candidates[candidates==0] <- 1
    
    for(i in 1:newfocineeded){
      localfoci <- rbind(localfoci,as.data.frame(raster::xyFromCell(tempFociRaster,which.max(candidates))))
      candidates[which.max(candidates)] <- 0
    }
  }
  
  allocation <- rmultinom(1,leave,c(1:UserParams$Dispersal$NFOCI))[,1]
  localfoci$indiv <- allocation 
  
  updatedStack <- UpdateAbund(focalxy=xy,stack=newStack,df=localfoci,DensRaster,FreqList)     # move new individuals to local foci
  
  ## update foci
  ndx <- raster::cellFromXY(newFociRaster,localfoci[,c(1,2)])
  newFociRaster[ndx] <- 1  
  
  #assign(x="newFociRaster",value=newFociRaster,envir=env)
  # plot(newFociRaster)
  retlist <- list()
  retlist$updatedStack <- updatedStack
  retlist$newFociRaster <- newFociRaster
  return(retlist)
}

generateStayFunction <- function(UserParams){
  force(UserParams)
  fun <-  function(x,y){
    disprate <- ifelse(y==1,UserParams$Dispersal$PLAGUE_DISPERSAL,UserParams$Dispersal$BASELINE_DISPERSAL)   # NOTE: this creates "donut holes" in colonies where plague has been...
    
    value <- ifelse(x>UserParams$Popbio$MAXDENS_HA,round(UserParams$Popbio$MAXDENS_HA*(1-disprate)),round(x*(1-disprate)))
    value2 <- ifelse(x<UserParams$Popbio$MINDENS_HA,0,value)
    #if(any(!is.na(c))) value[!is.na(c)] <-   #rpois(length(which(!is.na(c))),c[!is.na(c)])
    return(value2)
  }
  return(fun)
}
  


# stayFunction <- function(x,y,UserParams){
#   disprate <- ifelse(y==1,UserParams$Dispersal$PLAGUE_DISPERSAL,UserParams$Dispersal$BASELINE_DISPERSAL)   # NOTE: this creates "donut holes" in colonies where plague has been...
#   
#   value <- ifelse(x>UserParams$Popbio$MAXDENS_HA,round(UserParams$Popbio$MAXDENS_HA*(1-disprate)),round(x*(1-disprate)))
#   value2 <- ifelse(x<UserParams$Popbio$MINDENS_HA,0,value)
#   #if(any(!is.na(c))) value[!is.na(c)] <-   #rpois(length(which(!is.na(c))),c[!is.na(c)])
#   return(value2)
# }

####################
# MAIN DISPERSAL FUNCTION
####################

#t=which(plagueyear)[1]
doDispersal <- function(UserParams,DensRaster,PlagueRaster,newFociRaster,FreqList,BaseLandscape,DispList){

  # ## build up the results for dispersal
  # newStack <- stack(list(
  #   DensRaster = raster::reclassify(DensRaster,rcl=c(-Inf,Inf,0)),   # blank raster for filling in the dispersal
  #   FreqList = raster::reclassify(FreqList,rcl=c(-Inf,Inf,0))   # blank raster for filling in the plague survival
  # ))    # raster stack that is updated during the dispersal process...
  # 
  ## determine the cells from which to draw dispersers (cells with positive abundance...)

  newDensRaster <- raster::reclassify(DensRaster,rcl=c(-Inf,Inf,0))   # blank raster for filling in the dispersal
  newFreqList <- raster::reclassify(FreqList,rcl=c(-Inf,Inf,0))
  names(newFreqList) <- names(FreqList)
  
  newStack <- list(DensRaster=newDensRaster,
                   FreqList=FreqList)
  
  newStack <- raster::stack(newStack)
  
  # before doing dispersal, determine who stays in place...
  
  sF <- generateStayFunction(UserParams)   # use closure to force arguments to match
  stayRaster <- raster::overlay(DensRaster,PlagueRaster,fun=sF)
    # plot(PlagueRaster)
    # plot(DensRaster)
    # plot(stayRaster)
    #stayRaster[stayRaster<50]
  
  leaveRaster <- DensRaster-stayRaster   # everyone that didn't stay has to leave!
  # plot(leaveRaster)
  
  newStack[['DensRaster']] <- stayRaster   # update the new DensRaster with the known stayers. This deals with all "stay" individuals so that only leavers must be considered
  # update the current rasters with the known stayers
      #rasterNames <- rasterNames #names(newStack)
      #n = rasterNames[2]
  geneNames <- names(FreqList)
  
  n="gene1"              # add the genetics of staying individuals
  for(n in geneNames){
    newStack[[n]][stayRaster>=0] <- newStack[[n]][stayRaster>=0]*stayRaster[stayRaster>=0]     # stayers contribute to post-dispersal allele frequencies
      # plot(eval(parse(text=n)))
      # plot(eval(parse(text=n))*stayRaster)
      # plot(newStack[[n]])
      # extent <- drawExtent()
      # plot(raster::crop(stayRaster,extent))   # okay this is actually working
  }
 
  focalCells <- which(leaveRaster@data@values>0)   # identify cells with dispersers leaving
  xy_df <- data.frame(
    raster::xyFromCell(BaseLandscape$patchRaster,focalCells)
  )
  # reshuffle the order (make foci of expansion a bit more random...)
  xy_df <- xy_df[sample(c(1:nrow(xy_df))),]
  
  # head(xy_df)
  

  #plot(newStack[[rasterNames[2]]])
  counter <- 1  
  focalcell <- focalCells[counter]  # cell from which to draw potential dispersers 
  for(focalcell in focalCells){    # loop through focal cells [actually should only loop through cells that have positive abundance... ]
    xy <- xy_df[counter,] 
    
    plagueStatus <- ifelse(PlagueRaster[focalcell]==1,"plague","noPlague")   # determine if plague
    DispRates <- GetDispersalRates(plagueStatus,UserParams)  # get the current dispersal rates
    
    neighborhoodRasters <- makeNeighborhoodRasters(newStack,plagueStatus, xy,UserParams,DispList)  # characterize the possible sites to move to in the neighborhood
    
    thisAbund <- as.numeric(DensRaster[focalcell])  # abundance in the focal cell
    
    overcrowded <- thisAbund>UserParams$Popbio$MAXABUND # if overcrowded, then dispersal will unfold in a certain way: e.g., colony will expand outward
    undercrowded <- thisAbund<UserParams$Popbio$MINABUND  # if undercrowded, then dispersal can unfold in a different way: e.g., colony will restructure and group together at new focal areas
    
    stay <- stayRaster[focalcell]
    leave <- leaveRaster[focalcell]
    
    if((overcrowded)){   # if focal cell is overcrowded (and not currently experiencing plague)
      #stay <- stayRaster[focalcell] #MAXABUND*(1-disprate)  # number that should stay
      #newStack <- UpdateAbund(focalxy=xy,stack=newStack,df=as.data.frame(cbind(xy,stay)))   # keep the "stay" individuals in place
      #leave <- leaveRaster[focalcell] #max(0,thisAbund-MAXABUND) + MAXABUND*disprate 
      #newPop[focalcell] <- newPop[focalcell] + stay   # make sure that all staying individuals stay put... 
      
      if(plagueStatus=="plague") leave_snuggle <- 0  
      if(plagueStatus=="noPlague") leave_snuggle <- floor(leave*UserParams$Dispersal$SNUGGLE)   # these individuals will try to find a place to settle next door
      leave_kernel <- leave-leave_snuggle     # these individuals will obey the dispersal kernel
      # expand the colony
      newStack <- SpreadOut(newStack,DensRaster,FreqList,
                            leave_snuggle,
                            plagueStatus,xy,neighborhood_raster=neighborhoodRasters[['abundRaster']],
                            neighborhood_mask_raster=neighborhoodRasters[['maskRaster']],
                            freeSpace_raster=neighborhoodRasters[['freeSpaceRaster']],
                            UserParams,
                            DispList)
      
      # do long-distance dispersal (according to the dispersal kernel)
      newStack <- LongDistanceDispersal(newStack=newStack,DensRaster,FreqList,leave_kernel=leave_kernel,
                                        plagueStatus=plagueStatus,
                                        xy=xy,
                                        neighborhood_raster=neighborhoodRasters[['abundRaster']],
                                        freeSpace_raster=neighborhoodRasters[['freeSpaceRaster']],
                                        neighborhood_kernel_raster=neighborhoodRasters[['kernelRaster']],
                                        UserParams)
      
      
    }else if((!overcrowded)&(!undercrowded)){   # if focal cell is not overcrowded or undercrowded
      #stay <- stayRaster[focalcell] #round(thisAbund * (1-BASELINE_DISPERSAL))   # number of individuals staying in the focal cell
      # newStack <- UpdateAbund(focalxy=xy,stack=newStack,df=as.data.frame(cbind(xy,stay)))   # keep the "stay" individuals in place (and update the evolving factors)
      #leave <- leaveRaster[focalcell] thisAbund-stay  # how many are leaving?
      
      # just do long-distance dispersal (according to the dispersal kernel)
      newStack <- LongDistanceDispersal(newStack=newStack,DensRaster,FreqList,leave_kernel=leave,
                                        plagueStatus=plagueStatus,
                                        xy=xy,
                                        neighborhood_raster=neighborhoodRasters[['abundRaster']],
                                        freeSpace_raster=neighborhoodRasters[['freeSpaceRaster']],
                                        neighborhood_kernel_raster=neighborhoodRasters[['kernelRaster']],
                                        UserParams)
      
    }else if((undercrowded)){    # if focal cell is undercrowded
      
      #stay <- stayRaster[focalcell] # in this case, all animals leave  (revisit this!)
      #leave <- stayRaster[focalcell]   #thisAbund
      
      # consolidate into new focal populations
      temp <- ConsolidatePops(newStack,DensRaster,FreqList,leave=leave,
                                  plagueStatus,
                                  xy=xy,
                                  newFociRaster,
                                  neighborhood_raster=neighborhoodRasters[['abundRaster']],
                                  neighborhood_mask_raster=neighborhoodRasters[['maskRaster']],
                                  UserParams)
      newStack <- temp$updatedStack
      newFociRaster <- temp$newFociRaster
      
      
    }
    #if(counter%%100==0) cat(sprintf("%s...",counter))
    counter<-counter+1
  }  # end loop through focal cells
  
#   plot(newStack[["DensRaster"]])
#   plot(newFociRaster)
#   plot(DensRaster)
#   plot(FreqList[[1]])
#   plot(PlagueRaster)
#   plot(newStack[[rasterNames[2]]])     # plague survival
  #   plot(newStack[[rasterNames[3]]])
  
  ##########
  # Expand focal cells
  ##########
  
  # plot(newFociRaster)
  
  foci_ndx <- which((newFociRaster@data@values==1)&(newStack[["DensRaster"]]@data@values>0))
  totfoci <- length(foci_ndx)
  # focalcell <- foci_ndx[6]
  if(totfoci>0){
    names <- paste("gene",1:UserParams$Genetics$NGENES,sep="") 
    focalcell = foci_ndx[1]
    for(focalcell in foci_ndx){  #loop through foci and spread them out!
      
      xy <- as.data.frame(raster::xyFromCell(DensRaster,focalcell))
      
      plagueStatus <- ifelse(PlagueRaster[focalcell]==1,"plague","noPlague")   # determine if plague
      DispRates <- GetDispersalRates(plagueStatus,UserParams)  # get the current dispersal rates
      
      #browser()
      neighborhoodRasters <- makeNeighborhoodRasters(newStack,plagueStatus,xy,UserParams,DispList)  # characterize the possible sites to move to in the neighborhood
      
      thisAbund <- as.numeric(newStack[['DensRaster']][focalcell])  # abundance in the focal cell
      
      # compute the mean values for evolving layers for focal cells so that these values can be spread out
      i = names[2]
      for(i in names){
        if(i!="DensRaster"){
          temp <- newStack[[i]][focalcell]/thisAbund    # get the population mean for key evolving variables
          eval(parse(text=sprintf("FreqList[[\"%s\"]][focalcell]<-temp",i)))   # does this need a double arrow when in function? how to assign this to global environment
          #assign(x=i,value=eval(parse(text=i)), envir = env)   # assign the variable to the global environment
        }
      }
      
      #assign(x="FreqList",value=eval(parse(text="FreqList")), envir = env)   # assign the variable to the global environment
      
      overcrowded <- thisAbund>UserParams$Popbio$MAXABUND # if overcrowded, then dispersal will unfold in a certain way: e.g., colony will expand outward
      
      if((overcrowded)){   # if focal cell is overcrowded 
        stay <- UserParams$Popbio$MAXDENS_HA-thisAbund  # number that should stay (here, negative so that individuals are removed)
        newStack <- UpdateAbund(focalxy=xy,stack=newStack,df=as.data.frame(cbind(xy,stay)),DensRaster,FreqList)   # keep the "stay" individuals in place
        leave <- max(0,thisAbund-UserParams$Popbio$MAXABUND)  
        
        # expand the colony
        newStack <- SpreadOut(newStack,DensRaster,FreqList,
                              leave_snuggle=leave,
                              plagueStatus=plagueStatus,xy=xy,neighborhood_raster=neighborhoodRasters[['abundRaster']],
                              neighborhood_mask_raster=neighborhoodRasters[['maskRaster']],
                              freeSpace_raster=neighborhoodRasters[['freeSpaceRaster']],UserParams,DispList)
      }
    }
  }
  
  
  newFreqList <- UpdateStack(stack=newStack)   # make the new population, get the averages for the evolving layers and save them to env
 
  newPop <- newStack[['DensRaster']]
  
  ret <- list()
  ret$newPop <- newPop
  ret$newFociRaster <- newFociRaster
  ret$newFreqList <- newFreqList
  
  #  DensRaster <- newStack[['DensRaster']]
  return(ret)
}

##################
# FECUNDITY FUNCTIONS
##################

# c=c(1,2,3,4,NA,5)
demographicStoch <- function(c){
  value=c
  if(any(!is.na(c))) value[!is.na(c)] <- rpois(length(which(!is.na(c))),c[!is.na(c)])    # why does this throw errors sometimes?
  return(value)
}

doReproduce <- function(UserParams,DensRaster,PlagueRaster,deviate){
  thisPop <- DensRaster
  #thisFec <- rnorm(1,BASELINE_MEANFEC,CV_FECUNDITY*BASELINE_MEANFEC)
  thisFec <- UserParams$Popbio$BASELINE_MEANFEC + (UserParams$Popbio$CV_FECUNDITY*UserParams$Popbio$BASELINE_MEANFEC)*deviate
  thisFec <- max(0.1,thisFec)  # make sure fecundity is not zero
  #thisPop[(thisPop<MINABUND)] <- 0  # populations below the allee threshold cannot breed
  thisPop[(PlagueRaster==0)&(thisPop>UserParams$Popbio$MINABUND)] <- thisPop[(PlagueRaster==0)&(thisPop>UserParams$Popbio$MINABUND)] + DensRaster[(PlagueRaster==0)&(thisPop>UserParams$Popbio$MINABUND)]*thisFec
  thisPop[(PlagueRaster==1)&(thisPop>UserParams$Popbio$MINABUND)] <- thisPop[(PlagueRaster==1)&(thisPop>UserParams$Popbio$MINABUND)] + (DensRaster[(PlagueRaster==1)&(thisPop>UserParams$Popbio$MINABUND)]*thisFec)/2   # reduced fecundity under plague...
  thisPop <- raster::calc(thisPop,fun=demographicStoch)
  # plot(thisPop)
  return(thisPop)   
}

##################
# SURVIVAL FUNCTIONS
##################

#a=c(NA,NA,0.5,0.7)
#a=c(NA,NA,NA)
# a = c(NA,NA,NA,NA,NA)
getYearVariate <- function(a,deviate,cv){            # note: this function could work for fecundity too...
  sd <- a*cv
  value <- rep(NA,times=length(a))
  if(any(!is.na(sd))){ 
    value[!is.na(sd)] <- a[!is.na(sd)] + deviate*sd[!is.na(sd)]    # for now, spatially correlated... 
    value[value<SURVMIN_NOPLAGUE] <- SURVMIN_NOPLAGUE
    value[value>SURVMAX_NOPLAGUE] <- SURVMAX_NOPLAGUE
  }
  return(value)
}


# getSurvival <- function(resistanceStatus="susceptible",plagueStatus="plague"){
#   survival=0
#   if((plagueStatus=="noPlague")&(resistanceStatus=="susceptible")) survival = BASELINE_MEANSURV
#   if((plagueStatus=="noPlague")&(resistanceStatus=="resistant")) survival = BASELINE_MEANSURV - FITNESS_COST*(BASELINE_MEANSURV*BASELINE_PLAGUESURV_RESIST-SURVMIN_PLAGUE)
#   if((plagueStatus=="plague")&(resistanceStatus=="susceptible")) survival = BASELINE_PLAGUESURV
#   if((plagueStatus=="plague")&(resistanceStatus=="resistant")) survival = BASELINE_MEANSURV*BASELINE_PLAGUESURV_RESIST
#   return(survival)
# }


    # new "getsurvival" function returns a map
getSurvivalFunc <- function(resistanceStatus="susceptible",plagueStatus="plague",FCRaster,UserParams,BaseLandscape){
    #survival=raster::reclassify(BaseLandscape$patchRaster,rcl=c(-Inf,Inf,0))
  if((plagueStatus=="noPlague")&(resistanceStatus=="susceptible")) survival = raster::reclassify(BaseLandscape$patchRaster,rcl=c(-Inf,Inf,UserParams$Popbio$BASELINE_MEANSURV))
  if((plagueStatus=="noPlague")&(resistanceStatus=="resistant")) survival = UserParams$Popbio$BASELINE_MEANSURV - FCRaster*(UserParams$Popbio$BASELINE_MEANSURV*UserParams$Popbio$BASELINE_PLAGUESURV_RESIST-UserParams$Popbio$SURVMIN_PLAGUE)
  if((plagueStatus=="plague")&(resistanceStatus=="susceptible")) survival = survival=raster::reclassify(BaseLandscape$patchRaster,rcl=c(-Inf,Inf,UserParams$Popbio$BASELINE_PLAGUESURV))
  if((plagueStatus=="plague")&(resistanceStatus=="resistant")) survival = survival=raster::reclassify(BaseLandscape$patchRaster,rcl=c(-Inf,Inf,UserParams$Popbio$BASELINE_MEANSURV*UserParams$Popbio$BASELINE_PLAGUESURV_RESIST))
  return(survival)  
}

getMeanSurvival <- function(FCRaster,UserParams,BaseLandscape){
  meansurv <- list()
  for(i in c("resistant","susceptible")){
    meansurv[[i]] <-list()
    for(j in c("plague","noPlague")){
      meansurv[[i]][[j]] <- getSurvivalFunc(i,j,FCRaster,UserParams,BaseLandscape)
    }
  }
  return(meansurv)
}

getSurvival <- function(deviate,cv,UserParams,FCRaster,BaseLandscape){
  surv <- getMeanSurvival(FCRaster,UserParams,BaseLandscape)
  
  for(i in c("resistant","susceptible")){
    for(j in c("plague","noPlague")){
      surv[[i]][[j]] <- surv[[i]][[j]] + deviate*(surv[[i]][[j]]*cv)
      if(j=="plague"){  
        surv[[i]][[j]][surv[[i]][[j]]>UserParams$Popbio$SURVMAX_PLAGUE] <- UserParams$Popbio$SURVMAX_PLAGUE
        surv[[i]][[j]][surv[[i]][[j]]<UserParams$Popbio$SURVMIN_PLAGUE] <- UserParams$Popbio$SURVMIN_PLAGUE
      }
      if(j=="noPlague"){
        surv[[i]][[j]][surv[[i]][[j]]>UserParams$Popbio$SURVMAX_NOPLAGUE] <- UserParams$Popbio$SURVMAX_NOPLAGUE
        surv[[i]][[j]][surv[[i]][[j]]<UserParams$Popbio$SURVMIN_NOPLAGUE] <- UserParams$Popbio$SURVMIN_NOPLAGUE
      }
    }
  }
  return(surv)
}

     # note: maybe we should just break out poparray by resistance during the survival function... Otherwise just frequencies... 
     #     in that case, we need to update frequencies here too. This is where "enrichment" happens.
doSurvival <- function(DensRaster,PlagueRaster,FreqList,UserParams,BaseLandscape,deviate,cv){   # PopArray=PopArray
  #thisPop <- raster::getValues(PopArray)
  
  temp <- GetStructuredPop(DensRaster,FreqList,UserParams)    # break out population into resistant and non-resistant (and account for resistance factors)
  thisPop <- temp$Pop
  
  allelefreq <- temp$allelefreq
  
    # structuredFreq <- GetStructuredFreqList(DensRaster,FreqList)  # break out resistance factors into resistant and non-resistant 
  
  FCRaster <- FitnessCost(FreqList,UserParams)      # compute fitness costs
  
  surv <- getSurvival(deviate,cv,UserParams,FCRaster,BaseLandscape)    # get survival for all possible combinations of resistance and plague
  
  
  ###########
  #  PERFORM SURVIVAL
  ###########
  
  status="resistant"
  for(status in c("resistant","susceptible")){
    thisPop[[status]][PlagueRaster==1] <- thisPop[[status]][PlagueRaster==1]*surv[[status]][["plague"]][PlagueRaster==1]
    thisPop[[status]][PlagueRaster==0] <- thisPop[[status]][PlagueRaster==0]*surv[[status]][["noPlague"]][PlagueRaster==0]
  }
  
  ##########
  #  DEMOGRAPHIC STOCHASTICITY
  ##########
  
  for(status in c("resistant","susceptible")){
    thisPop[[status]] <- raster::calc(thisPop[[status]],fun=demographicStoch)
  }
  
  
  ##########
  #  REVERT TO UNSTRUCTURED POPULATION
  ##########
  
  temp <- GetUnstructuredPop(thisPop,UserParams,FreqList,allelefreq)    # break out population into resistant and non-resistant (and account for genes/resistance factors)
  FreqList <- temp$FreqList
  NewDensRaster <- temp$Dens
  
  ret <- list()
  ret$FreqList <- FreqList
  ret$NewDensRaster <- NewDensRaster
  return(ret)
}

##################
# DD SURVIVAL FUNCTION
##################

doDDSurvival <- function(DensRaster,UserParams){
  
  thisPop <- DensRaster
  #for(status in c("resistant","susceptible")){
    thisPop[DensRaster>(UserParams$Popbio$MAXDENS_HA*1.15)] <- UserParams$Popbio$MAXDENS_HA*1.15  # kill off all individuals in populations above the threshold
  #}
  return(thisPop)
  
}


##################
# ALLEE FUNCTION
##################

doAllee <- function(){
  
  thisPop <- DensRaster
  thisPop[DensRaster<MINABUND] <- 0  # kill off all individuals in populations below the threshold
  return(thisPop)
  
}



#####################
# FUNCTIONS FOR DETERMINING RESISTANCE STATUS FROM GENOTYPE
#####################


# this function takes gene frequencies and converts to "factor" frequencies...
Gene2Factor <- function(UserParams,FreqList){  # FreqList
  newlist <- list()
  gene=1
  for(gene in 1:UserParams$Genetics$NGENES){
    name = sprintf("gene%s",gene)
    name2 = sprintf("factor%s",gene)
    newlist[[name2]] =  FreqList[[name]]^2 * UserParams$Genetics$DOMINANCE[gene,1] +
      (2*FreqList[[name]]*(1-FreqList[[name]])) * UserParams$Genetics$DOMINANCE[gene,2]
  }
  newlist <- raster::stack(newlist)
  return(newlist)
}

# this function takes gene frequencies and abundances and converts to resistant freqs and susceptible freqs...
StrFreqFunc <- function(DensRaster,ResistRaster,FreqList,UserParams){
  reslist <- list()
  suslist <- list()
  #  temp <- DensRaster*F
  gene=1
  PerHet <- ((2/FreqList)-2)/(((2/FreqList)-2)+1)    # percent heterozygotes for dominant factors in resistant pool
  for(gene in 1:UserParams$Genetics$NGENES){
    name = sprintf("gene%s",gene)
    suslist[[name]] <- 2*DensRaster*FreqList[[name]]  # total resistance alleles in the population
    reslist[[name]] <- raster::reclassify(FreqList[[name]],rcl=c(-Inf,Inf,0))
    
    # number of resistance alleles in the resistant pool for each gene
    if(sum(UserParams$Genetics$DOMINANCE[gene,])==2){   # if fully dominant
      reslist[[name]] <- reslist[[name]] + (2*ResistRaster) * (1-PerHet[[name]])  +
        ResistRaster*PerHet[[name]]
      suslist[[name]] <- suslist[[name]] - reslist[[name]]  # number of resistance alleles remaining in the susceptible population    
    }else{
      reslist[[name]] <- reslist[[name]] + (2*ResistRaster) * UserParams$Genetics$DOMINANCE[gene,1] +
        (ResistRaster) * UserParams$Genetics$DOMINANCE[gene,2]
      suslist[[name]] <- suslist[[name]] - reslist[[name]]  # number of resistance alleles remaining in the susceptible population
    }  
    reslist[[name]] <- reslist[[name]]/(2*ResistRaster)  # allele frequency in res pop
    suslist[[name]] <- suslist[[name]]/(2*(DensRaster-ResistRaster))  # freq in sus pop
  }
  reslist <- raster::stack(reslist)
  suslist <- raster::stack(suslist)
  retlist <- list()
  retlist$reslist <- reslist
  retlist$suslist <- suslist
  return(retlist)
  # assign(x="reslist",value=reslist, envir = env)   # assign the variable to the global environment
  # assign(x="suslist",value=suslist, envir = env)   # assign the variable to the global environment
}

resistfunc <- function(ngenes=UserParams$Genetics$NGENES){      # this function assumes that all factors are needed for resistance                    
  nargs <- ngenes
  arguments <- paste("X",c(1:nargs),sep="")
  arguments2 <- paste(arguments, collapse=",")
  arguments3 <- paste(arguments, collapse="*")
  expression <- sprintf("function(%s) %s ",arguments2,arguments3)
  eval(parse(text=expression))
}


#UserParams$Genetics$RESISTANCE_SCENARIOS[[1]]


IsResistant <- function(DensRaster,FreqList,fungen,UserParams){
  FactorList <- Gene2Factor(UserParams,FreqList)  # FreqList
  temp <- raster::overlay(FactorList,fun=fungen(UserParams$Genetics$NGENES)) #round(InitDensRaster*InitFreq[["gene1"]])   # freq of resist for each grid cell
  ResistRaster <- raster::overlay(DensRaster,temp,fun=function(x,y) x*y)      # numbers of resistant individuals in each grid cell
  temp <- StrFreqFunc(DensRaster,ResistRaster,FreqList,UserParams) # returns "suslist" and "reslist" to global env
  
  
  suslist <- temp$suslist

  reslist <- temp$reslist
  ret <- list()
  ret$allelefreq <- list()
  ret$allelefreq$sus <- suslist
  ret$allelefreq$res <- reslist
  ret$ResistRaster <- ResistRaster
  return(ret)
}

# use a closure    
FCfunc <- function(ngenes,UserParams){   # assume that fitness costs are simply additive
  nargs <- ngenes
  force(UserParams)  #??
  arguments <- paste("X",c(1:nargs),sep="")
  arguments2 <- paste(arguments, collapse=",")
  arguments3 <- paste(paste("UserParams$Genetics$FITNESS_COST[",c(1:ngenes),"]*",arguments,sep=""),collapse="+")
  
  #FITNESS_COST[1]*X1 + FITNESS_COST[2]*X2 
  expression <- sprintf("function(%s) %s ",arguments2,arguments3)
  eval(parse(text=expression))
}

FitnessCost <- function(FreqList,UserParams){
  newfunc <- FCfunc(UserParams$Genetics$NGENES,UserParams)
  FCraster <- raster::overlay(FreqList,fun=newfunc)  # degree of fitness cost
  return(FCraster)
}


# function for computing the allele (factor) frequecies for the structured population (resistant and susceptible...)


# use information on frequencies of resistance factors to struture population into resistance categories
GetStructuredPop <- function(DensRaster,FreqList,UserParams){
  Pop <- list()
  temp <- IsResistant(DensRaster,FreqList,resistfunc,UserParams)      # structure by susceptible and resistant
  Pop[["resistant"]] <-  temp$ResistRaster
  Pop[["susceptible"]] <- DensRaster - Pop[["resistant"]]
  allelefreq <- temp$allelefreq
  Pop <- raster::stack(Pop)
  ret <- list()
  ret$Pop <- Pop
  ret$allelefreq <- allelefreq
  return(ret)
}


# use information on frequencies of resistance factors to struture population into resistance categories
GetUnstructuredPop <- function(PopArray,UserParams,FreqList,allelefreq){
  Dens <- raster::overlay(PopArray,fun=sum)    # get total population size
  i=1
  for(i in 1:UserParams$Genetics$NGENES){
    name <- sprintf("gene%s",i)
    ndx1 <- Dens@data@values>0
    ndx2 <- Dens@data@values==0
    FreqList[[name]][ndx1] <- (allelefreq$sus[[name]][ndx1]*PopArray[["susceptible"]][ndx1] + allelefreq$res[[name]][ndx1]*PopArray[["resistant"]][ndx1])/Dens[ndx1]
    FreqList[[name]][ndx2] <- 0
  }
  ret <- list()
  ret$FreqList <- FreqList
  #assign("FreqList",FreqList,envir=env)
  ret$Dens <- Dens
  return(ret)
}


#######################
# INITIALIZE ALLELE FREQUENCIES / RESISTANCE FACTORS [keep for now- will be multiple genes in the model somehow]
#######################
# NOTE: Some regions are more likely to evolve faster because they have greater percentages of those genes that can confer resistance. 

GetInitFreqs <- function(UserParams,BaseLandscape){

  InitFreqList <- list()
  
  #UserParams$Genetics$INITFREQ
  
  i=1
  for(i in 1:UserParams$Genetics$NGENES){
    name <- sprintf("gene%s",i)
    temp <- rnorm(BaseLandscape$nPatches,UserParams$Genetics$INITFREQ[1],UserParams$Genetics$INITFREQ_SD[1])
    temp <- ifelse(temp<0,0.01,temp)
    InitFreqList[[name]] <- raster::reclassify(BaseLandscape$patchIDRaster,rcl=cbind(c(1:BaseLandscape$nPatches),temp))
  }
  
  #plot(InitFreqList[["gene1"]])
  InitFreqList <- raster::stack(InitFreqList)
  
  #plot(InitFreqList)
  return(InitFreqList)
}



