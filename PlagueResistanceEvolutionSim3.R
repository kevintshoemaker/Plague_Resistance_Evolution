########################################
##   DISEASE RESISTANCE EVOLUTION MODEL
########################################

# time step is one year

# consider scalar demographic model for now, but with potential for matrix model

# 8/29/16: changed to consider a binary resistance trait that is potentially a function of several gene loci... 

############
## CLEAR WORKSPACE
############

rm(list=ls())


############
## LOAD PACKAGES
############

library("raster")
library("secr")
library("igraph")

############
## GLOBAL VARIABLES: DEFINE RANGE OF POSSIBLE SCENARIOS  [not implemented yet]
############

# Bg : "background"  - for looking at intermittent selection...

BgDiseaseSpatialExtents <- c("gridcell","subpatch","patch","regional","landscape")   
BgDiseaseFrequencies <- c(1,2,3,5,10,20,30,40,50)

  # ProbTransmissions <- seq(0,1,by=0.1)   # for now, no transmission


############
## USER-DEFINED VARIABLES
############

          # define the gridded landscape
NROWS <- 50
NCOLS <- 50

CELLAREA <- 1     # in hectares
CELLAREA_M2 <- CELLAREA*10000   # in square meters

CELLWIDTH <- sqrt(CELLAREA_M2)    # in meters
HALFCELLWIDTH <- CELLWIDTH/2

          # define the clustering (degree to which species prefers to establish residence near members of its own kind)
               # 1 is complete tendency to cluster in space. 0 is agnostic to members of its own kind. -1 is tendency to avoid members of its own kind
SNUGGLE <- 0.75

          # define the percent of the landscape that is suitable
PER_SUITABLE <- 0.4

          # define the maximum per-cell number of individuals
MAXDENS <- 100
MAXABUND <- MAXDENS * CELLAREA

         # define the minimum per-cell number of individuals (provide a simple hard Allee effect)
MINDENS <- 15
MINABUND <- MINDENS*CELLAREA

         # maximum survival under plague (limit to resistance)
#   RESISTANCE_LIMIT <- 0.75   # deprecated  

         # baseline survival in a non_plague year (baseline survival for a naive population under no plague)
BASELINE_MEANSURV <- 0.6

         # Variation in survival among years, expressed as CV
CV_SURVIVAL <- 0.2 

         # fecundity in a non-plague year (baseline, number of offspring per adult, not sex structured)
BASELINE_MEANFEC <- 3.2

         # temporal variation in fecundity, expressed as a CV
CV_FECUNDITY <- 0.5 

         # minimum survival under plague (completely naive population)
BASELINE_PLAGUESURV <- 0.05

         # survival under plague for resistant individuals
BASELINE_PLAGUESURV_RESIST <- 0.5

         # mean change in the ability to withstand plague among the survivor population. (accounts for limited heritability)
#SURVCHANGE_NEXTPLAGUE <- 0.05  # deprecated

         # variation in the ability to survival plague among survivors, expressed as a SD [variation in standing genetic propensity to evolve resistance]
#SD_SURVCHANGE_NEXTPLAGUE <- 0.05    # deprecated

         # plague survival rate as a function of the resistance factors?? [TODO] 

         # change in the ability to survive in a non-plague year, as a percentage of the fitness benefit in a plague year (fitness costs to resistance)
#FITNESS_COST <- 0.1    # now means the degree to which survival is reduced for resistant individuals in nonplague years

NGENES <- 2

NWAYS_OF_GETTING <- 1   # ways of getting resistance
RESISTANCE_SCENARIOS <- list()
RESISTANCE_SCENARIOS[[1]] <- c("factor1","factor2")   # for now, needs both factors! 
names(RESISTANCE_SCENARIOS[[1]]) <- c("AND","AND")   # names follow boolean conventions. In this case, both factors are required for resistance

FITNESS_COST <- numeric(NGENES)

FITNESS_COST[1] <- 0.1     # fitness cost of the first gene
FITNESS_COST[1] <- 0.05     # fitness cost of the second gene

INITFREQ <- numeric(NGENES)
INITFREQ[1] <- 0.09
INITFREQ[2] <- 0.1

INITFREQ_SD <- 0.03     # degree of variation in initial frequency of resistance.

         # maximum annual dispersal distance (m)  
MAXDISPERSAL <- 500

         # number of colony foci to establish (within the dispersal range) after colonies get below the low-density threshold
NFOCI <- 1

         # rate of transmission (per-disperser probability of initiating an outbreak in recipient population) (under optimal plague conditions?)
PROB_TRANSMISSION <- 0  # [not implemented yet...] [maybe don't need...] [but critical if looking at tradeoffs]

         # dispersal rate (independent of density)
BASELINE_DISPERSAL <- 0.05

         # dispersal distance (m) for plagued-out populations
MAXDISPERSAL_PLAGUE <- 1000   # individuals from plagued-out populations might move farther than normal population

         # dispersal rate for plagued-out populations
PLAGUE_DISPERSAL <- 0.95     # individuals from plagued-out populations might have a higher tendency to move than normal populations- this can affect the spread of plague and the spread of plague resistance genes... 

         # minimum survival for non-plague populations
SURVMIN_NOPLAGUE <- 0.1

         # maximum survival for non-plague populations
SURVMAX_NOPLAGUE <- 0.9

         # minimum survival for plague populations
SURVMIN_PLAGUE <- 0.01

         # maximum survival for plague populations (circumscribing env stochasticity at the population)
SURVMAX_PLAGUE <- 0.75



############
## SIMULATION CONTROLS
############

NYEARS <- 20

############
## SET BASE DIRECTORY
############

KEVIN_LAPTOP <- TRUE # FALSE #   
KEVIN_OFFICEPC <- FALSE # TRUE # 

if(KEVIN_LAPTOP) BASE_DIR <- "C:\\Users\\Kevin\\Dropbox\\PlagueModeling\\ResistanceEvolution"
if(KEVIN_OFFICEPC) BASE_DIR <- "E:\\Dropbox\\PlagueModeling\\ResistanceEvolution"

if(KEVIN_LAPTOP) GIT_DIR <- "C:\\Users\\Kevin\\GIT\\Plague_Resistance_Evolution"
if(KEVIN_OFFICEPC) BASE_DIR <- "E:\\GI\\Plague_Resistance_Evolution"

############
## SET UP WORKSPACE (define global variables)
############

 # RSCRIPT_DIR <- sprintf("%s\\Rscripts",GIT_DIR)
DATA_DIR <- sprintf("%s\\Data",BASE_DIR)

FIGS_DIR <- sprintf("%s\\RawFigs",BASE_DIR)

setwd(DATA_DIR)

#####################
# LOAD FUNCTIONS
#####################

setwd(GIT_DIR)
source("PlagueResistanceEvolution_FUNCTIONS.R")

#####################
# INITIALIZE DISPERSAL   (for both plague and no plague... )   
#####################

InitializeDispersal()

# plot(tempmask)
# ?make.grid
# ?make.mask
# ?randomHabitat
# 
# plot(temppatches)
# attributes(temppatches)
# class(temppatches)
# names(temppatches)
# covariates(temppatches)$habitat

########################
# INITIALIZE LANDSCAPE
########################

InitializeLandscape(solid=F)   # generate patchmaps etc.

plot(patchRaster)
plot(patchIDRaster)

########################
# GET PLAGUE MODEL
########################

PlagueModel <- GetPlagueModel()    # for now, use fake plague model- will be a statistical model

######################################
#########################
######  INITIALIZE!
#########################


########################
# INITIALIZE POPULATION
########################

InitDensRaster <- KRaster    # initialize abundance at carrying capacity
plot(InitDensRaster)


#######################
# INITIALIZE RESISTANCE FACTORS [keep for now- will be multiple genes in the model somehow]
#######################

# NFACTORS = 1
# 
# ResistanceFactors <- list()
# 
# for(i in 1:NFACTORS){
#   temp <- rnorm(nPatches,SURVCHANGE_NEXTPLAGUE,SD_SURVCHANGE_NEXTPLAGUE)
#   temp <- ifelse(temp<0,0,temp)
#   PlagueResistancePotentialRaster <- reclassify(patchIDRaster,rcl=cbind(c(1:nPatches),temp))
#   # plot(PlagueResistancePotentialRaster)
#   ResistanceFactors[[i]] <- PlagueResistancePotentialRaster   # now this represents something like allele frequency
# } 
# 
# ResistanceFactors <- stack(ResistanceFactors)

#######################
# INITIALIZE PLAGUE SURVIVAL [deprecated]
#######################

# NextPlagueSurvRaster <- patchRaster*BASELINE_PLAGUESURV
# plot(NextPlagueSurvRaster)

#######################
# INITIALIZE NORMAL SURVIVAL [deprecated?]
#######################

# NextNormalSurvRaster <- patchRaster*BASELINE_MEANSURV
# plot(NextNormalSurvRaster)

#######################
# STRUCTURED SURVIVAL?  
#######################

#  NOTE: ultimately, each gene will have its own fitness cost?? [like it!]

#######################
# INITIALIZE ALLELE FREQUENCIES
#######################
# NOTE: Some regions are more likely to evolve faster because they have greater percentages of those genes that can confer resistance. 

# temp <- rnorm(nPatches,SURVCHANGE_NEXTPLAGUE,SD_SURVCHANGE_NEXTPLAGUE)
# temp <- ifelse(temp<0,0.01,temp)
# PlagueResistancePotentialRaster <- reclassify(patchIDRaster,rcl=cbind(c(1:nPatches),temp))    # deprecate? change to initial frequency?
# plot(PlagueResistancePotentialRaster)  # now this can be an allele frequency map

InitFreqList <- list()

for(i in 1:NGENES){
  name <- sprintf("gene%s",i)
  temp <- rnorm(nPatches,INITFREQ[1],INITFREQ_SD[1])
  temp <- ifelse(temp<0,0.01,temp)
  InitFreqList[[name]] <- reclassify(patchIDRaster,rcl=cbind(c(1:nPatches),temp))
}

plot(InitFreqList[["gene1"]])

InitFreqList <- stack(InitFreqList)

plot(InitFreqList)

#####################
# FUNCTION FOR DETERMINING RESISTANCE STATUS FROM GENOTYPE
#####################

# ultimately this will be a function of underlying binary genotypes. 
     # e.g., if factor A is present in 20% of the pop and factor B is present in 5% of the pop and 
     # both factors are necessary for resistance, then 1% of the population will be resistant. This function could be complex.
     #  i.e., either factors (A and B) OR (C and D) lead to resistance.

# this can be simple as long as we don't assume linkage or some such thing!

# t <- c(NA,NA,)
# IsResistantFunc <- function(abundraster,freqstack){
#   
#       # for now, assume you need both genes to be resistant
# }

    # use a closure?

# TODO: account for multiple different ways of getting resistance
# TODO: keep track of frequencies of each factor within the resistant and susceptible pools (so they can be enriched appropriately)
resistfunc <- function(ngenes){      # this function assumes that all factors are needed for resistance                    
  nargs <- ngenes
  arguments <- paste("X",c(1:nargs),sep="")
  arguments2 <- paste(arguments, collapse=",")
  arguments3 <- paste(arguments, collapse="*")
  expression <- sprintf("function(%s) %s ",arguments2,arguments3)
  eval(parse(text=expression))
}


IsResistant <- function(DensRaster=InitDensRaster,FreqList=InitFreqList,fungen=resistfunc){
  temp <- overlay(FreqList,fun=fungen(NGENES)) #round(InitDensRaster*InitFreq[["gene1"]])   # freq of resist for each grid cell
  ResistRaster <- overlay(DensRaster,temp,fun=function(x,y) x*y)      # numbers of resistant individuals in each grid cell
  return(ResistRaster)
}

    # use a closure    
FCfunc <- function(ngenes){   # assume that fitness costs are simply additive
  nargs <- ngenes
  arguments <- paste("X",c(1:nargs),sep="")
  arguments2 <- paste(arguments, collapse=",")
  arguments3 <- paste(paste("FITNESS_COST[",c(1:ngenes),"]*",arguments,sep=""),collapse="+")
  
      #FITNESS_COST[1]*X1 + FITNESS_COST[2]*X2 
  expression <- sprintf("function(%s) %s ",arguments2,arguments3)
  eval(parse(text=expression))
}

FitnessCost <- function(FreqList=InitFreqList){
  FCraster <- overlay(FreqList,fun=FCfunc(NGENES))  # degree of fitness cost
  return(FCraster)
}


   # function for computing the allele (factor) frequecies for the structured population (resistant and susceptible...)



#####################
# INITIALIZE POPULATION
#####################


### Code block for pop starting from small loci

InitDensRaster2 <- reclassify(patchIDRaster,rcl=c(-Inf,Inf,0))    # for testing
ndx <- sample(which(!is.na(InitDensRaster2@data@values)),size=3)
InitDensRaster2[ndx] <- 1000   # initialize population in random locations
InitDensRaster <- InitDensRaster2
#PopArray2 <- InitDensRaster   # copy, for dispersal algorithm... 


GetStructuredPop <- function(DensRaster=InitDensRaster,FreqList=InitFreqList){
  Pop <- list()
  Pop[["resistant"]] <- IsResistant(DensRaster,FreqList)      # structure by susceptible and resistant. 
  Pop[["susceptible"]] <- DensRaster - Pop[["resistant"]]
  Pop <- stack(Pop)
  return(Pop)
}

        # use information on frequencies of resistance factors to struture population into resistance categories
GetStructuredPop <- function(DensRaster=InitDensRaster,FreqList=InitFreqList){
  Pop <- list()
  Pop[["resistant"]] <- IsResistant(DensRaster,FreqList)      # structure by susceptible and resistant. 
  Pop[["susceptible"]] <- DensRaster - Pop[["resistant"]]
  Pop <- stack(Pop)
  return(Pop)
}


PopArray <- GetStructuredPop(InitDensRaster)

plot(PopArray)
  

# keep track of resistance frequencies in each 
GetStructuredPop <- function(DensRaster=InitDensRaster,FreqList=InitFreqList){
  Freq <- list()
  Freq[["resistant"]] <- IsResistant(DensRaster,FreqList)      # structure by susceptible and resistant. 
  Pop[["susceptible"]] <- DensRaster - Pop[["resistant"]]
  Pop <- stack(Pop)
  return(Pop)
}

######################
# INITIALIZE PLAGUE PROCESS   [KTS: moving away from this and towards a statistical model]
######################
#  for now, assume that plague hits at the patch level, and is a random process.

# PROB_PLAGUE_YEAR <- 0.5 # probability that a plague event hits (landscape level?)   # for now, plague only hits one patch in a plague year
# 
# plagueyear = as.logical(rbinom(NYEARS,1,PROB_PLAGUE_YEAR))
# 
# plagueNow = floor(runif(NYEARS,1,nPatches+1)) * as.numeric(plagueyear)    ## which patch plagues out?

PlagueRaster_template <- reclassify(patchIDRaster,rcl=c(-Inf,Inf,0))   

PlagueRaster <- doPlague(PlagueRaster=PlagueRaster_template, PopArray=reclassify(patchIDRaster,rcl=c(-Inf,Inf,0)))


plot(PlagueRaster)


####################
# INITIALIZE SURVIVAL
####################

meansurv <- matrix(0,nrow=2,ncol=2)    # survival matrix (mean)
rownames(meansurv) <- c("resistant","susceptible")
colnames(meansurv) <- c("plague","noPlague")
meansurv["resistant","noPlague"] <- getSurvival("resistant","noPlague")
meansurv["resistant","plague"] <-   getSurvival("resistant","plague")
meansurv["susceptible","plague"] <- getSurvival("susceptible","plague")
meansurv["susceptible","noPlague"] <- getSurvival("susceptible","noPlague")
meansurv

####################
# START LOOP THROUGH YEARS
####################


    # names of important raster maps to save to file etc...
rasterNames  <- c("PopArray","NextPlagueSurvRaster","NextNormalSurvRaster","PlagueResistancePotentialRaster")   # Deprecate?
  
# t=which(plagueyear)[1]
t=0

t=t+1
for(t in 1:(NYEARS)){
  deviate <- rnorm(1)   #determine if this is a good year or a bad year (for now, survival and fecundity are perfectly correlated)
  cv=CV_SURVIVAL   # set up for using the getYearVariate function
  
#   ###################
#   # EXTERNAL PLAGUE PRESSURE
#   ###################
#   if(plagueyear[t]){
#     PlagueRaster <- PlagueRaster_template
#     PlagueRaster[patchIDRaster==plagueNow[t]] <- 1
#     plot(PlagueRaster)
#   }else{
#     PlagueRaster <- PlagueRaster_template
#   }
  

  ##################
  # DENSITY INDEPENDENT SURVIVAL (including plague survival)
  ##################
  
  PopArray <- doSurvival(PopArray=PopArray,PlagueRaster=PlagueRaster)
  # plot(PopArray)
#   plot(NextPlagueSurvRaster)
#   plot(NextNormalSurvRaster)
# plot(PlagueResistancePotentialRaster)
  
  ################
  # REPRODUCTION
  ################
  
  PopArray <- doReproduce(PlagueRaster = PlagueRaster)
  # plot(PopArray)
  
  
  ###############
  # DISPERSAL: Move individuals around the landscape (this takes a while!)
  ###############
  PopArray <- doDispersal(t=t,PlagueRaster=PlagueRaster)
  # plot(PopArray)    # good in t=1, not so much in t=2
  
  ###############
  # ALLEE EFFECT: REMOVE POPULATIONS BELOW A MINIMUM ABUNDANCE THRESHOLD
  ###############
    #PopArray <- doAllee()
  #if(MINABUND>0) PopArray <- doAllee()   # don't need this! all low-dens individuals move out anyway...
  # $
  
  ###############
  # CLEAR excess individuals from cells (DD)
  
  PopArray <- doDDSurvival()
  # plot(PopArray) 
  
  ##################
  # PLAGUE PRESSURE
  ##################
  
  PlagueRaster <- doPlague(PlagueRaster=PlagueRaster,PopArray = PopArray)
  # plot(PlagueRaster)
  
  
  
  
  # # plot(PopArray)
  # 
#   plot(PopArray)    # okay
#   plot(reclassify(NextNormalSurvRaster,rcl=c(NA,NA,0)))   # okay in t=1,2
#   plot(NextNormalSurvRaster)
#   plot(reclassify(NextPlagueSurvRaster,rcl=c(NA,NA,0)))    # okay t=1,2
#   
#   plot(newFociRaster)   # okay t=1
#   
#   plot(PlagueRaster)
  
  ################
  # MAKE PLOTS
  
  width = 500
  height= 500
  
      # abundance figure
  setwd(FIGS_DIR)
    file = sprintf("AbundanceFig_year%s.tif",t)
    tiff(file, width=width,height=height)
    plot(patchRaster,col=gray(0.7),legend=F)
     #col = colorRampPalette(c("red","red"))(1)
    col = rgb(0,seq(0,1,length=10),0)
    plot(reclassify(PopArray,rcl=c(-Inf,5,NA)),add=T,legend=T)
    #plot(reclassify(PlagueRaster,rcl=c(-Inf,0.01,NA)),col=rgb(1,0,0),add=T,alpha=0.5,legend=F)
    #plot(reclassify(NextPlagueSurvRaster,rcl=c(-Inf,0.001,NA)),col=heat.colors(10),add=T,legend=T)
  dev.off()

      # evolution figure
  setwd(FIGS_DIR)
    file = sprintf("PlagueSurvFig_year%s.tif",t)
    tiff(file, width=width,height=height)
    plot(patchRaster,col=gray(0.7),legend=F)
    #col = colorRampPalette(c("red","red"))(1)
    col = rgb(0,seq(0,maxValue(NextPlagueSurvRaster),length=10),0)
    plot(reclassify(NextPlagueSurvRaster,rcl=c(-Inf,0.001,NA)),add=T,col=col,legend=T)
    #plot(reclassify(PlagueRaster,rcl=c(-Inf,0.01,NA)),col=rgb(1,0,0),add=T,alpha=0.5,legend=F)
    #plot(reclassify(NextPlagueSurvRaster,rcl=c(-Inf,0.001,NA)),col=heat.colors(10),add=T,legend=T)
  dev.off()  
  
    # plague figure
  setwd(FIGS_DIR)
  file = sprintf("PlagueFig_year%s.tif",t)
  tiff(file, width=width,height=height)
  plot(patchRaster,col=gray(0.7),legend=F)
  #col = colorRampPalette(c("red","red"))(1)
  col = rgb(seq(0,maxValue(PlagueRaster),length=10),0,0)
  plot(reclassify(PlagueRaster,rcl=c(-Inf,0.01,NA)),add=T,col=col,legend=T)
  #plot(reclassify(PlagueRaster,rcl=c(-Inf,0.01,NA)),col=rgb(1,0,0),add=T,alpha=0.5,legend=F)
  #plot(reclassify(NextPlagueSurvRaster,rcl=c(-Inf,0.001,NA)),col=heat.colors(10),add=T,legend=T)
  dev.off()
  
  
  #gray.colors(10)
  
  # setwd(FIGS_DIR)
  # date <- Sys.Date()
  # for(i in rasterNames){
  #   filename = sprintf("%s_time%s_%s.tif",i,t,date)
  #   thisRaster <- eval(parse(text=i))
  #   writeRaster(thisRaster,filename=filename,format='GTiff',overwrite=T)
  # }
  
  # test <- raster(filename)
  # plot(test)
   
}


?rgb



### strategy. As data, there will be a host of genes that are identified as potentially related to resistance (candidate loci). 
# some of these genes may be apparently related to resistance in the field (rangewide screeing) and lab (challenge experiments... )
# some of these genes may be apparently related to resistance in the field on the basis of the SPV trials



### To model this is tricky. We are not talking about just one gene probably. We are talking about a messier situation
###    Some of these resistance alleles may be dominant, others recessive. That is, the heterozygote may in some cases be killed by plague just like susceptibles.
###       in other cases, the heterozygote may survive plague similar to a homozygote "resistant"
###       in other cases, the heterozygote may be resistant, while both homozygotes may be susceptible (overdominance)
###       in other cases, the heterozygote may be killed, while the homozygotes may survive (underdominance)
###       in other cases, there may be multiple alleles at certain key loci at a population level, each of which has different effect sizes and dominance.
###       in other cases, genes may act in concert to confer resistance (epistasis)
###     ***dominance of resistance can have an important role to play in influencing the rate of resistance evolution

###  Needless to say, this is complex! To address this issue, we assume that 
# there is a maximum resistance level (biological limit to plague survival given the current gene pool) (assume no mutation)
# The mixture of dominance and epistasis and linkage means that the degree to which resistance factors are passed from generation to generation are limited... 
# There will be much uncertainty as to epistasis and dominance. 
# no indiv will have all resistance factors- so gene interactions and such will be uncertain. We won't have this information
# also, the results from the lab challenges may differ substantially from the results from the field screening. The genes that appear to evolve in the
#   screening may differ from the genes that appear to affect survival in lab challenges.
# This difference (i.e., between lab and field) could be very telling and an interesting thing to make hypotheses about...
# keep track of abundance for all possible genotypes? BUT this gets very computationally intensive for large numbers of genes.
# compute survival for all possible genotype? we won't have this information.


##  NEW STRATEGY:
# for each cell keep track of the current survival rate in the hypothetical case of plague and the current survival rate in the hypothetical case of no plague.  

 

## new strategy will better accommodate the data we have. We have multi-locus allele frequencies before and after plague, and we have
# the degree to which the resulting gene pool might improve survival. We possibly have survival over several plague outbreaks (see if the
# survival rate increases with each plague outbreak?). We have information on colony collapse for many years, which may help to indicate
# how fast resistance is evolving in this system...


## HYPOTHESES

# Endemic plague can play a strong role in maintaining resistance during periods of no plague.   [dynamics in endemic vs naive?]
# Intermittent selection coupled with fitness costs can explain the plague equilibrium in the Asian case
# Effective evolution of resistance depends on large-scale outbreaks- beyond what is possible if pdogs were the primary vector.
# Is the asian system evidence of fitness costs- is there any other way to make that happen?

## QUESTIONS

# Is it possible that we could look at other factors influencing disease severity, like other diseases
# can we evaluate the breakdown of social units under plague and how that contributes to disease spread and the evolution of resistance?

# can we model alternative hosts? Pool all rodents together?? 
# what about alternative dispersers of infected flea vectors? What does that do?

# how does the degree of resistance in the alternative hosts influence the rate of evolution?

# should we include a flea population also? Three populations: flea, rodent, and focal host? 

# either resistance can increase flea load and host competence or decrease flea load and host competence.. 

# for plague to spread, hosts need to transmit to fleas, which transmit to other 

# focus on novel vector-borne disease?

# do we need an IBM to complement this landscape model? - influence of resistance on transmission rate?
 # - influence of pdog social behavior on transmission rate? Breakdown of sociality?
 # - clipping? grooming?


# we need to keep track of individual genes because that's the data we will have. (how to link this data with this model?)

# should we account for the percentage of dispersers that establish and mate in the recipient population??

# should plague resistance potential decrease over time??  Simulating the loss of evolutionary potential under strong selection?

# effect of temporal/spatial variance in selection pressures / potential for resistance / fitness costs ??

# what is the role of Allee effects?


# can we verify the model of colonization and colony expansion using the rangewide data synthesis?  Probably, right??

##### NOTES

# might include strong allee effect. At low densities, more susceptible to predators etc. At high densities, populations actually do better, until they reach the threshold.. 

# should resistance acquisition be proportional to the strength of the selection force?



#### New idea!   a prevalence map for each resistance factor. Those maps get smoothed during dispersal and amplified during survival/selection
  # these factors can influence survival to plague as an addidive or non-additive relationship... 


# we can make hypotheses about the additivity, about the spatial distribution of the factors and prevalence etc.  
# can some factors be lost to the population if they are very low frequency?

# this could be testable, or more testable that is... 



###### from discussion on 8/5/16

# resistance to plague reduces the transmission of plague. How does that affect this model?
# is there some threshold by which transmission among colonies is notably diluted? At what level do you need resistance before it affects transmission??
# what are the effects of the reduction of transmission with resistance on the evolution of resistance itself?

# no such thing as chronic disease in the plague world. 
# should include a plague reservoir (other rodent species, pooled)
#   every once in a while plague jumps through the reservoir
#   look at avialable literature on plague in other rodents- and the env drivers especially
# pdogs serve as plague amplifiers but not the primary reservoir
#    given that, how do pdogs really affect the plague cycle?
#    can we look at potential resistance within the reservoir (and the type of resistance?)
#    
# stay away from enzootic vs epizootic (but I'm not so sure this isn't interesting)

# what level of resistance is necessary to allow for large pdog complex?

# [my thought] what about evolution of social behaviors? lower densities? Decreased tendency to form large contiguous colonies?

# this is a good model system because pdogs were truly naive to plague, and near 100% susceptible

# can we think about getting genome scan for reservoir species?

# two phases of transmission- among individuals and among colonies. we focus on the second

## Notes from Katie...

# resistance as a mixed trait- encourages endemicity? Some of the population is always susceptible and capable of transmitting the disease. 
# some of the population is resistant and incapable of transmitting the disease. Why encourage endemicity? Because some susceptible individuals are
# always present

# It will be difficult to incorporate the change in transmission with resistance because we don't have the data

# It should be fine to think about transmission as the probability of getting plague (especially as a function of pop density, size and number of 
    # surrounding infected populations)

# Tonie really wants to focus on coevolution, and has historical plague samples

# I think the stability of plague makes this a good model system for studying the host evolution- isolating that as the main evolutionary process

# Plague is a good model system because it is an inherently virulent disease. Novel emerging diseases are likely to exhibit extreme virulence and 
#  influence host evolution in a similar way

# Focus on novel vector-borne diseases. Vector-borne diseases are likely to become ever more prevalent and new diseases are likely to keep emerging
# in a changing world.    (focus on climate change)

# Plague is a good model system because we have data for a range of sites in which hosts are naive and have had many generations with which to evolve. 


# still, it seems that including data on the asian system would be very instructive... 


# Our modeling framework integrates simulation and statistical modeling. An appropriate way to model this system

# hot epidemic with spatial "popping"- this is typical of a virulent disease- and the resulting intermittent selection pressures could 
#   have a huge impact on host resistance evolution







