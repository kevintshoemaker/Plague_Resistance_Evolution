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
## SIMULATION CONTROLS
############

NYEARS <- 20

############
## SET BASE DIRECTORY
############

KEVIN_LAPTOP <- FALSE #  TRUE #  
KEVIN_OFFICEPC <- TRUE # FALSE # 

if(KEVIN_LAPTOP) BASE_DIR <- "C:\\Users\\Kevin\\Dropbox\\PlagueModeling\\ResistanceEvolution"
if(KEVIN_OFFICEPC) BASE_DIR <- "E:\\Dropbox\\PlagueModeling\\ResistanceEvolution"

if(KEVIN_LAPTOP) GIT_DIR <- "C:\\Users\\Kevin\\GIT\\Plague_Resistance_Evolution"
if(KEVIN_OFFICEPC) GIT_DIR <- "E:\\GIT\\Plague_Resistance_Evolution"

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

############
## USER-DEFINED VARIABLES
############

UserParams <- DefineUserParams()

#####################
# INITIALIZE DISPERSAL   (for both plague and no plague... )   
#####################

InitializeDispersal(UserParams)

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
######  INITIALIZE THE SYSTEM!
#########################


########################
# INITIALIZE POPULATION
########################

InitDensRaster <- KRaster    # initialize abundance at carrying capacity
plot(InitDensRaster)


#######################
# INITIALIZE ALLELE FREQUENCIES / RESISTANCE FACTORS [keep for now- will be multiple genes in the model somehow]
#######################
# NOTE: Some regions are more likely to evolve faster because they have greater percentages of those genes that can confer resistance. 

InitFreqList <- GetInitFreqs(UserParams)

plot(InitFreqList)


#####################
# INITIALIZE POPULATION
#####################


### Code block for pop starting from small loci

InitDensRaster2 <- reclassify(patchIDRaster,rcl=c(-Inf,Inf,0))    # for testing
ndx <- sample(which(!is.na(InitDensRaster2@data@values)),size=3)
InitDensRaster2[ndx] <- 1000   # initialize population in random locations
InitDensRaster <- InitDensRaster2
#PopArray2 <- InitDensRaster   # copy, for dispersal algorithm... 


PopArray <- GetStructuredPop(InitDensRaster)

plot(PopArray)
  

######################
# INITIALIZE PLAGUE PROCESS   [KTS: moving away from this and towards a purely statistical model]
######################
#  for now, assume that plague hits at the patch level, and is a random process.

# PROB_PLAGUE_YEAR <- 0.5 # probability that a plague event hits (landscape level?)   # for now, plague only hits one patch in a plague year
# 
# plagueyear = as.logical(rbinom(NYEARS,1,PROB_PLAGUE_YEAR))
# 
# plagueNow = floor(runif(NYEARS,1,nPatches+1)) * as.numeric(plagueyear)    ## which patch plagues out?

PlagueRaster_template <- reclassify(patchIDRaster,rcl=c(-Inf,Inf,0))   

PlagueRaster <- doPlague(UserParams,PlagueRaster=PlagueRaster_template, PopArray=reclassify(patchIDRaster,rcl=c(-Inf,Inf,0)))


plot(PlagueRaster)



####################
# START LOOP THROUGH YEARS
####################


    # names of important raster maps to save to file etc...
#rasterNames  <- c("PopArray","NextPlagueSurvRaster","NextNormalSurvRaster","PlagueResistancePotentialRaster")   # Deprecate?


# t=which(plagueyear)[1]
t=0
t=t+1
for(t in 1:(NYEARS)){
  deviate <- rnorm(1)   #determine if this is a good year or a bad year (for now, survival and fecundity are perfectly correlated)
  cv=UserParams$Popbio$CV_SURVIVAL   # set up for using the getYearVariate function
  
  if(t==1){ FreqList=InitFreqList; DensRaster=InitDensRaster; newFociRaster <- reclassify(PopArray,rcl=c(-Inf,Inf,0)) }    # initial conditions
  
  ##################
  # DENSITY INDEPENDENT SURVIVAL (including plague survival)
  ##################
  
  #PopArray <- GetStructuredPop(DensRaster)
  
  DensRaster <- doSurvival(UserParams,DensRaster,PlagueRaster,FreqList)
  # plot(DensRaster)
  
  ################
  # REPRODUCTION
  ################
  
  DensRaster <- doReproduce(UserParams,DensRaster,PlagueRaster)    # TODO: make specific to each resistance type...?
  # plot(DensRaster)
  
  
  ###############
  # DISPERSAL: Move individuals around the landscape (this takes a while!)
  ###############
  PopArray <- doDispersal(UserParams,PlagueRaster)
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




########## DEPRECATED CODE

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

# ####################
# # INITIALIZE SURVIVAL  (deprecate?)
# ####################
# 
# meansurv <- matrix(0,nrow=2,ncol=2)    # survival matrix (mean)
# rownames(meansurv) <- c("resistant","susceptible")
# colnames(meansurv) <- c("plague","noPlague")
# meansurv["resistant","noPlague"] <- getSurvival("resistant","noPlague")
# meansurv["resistant","plague"] <-   getSurvival("resistant","plague")
# meansurv["susceptible","plague"] <- getSurvival("susceptible","plague")
# meansurv["susceptible","noPlague"] <- getSurvival("susceptible","noPlague")
# meansurv

