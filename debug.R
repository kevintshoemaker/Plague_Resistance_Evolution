
############
## CLEAR WORKSPACE
############

rm(list=ls())


############
## SIMULATION CONTROLS
############

NYEARS <- 25

############
## SET GLOBAL VARS
############

KEVIN_LAPTOP <- TRUE #  FALSE #  
KEVIN_OFFICEPC <- FALSE # TRUE # 

if(KEVIN_LAPTOP) GIT_DIR <- "C:\\Users\\Kevin\\GIT\\Plague_Resistance_Evolution"
if(KEVIN_OFFICEPC) GIT_DIR <- "E:\\GIT\\Plague_Resistance_Evolution"

#####################
# LOAD FUNCTIONS
#####################

setwd(GIT_DIR)
source("PlagueResistanceEvolution_FUNCTIONS.R")


############
## SET UP WORKSPACE AND LOAD PACKAGES
############

dirs <- SetUpDirectories()
num_cores <- detectCores() - 1   # for setting up cluster... leave one core free for windows background processes?

############
## SAMPLE FROM LATIN HYPERCUBE
############

N_LHS_SAMPLES <- 20

masterDF <- MakeLHSSamples(nicheBreadthDir=dir,NicheBreadth)

###########
##  START A PARALLEL FOR LOOP
###########

cl <- makeCluster(num_cores,outfile="LOG.TXT")
registerDoParallel(cl=cl)    # make the cluster


MakeWorker(NYEARS, masterDF, dirs)(1)






