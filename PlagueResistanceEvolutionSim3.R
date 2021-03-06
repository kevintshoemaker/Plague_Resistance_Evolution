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
## SIMULATION CONTROLS
############

NYEARS <- 56    # 1980 to 2036

############
## SET GLOBAL VARS
############

KEVIN_LAPTOP <- FALSE #  
KEVIN_OFFICEPC <- FALSE # FALSE #  TRUE # 
KEVIN_LAPTOP2 <- FALSE # TRUE #
KEVIN_LAPTOP3 <- TRUE #
MIRANDA <- FALSE # TRUE

if(KEVIN_LAPTOP) GIT_DIR <- "C:\\Users\\Kevin\\GIT\\Plague_Resistance_Evolution"
if(KEVIN_OFFICEPC) GIT_DIR <- "E:\\GIT\\Plague_Resistance_Evolution"
if(KEVIN_LAPTOP2) GIT_DIR <- "C:\\Users\\KevinT_Kevin\\Documents\\GitHub\\Plague_Resistance_Evolution"
if(KEVIN_LAPTOP3) GIT_DIR <- "C:\\Users\\Kevin\\Documents\\GitHub\\Plague_Resistance_Evolution"
if(MIRANDA) GIT_DIR <- "C:\\Users\\kshoemaker\\Documents\\GitHub\\Plague_Resistance_Evolution"

#####################
# LOAD FUNCTIONS
#####################

setwd(GIT_DIR)
source("PlagueResistanceEvolution_FUNCTIONS.R")


############
## SET UP WORKSPACE AND LOAD PACKAGES
############

dirs <- SetUpDirectories()

rm(GIT_DIR)

num_cores <- parallel::detectCores() - 2   # for setting up cluster... leave one core free for windows background processes?

#num_cores <- 10


############
## SAMPLE FROM LATIN HYPERCUBE
############

N_LHS_SAMPLES <- 30  # 200

masterDF <- MakeLHSSamples(add=FALSE)    # add=FALSE
 
#nrow(masterDF)

rep=1
fake=T


oops <- FALSE    # if needed, recover the previous master dataframe...
if(oops){
  setwd(dirs$DATA_DIR)
  masterDF2 <- read.csv("masterDF_prelim2017-07-14.csv")   # masterDF_prelim2017-07-14.csv
  masterDF <- masterDF2
  rm(masterDF2)
}


########### Read in (simulated) genetics data

all_loci <- read.table("allele-freqs_10pops.txt",header=TRUE)
key_markers <- paste("Marker",c(10,32,43,50,88,19,81,96,115,127,6,3,48,102,111),sep="")
bg <- sample(setdiff(colnames(all_loci[,-c(1,2)]),key_markers),5,replace=FALSE)
resistance_fingerprint <- all_loci[all_loci$Allele==2,c(key_markers,bg)]
row.names(resistance_fingerprint) <- all_loci$pop[as.numeric(row.names(resistance_fingerprint))]
resistance_fingerprint


###############
# code the rules for resistance

rules_for_resistance <- list()
rules_for_resistance[[1]] <- list()
rules_for_resistance[[1]]$loci <- paste("Marker",c(10,32,43,50,88),sep="")
rules_for_resistance[[1]]$rules <- do.call(sprintf,c(list("1-((1-%s)*(1-%s)*(1-%s)*(1-%s)*(1-%s))"),as.list(rules_for_resistance[[1]]$loci)))

rules_for_resistance[[2]] <- list()
rules_for_resistance[[2]]$loci <- paste("Marker",c(6,19,81,96,115,127),sep="")
rules_for_resistance[[2]]$rules <- do.call(sprintf,c(list("(%s)*(1-((1-%s)*(1-%s)*(1-%s)*(1-%s)*(1-%s)))"),as.list(rules_for_resistance[[2]]$loci)))

rules_for_resistance[[3]] <- list()
rules_for_resistance[[3]]$loci <- paste("Marker",c(3,48,102,111),sep="")
rules_for_resistance[[3]]$rules <- do.call(sprintf,c(list("(%s*%s*%s*%s)"),as.list(rules_for_resistance[[3]]$loci)))

rules_for_resistance

###############
# code the dominance structure

temp <- matrix(0,nrow=ncol(resistance_fingerprint),ncol=3)
colnames(temp) <- c("2x(rr)","1x(rs)","0x(ss)")
rownames(temp) <- colnames(resistance_fingerprint)
temp

atleastone <- paste("Marker",c(3,48,102,111,6),sep="") 
needstwo <- setdiff(colnames(resistance_fingerprint),atleastone)
alo_ndx <- match(atleastone,colnames(resistance_fingerprint))
nt_ndx <- match(needstwo,colnames(resistance_fingerprint))

temp[alo_ndx,] <- matrix(rep(c(1,1,0),times=length(alo_ndx)),ncol=3,byrow = T)
temp[nt_ndx,] <- matrix(rep(c(1,0,0),times=length(nt_ndx)),ncol=3,byrow = T)

dmat <- temp

# Here are some rulesets that guided my thinking:
#   1) If an individual is homozygous for allele 2 at locus 10, 32, 43, 50, or 88, that individual will be resistant
# 2) If an individual is homozygous for allele 2 at locus 19, 81, 96, 115, or 127 and has at least one copy of allele 2 at locus 6, that individual will be resistant
# 3) If an individual has at least one copy of allele 2 at locus 3, 48 102 and 111, that individual will be resistant


## note masterdf is written to data directory

# setwd(dirs$DATA_DIR)
# masterDF <- read.csv("masterDF_prelim2017-01-26.csv",header=T)

###########
##  START A PARALLEL FOR LOOP
###########

library(parallel)
library(doParallel)

cl <- parallel::makeCluster(num_cores,outfile="LOG.TXT")
doParallel::registerDoParallel(cl=cl)    # make the cluster


#######################
## objects to export to each node in the cluster

# functionlist <- c()   # , 'mp.write'
# filelist <- c()  # 'masterDF', 'MP_DIRECTORY','template','GENTIME','humanArrival.df','EXE_DIRECTORY','DLL_FILENAME','dispersalFunc.df','DistClasses','NPOPS','DistBins',
# 
# objectlist <- c(functionlist,filelist)   # full list of objects to export
# 
# 
# packagelist <- c("secr","igraph","raster")

allsamples <- foreach(i = 101: nrow(masterDF)    # (300+1): nrow(masterDF)
                      # .export=objectlist,
                      # .packages = packagelist,
                      # .errorhandling=c("pass")
) %dopar% {   
  
  # #####################
  # # LOAD FUNCTIONS
  # #####################
  # 
  # setwd(GIT_DIR)
  # source("PlagueResistanceEvolution_FUNCTIONS.R")
  # 
  
  ############
  ## SET UP WORKSPACE AND LOAD PACKAGES
  ############
  
  #SetUpWorkspace()
  #num_cores <- detectCores() - 1   # for setting up cluster... leave one core free for windows background processes?
  
  temp <- MakeWorker(NYEARS, masterDF, dirs, fake)(i)  #DoSimulateResistance(rep=i)    # simulate for these params... 
  
}     ## end parallel for loop


#MakeWorker(NYEARS, masterDF, dirs)(1)

###################
# CLOSE CLUSTER
###################

if(!is.null(cl)) {
  parallel::stopCluster(cl)
  cl <- c()
}


?rgb


########################
# MAKE MOVIES POST HOC (if necessary)
########################

makemov <- FALSE

if(makemov){
  rep=2
  for(rep in 1:nrow(masterDF)){
    if(rep%%2==0) MakeMovie(rep)
  }
}

########################
# HARVEST DATA INTO A CONVENIENT FORMAT FOR ANALYSIS
########################

masterDF2 <- HarvestData(masterDF,dirs)

#masterDF2[131,]

########################
# ANALYZE DATA
########################

### use random forest?

#### Read in the script from github

source_github("https://raw.githubusercontent.com/kevintshoemaker/Random-Forest-Functions/master/","RF_Extensions.R")

#library(coin)
library(party)

#Define catagorical variables as factors


masterDF2 <- na.omit(masterDF2)

masterDF2$ISRES <- as.logical(masterDF2$ISRES)    # as.factor(masterDF2$ISRES)
masterDF2$ISEXT <- as.logical(masterDF2$ISEXT) # as.factor

masterDF2$DOMINANCE <- as.factor(masterDF2$DOMINANCE)


df <- na.omit(masterDF2)

############### NAMING VARIABLES ############

predictorNames <- c(  "% Habitat",       # nice readable names
                      "Conspecific attraction",
                      "Plague colony disruption",
                      "Maximum density per ha",
                      "Plague survival, susceptible",
                      "Mean fecundity",
                      "Fitness cost of resistance",
                      "Initial resistance allele freq.",
                      "Dominance of resistance factors"
)

pred.names=c(  "PER_SUITABLE",      
               "SNUGGLE",
               "PLAGUE_DISPERSAL",
               "MAXDENS",
               "BASELINE_PLAGUESURV",
               "BASELINE_MEANFEC",
               "FITNESS_COST",
               "INITFREQ",
               "DOMINANCE"
)


#name check
cbind(pred.names,predictorNames)


#### Define response variable

response=  "ISRES"    # "ISEXT"   #   "ISRES"    # 

#### Define our formula (response ~ predictors)

formula1 <- as.formula(paste(response,"~",paste(pred.names,collapse="+")))




#source("C:\\Users\\Kevin\\GIT\\Random-Forest-Functions\\RF_Extensions.R")   # change to your script locations

##### CONDITIONAL INFERENCE TREE  ##################

res.tr <- ctree(formula=formula1, data=df, controls = ctree_control(mincriterion = 0.85,maxdepth = 3))

graphics.off()
plot(res.tr)

summary(res.tr)

###########################################################
###############  CFOREST #################

cforestControl <- cforest_unbiased(ntree=1000,mtry=5)   # change back to 500!!
cforestControl@fraction <- 0.6

cforestControl@gtctrl@mincriterion <- 0.85

rf_model1 <- cforest(formula1, controls=cforestControl, data=df)

# get the importance values
model1_importance<-varimp((rf_model1), conditional= FALSE)

graphics.off()
lengthndx <- length(model1_importance)
#par(mai=c(0.95,3.1,0.6,0.4))
par(mai=c(1.4,3.4,0.6,0.9))
col <- rainbow(lengthndx, start = 3/6, end = 4/6)      # rep(brewer.pal(6,"Blues"),each=2)
barplot(height=model1_importance[order(model1_importance,decreasing = FALSE)],
        horiz=T,las=1,main="Order of Importance of Predictor Variables",
        xlab="Index of overall importance",col=col,           
        names.arg=predictorNames[match(names(model1_importance),pred.names)][order(model1_importance,decreasing = FALSE)])



##### Make univariate plots of the relationships- plot all relationships at once

graphics.off()
RF_UnivariatePlots(object=rf_model1, varimp=model1_importance, data=df,  #   
                   predictors=pred.names, labels=predictorNames, allpredictors=pred.names,plot.layout=c(2,2))



##### Make univariate plots of the relationships- plot one relationship at a time

graphics.off()
RF_UnivariatePlots(object=rf_model1, varimp=model1_importance, data=df,  #   
                   predictors=pred.names[7], labels=predictorNames[7], allpredictors=pred.names,plot.layout=c(1,1))


# return the data for plotting
PlotData <- RF_UnivariatePlots(object=rf_model1, varimp=model1_importance, data=df,  #   
                               predictors=pred.names, labels=predictorNames, allpredictors=pred.names, plot.layout=c(1,1),plot=F)



####################################
#######################   RANDOM FOREST FIND AND PLOT INTERACTIONS

# NOTE: this one can take a very long time   ...
rf_findint <- RF_FindInteractions(object=rf_model1,data=df,predictors=pred.names)

# display and plot out interactions...
rf_findint$interactions1

rf_findint$rank.list1

### plot interaction strength
graphics.off()
lengthndx <- min(9,nrow(rf_findint$rank.list1))
par(mai=c(0.95,3.1,0.6,0.4))
#ndx <- ndx <- which(predictors%in%pred.names)
barplot(height=(rf_findint$rank.list1[c(1:min(9,nrow(rf_findint$rank.list1))),5][c(lengthndx:1)]),
        horiz=T,las=1,main=paste(response, sep=""),
        xlab="Index of interaction strength",col=brewer.pal(lengthndx,"Blues"),           
        names.arg=paste("",predictorNames[match(rf_findint$rank.list1[,2][c(lengthndx:1)],pred.names)],"\n",predictorNames[match(rf_findint$rank.list1[,4][c(lengthndx:1)],pred.names)],sep="") )

graphics.off()


rf_findint$rank.list1



fam="gaussian"
graphics.off()
#svg(filename = "IntFig2.svg",
# width = 7, height = 7, pointsize = 12,
# onefile = TRUE, family = "sans", bg = "white")

#### visualize the interactions

RF_InteractionPlots(x=3,y=8,object=rf_model1,data=df,predictors=pred.names,family=fam) 

RF_InteractionPlots(x=9,y=1,object=rf_model1,data=df,predictors=pred.names,family=fam)

RF_InteractionPlots(x=3,y=5,object=rf_model1,data=df,predictors=pred.names,family=fam)

dev.off()
graphics.off()





###################################
#################### CROSS VALIDATION CODE

n.folds = 10
foldVector = rep(c(1:n.folds),times=floor(length(newX_listed$Listed)/9))[1:length(newX_listed$Listed)]
#n.folds = length(newX_listed$resp_factor)
#foldVector <- c(1:length(newX_listed$resp_factor))

counter = 1
CVprediction <- numeric(nrow(newX_listed))
CVobserved <- numeric(nrow(newX_listed))
realprediction <- numeric(nrow(newX_listed))
realdata <- numeric(nrow(newX_listed))

predictCols <- which(names(newX_listed)%in%pred.names)

data.controls = cforest_unbiased(ntree=50)
counter=1
response="Listed"    #"resp_factor"

#test <- numeric(nrow(newX_listed))
for(i in 1:n.folds){
  model <- cforest(formula1, data = newX_listed[which(foldVector!=i),], controls=data.controls) 
  predict_CV  <- predict(model,newdata=newX_listed[which(foldVector==i),],type="prob") 
  predict_real  <-  predict(rf_model1,newdata=newX_listed[which(foldVector==i),],type="prob")
  REAL <- newX_listed$Listed[which(foldVector==i)]
  for(j in 1:length(which(foldVector==i))){
    CVprediction[counter] <- as.numeric(predict_CV[[j]][,2])
    CVobserved[counter] <-  REAL[j]      
    realprediction[counter] <- as.numeric(predict_real[[j]][,2])   
    realdata[counter] <- REAL[j]         
    counter = counter + 1  
  }
}

fact=TRUE
if(fact){
  CVobserved = CVobserved-1
  realdata=realdata-1
}

CV_RMSE = sqrt(mean((CVobserved-CVprediction)^2))       # root mean squared error for holdout samples in 10-fold cross-validation ...
real_RMSE = sqrt(mean((CVobserved-realprediction)^2))  # root mean squared error for residuals from final model

# print RMSE statistics
CV_RMSE 
real_RMSE   

binaryresponse=TRUE

if(binaryresponse){
  graphics.off()
  par(mfrow=c(2,1))
  pred <- prediction(CVprediction,CVobserved)     # for holdout samples in cross-validation
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  plot(perf)
  text(.9,.1,paste("AUC = ",round(auc@y.values[[1]],2),sep=""))
  
  pred <- prediction(realprediction,CVobserved)     # for final model
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  plot(perf)
  text(.9,.1,paste("AUC = ",round(auc@y.values[[1]],2),sep=""))
}

# COHEN KAPPA statistics

graphics.off()
par(mfrow=c(2,1))
thresholds <- seq(0.01,0.99,length=101)   # "artificial" extinction thresholds across which to examine performance
kappa <- numeric(length(thresholds))
for(i in 1:length(thresholds)){
  trueLabels <- CVobserved
  predLabels <- ifelse(CVprediction>=thresholds[i],1,0)
  tot <- length(CVobserved)
  tp <- length(which((trueLabels==1)&(predLabels==1)))  
  tn <- length(which((trueLabels==0)&(predLabels==0)))
  fp <- length(which((trueLabels==0)&(predLabels==1)))
  fn <- length(which((trueLabels==1)&(predLabels==0)))
  pr_agree <- (tp+tn)/tot    # overall agreement, or accuracy
  pr_agree_rand <- ((tp+fn)/tot)*((tp+fp)/tot)+((fn+tn)/tot)*((fp+tn)/tot)
  kappa[i] <- (pr_agree-pr_agree_rand)/(1-pr_agree_rand)
}
plot(thresholds,kappa,type="l",xlab="Threshold", ylab="Cohen's Kappa", main="Holdout sample performance")

# find threshold value associated with highest Kappa for C-V data

cutoff <- thresholds[which.max(kappa)]
cutoff


kappa <- numeric(length(thresholds)) 
for(i in 1:length(thresholds)){
  trueLabels <- CVobserved
  predLabels <- ifelse(realprediction>=thresholds[i],1,0)    
  tot <- length(CVobserved)
  tp <- length(which((trueLabels==1)&(predLabels==1)))  
  tn <- length(which((trueLabels==0)&(predLabels==0)))
  fp <- length(which((trueLabels==0)&(predLabels==1)))
  fn <- length(which((trueLabels==1)&(predLabels==0)))
  pr_agree <- (tp+tn)/tot    # overall agreement, or accuracy
  pr_agree_rand <- ((tp+fn)/tot)*((tp+fp)/tot)+((fn+tn)/tot)*((fp+tn)/tot)
  kappa[i] <- (pr_agree-pr_agree_rand)/(1-pr_agree_rand)
}
plot(thresholds,kappa,type="l",xlab="Threshold", ylab="Cohen's Kappa", main="Performance: full model")



### display confusion matrix and kappa for a single threshold
trueLabels <- CVobserved
predLabels <- ifelse(CVprediction>=cutoff,1,0)    
tot <- length(CVobserved)
tp <- length(which((trueLabels==1)&(predLabels==1)))  
tn <- length(which((trueLabels==0)&(predLabels==0)))
fp <- length(which((trueLabels==0)&(predLabels==1)))
fn <- length(which((trueLabels==1)&(predLabels==0)))
pr_agree <- (tp+tn)/tot    # overall agreement, or accuracy
pr_agree_rand <- ((tp+fn)/tot)*((tp+fp)/tot)+((fn+tn)/tot)*((fp+tn)/tot)
kappa[i] <- (pr_agree-pr_agree_rand)/(1-pr_agree_rand)
kappa[i]
matrix(c(tp,fp,fn,tn),nrow=2,ncol=2)
sensitivity <- tp/(tp+fn)
specificity <- tn/(tn+fp)
toterror <- (fn+fp)/tot
sensitivity
specificity
toterror

if(binaryresponse){
  CVprediction[which(CVprediction==1)] <- 0.9999
  CVprediction[which(CVprediction==0)] <- 0.0001
  realprediction[which(realprediction==1)] <- 0.9999
  realprediction[which(realprediction==0)] <- 0.0001
}


realdata = CVobserved
fit_deviance_CV <- mean((CVobserved-CVprediction)^2)
if(binaryresponse) fit_deviance_CV <- mean(-2*(dbinom(CVobserved,1,CVprediction,log=T)-dbinom(realdata,1,realdata,log=T)))
fit_deviance_real <- mean((CVobserved-realprediction)^2)
if(binaryresponse) fit_deviance_real <- mean(-2*(dbinom(CVobserved,1,realprediction,log=T)-dbinom(realdata,1,realdata,log=T)))
null_deviance <- mean((CVobserved-mean(CVobserved))^2)
if(binaryresponse) null_deviance <- mean(-2*(dbinom(CVobserved,1,mean(CVobserved),log=T)-dbinom(realdata,1,realdata,log=T)))
deviance_explained_CV <- (null_deviance-fit_deviance_CV)/null_deviance   # based on holdout samples
deviance_explained_real <- (null_deviance-fit_deviance_real)/null_deviance   # based on full model...

deviance_explained_CV
deviance_explained_real



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

