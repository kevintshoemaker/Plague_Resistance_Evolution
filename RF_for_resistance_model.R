###########################################
#PDOG DATA FOR RESISTANCE MODEL		#
###########################################

rm(list=ls())

NYEARS <- 50

############
## SET UP
############


library(rgeos)
library(rgdal)
library(raster)
library(gdistance)
library(stringr)

KEVIN_LAPTOP <- FALSE #  
KEVIN_OFFICEPC <- FALSE # 
KEVIN_LAPTOP2 <- TRUE #  

if(KEVIN_LAPTOP) GIT_DIR <- "C:\\Users\\Kevin\\GIT\\Plague_Resistance_Evolution"
if(KEVIN_OFFICEPC) GIT_DIR <- "E:\\GIT\\Plague_Resistance_Evolution"
if(KEVIN_LAPTOP2) GIT_DIR <- "C:\\Users\\KevinT_Kevin\\Documents\\GitHub\\Plague_Resistance_Evolution"

setwd(GIT_DIR)
source("PlagueResistanceEvolution_FUNCTIONS.R")

dirs <- SetUpDirectories()

rm(GIT_DIR)

###################################
##########  GET PLAGUE MODEL    # DONE

setwd(dirs$plaguemod$DataDir)

load("randomForestModel_2017-05-12.RData")    # load rfp1    # DONE 

#Inspect object
rfp1           #Random Forest (party) object

##################################
##########  READ-IN ENVIRONMENTAL COVARIATES     # DONE

# Environmental covariates for Pawnee NG (just a small snippet)

setwd(dirs$plaguemod$CovDir)
load("EnvCov_smallPawnee.RData")

##### DEFINE ENV COVARIATES

plot=F

if(plot){
  plot(ENV_COVARS$lat.c)
  plot(ENV_COVARS$long.c)
  plot(ENV_COVARS$NED.c)     ## elevation
  plot(ENV_COVARS$prcp.SumFal)   # time series of maps
  plot(ENV_COVARS$prcp.WinSpr)   # time series
  plot(ENV_COVARS$prcp.year)     # ts
  plot(ENV_COVARS$sand0.c)     #??
  plot(ENV_COVARS$sand2.c)
  plot(ENV_COVARS$slope.c)
  plot(ENV_COVARS$tmax)   # ts  
} 


###################
# CONVERT TO SINGLE OBJECT (if needed)
###################

convert=F 

if(convert){
  ENV_COVARS <- list()
  
  ENV_COVARS$lat.c <- lat.c
  ENV_COVARS$long.c <- long.c
  ENV_COVARS$NED.c <- NED.c    
  ENV_COVARS$prcp.SumFal <- prcp.SumFal
  ENV_COVARS$prcp.WinSpr <- prcp.WinSpr
  ENV_COVARS$prcp.year <- prcp.year
  ENV_COVARS$sand0.c <- sand0.c
  ENV_COVARS$sand2.c <- sand2.c
  ENV_COVARS$slope.c <- slope.c
  ENV_COVARS$tmax <- tmax
  
  setwd(dirs$plaguemod$CovDir)
  save(ENV_COVARS, file="EnvCov_smallPawnee.RData")   # was EnvCov_stacks_smallPawnee.RData
  
  rm(lat.c,long.c,NED.c,prcp.SumFal,prcp.WinSpr,prcp.year,sand0.c,sand2.c,slope.c,tmax)
}


#######
# MAKE TEMPLATE    # DONE
#######

## create template    ## DONE
r <- raster::raster()   # template to create all other rasters for landscape of interest
extent(r) <- raster::extent(ENV_COVARS$lat.c) #Using smaller extent from covariates (snippet of Pawnee)
raster::res(r) <- c(100,100)



#######
# LOAD COLONY BOUNDARIES (max colonizable area)     ## DONE
#######

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
mcP <- raster::crop(master.colony.Pawne, raster::extent(r))   # do I need this?

patchRaster <- raster::rasterize(mcP, r, background=NA, field=1)

  # patchIDRaster <- raster::rasterize(mcP, r, background=NA, field=mcP@data[,])

plot(patchRaster)

#######
# FOR TESTING ONLY
#######

  ### number of years of colony data

site.years <- as.numeric(substr(names(master.colony.Pawne@data)[grep("y",names(master.colony.Pawne@data))],2,5))


################
  ### loading year-to-year shapefiles and turn into rasters

setwd(dirs$plaguemod$DataDir)

fgdb <- "ALL_shp_1.gdb"
fgdb_list <- sort(ogrListLayers(fgdb))

#List of shapefiles for given site
fgdb_list_site <- fgdb_list[substr(fgdb_list, 1, 5) == "Pawne"]   # get all pawnee colonies

Pawn <- stack(r) 
years <- seq(from=min(site.years), to=max(site.years)) 
i=years[1]
for (i in years){
  temp.yrs <- fgdb_list_site[which(str_sub(fgdb_list_site, -4,-1)==i)]
  rast.yrs <- stack(r)
  shp = temp.yrs[1]
  for(shp in temp.yrs){
    temp <- readOGR(dsn=fgdb, layer=shp)
    shp.fix <- gBuffer(temp, byid=TRUE, width=0)
    rast.temp <- rasterize(shp.fix, r)
    rast.temp[is.na(rast.temp)] <- 0 #Set NAs to 0 (useful later on)
    rast.yrs <- stack(rast.yrs, rast.temp)
    rast.yrs <- stackApply(rast.yrs, 1, sum)
  }
  Pawn <- stack(Pawn, rast.yrs)
}
names(Pawn) <- paste("P", years, sep="")
values(Pawn) <- ifelse(values(Pawn)>0,1,0)

plot(Pawn)


##################
#### use the year-year rasters to determine switches in occupancy status year to year (one fewer year than colony maps)

#For each year, we need to know if a pixel was active/inactive the year before.
status <- stack(r)
i=2
for (i in 2:nlayers(Pawn)){
  #Status change: subtract the previous year from the first year: +2 means the pixel became active, -1 means pixel became inactive
  chng <- Pawn[[i]] - Pawn[[i-1]]
  values(chng) <- ifelse(values(chng)>0, 2, ifelse(values(chng)<0, -1, 0))
  
  #Status maintained: multiply the two years together: 1 means active in both years
  maint <- Pawn[[i]] * Pawn[[i-1]]
  
  #Add together
  status.temp <- chng + maint
  
  #Change 0s to NAs
  values(status.temp) <- ifelse(values(status.temp)==0, NA, values(status.temp))
  
  #Report for each year
  status <- stack(status, status.temp)
}
names(status) <- paste("Stat", years[2:length(years)], sep="")


plot(status)


### determine plague status of each pixes on the basis of massive dieoffs and prev colony size from master shapefile

#If a pixel becomes inactive, was it attributable to plague?  

plagued <- stack(r)
for (i in 1:(nlayers(Pawn)-1)){  #loop through the years
  #Find the colonies that plagued out in each year and rasterize them
  col.numb <- i+(ncol(mcP@data)-(length(site.years)-1))  #"p" columns      
  temp.rows <- row(mcP@data)
  temp <- mcP@data[,col.numb]==1
  x <- temp.rows[which(temp==TRUE)]
  z <- status[[i]] #A template raster for later
  
  #If we actually have plagued out colonies...
  if(!is.na(x[1])){
    plagued.t <- mcP[x,]
    plagued.r <- rasterize(plagued.t, r, field=plagued.t@data[,col.numb])
    plagued.r[is.na(plagued.r)] <- 0
    #...Find the pixels that became inactive and were within a plagued out colony
    values(z) <- ifelse(!is.na(values(z)), 0, NA)
    values(z) <- ifelse(values(plagued.r)==1 & values(status[[i]])==-1, 1, values(z)) }
  #ifelse(!is.na(values(status[[i]])) & values(z)==1, 0, NA)) }
  
  #...otherwise, set all colony areas to 0
  if(is.na(x[1])){
    values(z) <- ifelse(!is.na(values(status[[i]])), 0, NA)}
  
  plagued <- stack(plagued, z)
}
#Now each plagued layer is a raster of pixels that had plague in each year (2:nyears) (value=1 for plague, value=0 for no plague)
names(plagued) <- paste("Plagued", years[2:length(years)], sep="")

plot(plagued)


#######
# COMPUTE INTERNAL COVARIATES/STATE VARIABLES    (for now, just a snippet of pawnee data)
#######

####INTERNAL COVARIATES


#OCCUPIED AREA using the master colony dataset, using the previous year's size as the raster field (NOT USED RIGHT NOW)

area <- stack(r)
for (i in 1:(nlayers(Pawn)-1)){
	temp <- rasterize(mcP, r, field=mcP@data[,i+1])   # use "area" field in master colony data
	mask.c <- status[[i]]
	temp.mask <- mask(temp, mask.c)
	area <- stack(area, temp.mask)
}
names(area) <- paste("area", years[2:length(years)], sep="")

plot(area)


#PATCH SIZE of previous year (this is the variable currently used in the model, and not occupied area)
psize <- stack(r)
i=1
for (i in 1:(nlayers(Pawn)-1)){
	temp.patch <- clump(Pawn[[i]], directions=8, gaps=FALSE)
	cells.patch <- freq(temp.patch)     # number of patches in clump
	cells.patch <- as.data.frame(cells.patch[1:(nrow(cells.patch)-1),,drop=F])
	temp.psize <- subs(temp.patch, cells.patch, by=1, which=2)  # number of cells in patch, for each pixel
	psize <- stack(psize, temp.psize)
}
names(psize) <- paste("psize", years[2:length(years)], sep="")

plot(psize)

#YEARS SINCE LAST PLAGUE
years.plague <- Pawn[[1]]
values(years.plague) <- NA
for (i in 2:(nlayers(Pawn))){
	temp.reset <- plagued[[i-1]]
	temp.reset[is.na(temp.reset)] <- 0
	temp.yrs <- years.plague[[i-1]]
	values(temp.yrs) <- ifelse(values(temp.reset)==1, 0, ifelse(is.na(values(temp.yrs)), NA, (values(temp.yrs) + 1)))
	years.plague <- stack(years.plague, temp.yrs)
}
years.plague <- years.plague[[2:length(years)]]
names(years.plague) <- paste("yrs.plg", years[2:length(years)], sep="")

plot(years.plague)


#SLOPE COST TO PLAGUE
slope.trans <- transition(ENV_COVARS$slope.c, mean, 8)   # cost to move to neighboring cells
slope.cost.plague <- stack(r)
dist.noplague <- Pawn[[1]]
values(dist.noplague) <- 50000    # 500000? #Need a raster for distance if there's no plague in the area (put distance at 500km)
for (i in 2:(nlayers(Pawn)-1)){
	#If there's no plague, give a really high cost-distance value (same as dist.noplague)
	temp <- plagued[[i-1]]
	values(temp) <- ifelse(values(temp)==0, NA, values(temp))
	if(any(!is.na(values(temp))))
	{plague.pts <- rasterToPoints(plagued[[i-1]], fun=function(x){x==1})[,1:2]  ##Need plagued areas as points
	slope.cost <- accCost(slope.trans, plague.pts)
	slope.cost <- resample(slope.cost, Pawn[[1]])} else {slope.cost <- dist.noplague}
	slope.cost.plague <- stack(slope.cost.plague, slope.cost)
}
names(slope.cost.plague) <- paste("slopecost", years[3:length(years)], sep="")


plot(slope.cost.plague)


#######
#Creating a dataframe will show you how the year matching works (but not all of the items in this df are read in with this script)

#######
#Create a dataframe/array where each row is a raster pixel, each column is an attribute (e.g., status, colony size), each dim+1 is a year
df <- array(NA, dim=c(dim(Pawn)[1]*dim(Pawn)[2], 30, dim(Pawn)[3]-2))
i=1
for (i in 1:(dim(Pawn)[3]-2)){
	df[,1,i] <- values(plagued[[i+1]])
	df[,2,i] <- values(status[[i+1]]) #Add also previous year's status as a predictor?
	df[,3,i] <- values(first[[i+1]])
	df[,4,i] <- values(area[[i+1]]) #Occupied area of the colony in previous year
	df[,5,i] <- values(psize[[i+1]]) #Patch size of previous year
	df[,6,i] <- values(n.neighbors[[i+1]])
	df[,7,i] <- values(n.neigh.plague[[i]])
	df[,8,i] <- values(dist.plague[[i]])
	df[,9,i] <- values(slope.cost.plague[[i]])
	df[,10,i] <- values(prcp.WinSpr[[i]]) #Winter/Spring precipitation of 2 previous (t-2) year
	df[,11,i] <- values(prcp.WinSpr[[i+1]]) #Winter/Spring precip of previous year
	df[,12,i] <- values(prcp.WinSpr[[i+2]]) #Winter/Spring precip of current year
	df[,13,i] <- values(prcp.SumFal[[i]]) #Summer/Fall precipitation 2 previous (t-2) year
	df[,14,i] <- values(prcp.SumFal[[i+1]]) #Summer/Fall precip of previous year
	df[,15,i] <- values(prcp.SumFal[[i+2]]) #Summer/Fall precip of current year
	df[,16,i] <- values(prcp.year[[i]]) #Whole year precip of 2 previous (t-2) year (Jan-Dec)
	df[,17,i] <- values(prcp.year[[i+1]]) #Whole year precip of previous year (Jan-Dec)
	df[,18,i] <- values(prcp.WinSpr[[i+2]] + prcp.SumFal[[i+2]]) #Whole year precip of current year (Jan-Sep)
	df[,19,i] <- values(tmax[[i]]) #Temp of 2 previous (t-2) year
	df[,20,i] <- values(tmax[[i+1]]) #Temp of previous year
	df[,21,i] <- values(tmax[[i+2]]) #Temp of current year
	df[,22,i] <- values(sand0.c)
	df[,23,i] <- values(sand2.c)
	df[,24,i] <- values(NED.c)
	df[,25,i] <- values(lat.c)
	df[,26,i] <- values(long.c)
	df[,27,i] <- values(years.plague[[i]])
	df[,28,i] <- values(year[[i+2]])
	df[,29,i] <- values(poison[[i+2]])
	df[,30,i] <- values(colonyID) #Colony ID needs to be last variable
}
#All the "non-active" cells are NAs




#######
# OLD CODE




























#