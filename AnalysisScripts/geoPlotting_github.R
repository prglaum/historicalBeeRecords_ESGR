###____________________###
##Species mapping - GBIF##
###____________________###

## Contents -----------
## Load libraries/packages - Line 13
## Load necessary data - Line 32
## Create basic species occurance plots - Line 39
## Create centroid using all gbif species entries per species in each major sample period - Line 97
## Create centroid of centroids weighted by species abundance in ESGR - Line 138
## Plot centroids from each centroid creation methodolgy - Line 172

##Load libraries -----------------------
library(stringr) #for string manipulation
library(ggplot2) #plotting
library(readr) #for loading
library(dplyr) #for data wrangling
library(readxl) ##for loading xlsx files
#library(tidyverse)

library(sf) # the main geo-plotting software package...
library(diagram)
library(units) #used for controlling the hull area's units to be km^2
library(moments) #used to calculate skewness and kurtosis
library(nngeo) #calculating nearest neighbor on maps
library(mapview) #For quick plotting maps
library(spData) # this is where the world data comes from. Used to define continents

library(leaflet) #basic leaflet maps, used to plot labeled points.
library(leaflet.extras2) #Needed plug in for extra functionality. Used for arrows connecting points

##Load data ------------------------
##Trait data by species. Also extirpation indicators from 70s to 2010s.
traits=read_excel("~/Grahametal2023_BeeTraits.xlsx",sheet="TraitData",na = "NA")
##Full GBIF occurrence data for our species.
##The data is large and is housed at .... (will add soon)
gbif=readr::read_csv("~/RawGBIFdata.csv")

##BASIC SPECIES OCCURANCE PLOTS (See Fig S2 or Anthophora_terminalis_Example_Map in AnalysisScripts folder).----------------------

##First load the world data
world <- spData::world
##Pull out the continent level shape files
continent=group_by(world, continent)%>%
  summarise(pop=sum(pop, na.rm=T))
##We don't need the population data, remove it
continent=select(continent,-pop)
#Load mapview package if you haven't already
library(mapview)
##plot the continents layer
mapview(continent)

##Choose a species to focus on:
spDF=subset(gbif,gbif$species==traits$species[44]) ##Choose a species to plot

##Some latitude and longitude entries in GBIF are 0.0 (i.e., no data)...scrubbing them here.
Latitudes=spDF$decimalLatitude; Latitudes[Latitudes==0] <- NA
Longitudes=spDF$decimalLongitude; Longitudes[Longitudes==0] <- NA

species.sf <- as.data.frame(cbind(Longitudes=Longitudes,Latitudes=Latitudes)) %>%
  st_as_sf(coords = c("Longitudes", "Latitudes"), crs = 4326)

##Finding nearest neighbor (NN) distance
##Hueristically chosen cutoff distance used to eliminate singletons across oceans.
cutOffDist=10030726;
##
distances <- nngeo::st_nn(species.sf, species.sf, k = 5, returnDist = T);
##Assigning the sum of the 5 nearest neighbor distances to species.sp
species.sf$dst <- distances$dist %>% lapply(., FUN = sum) %>% unlist();
##Only include species inside the cut off distance for outliers from nearest neighbors
species.sf=subset(species.sf,dst<cutOffDist);

species.sf = species.sf %>% st_join(continent);
species.sf = subset(species.sf,!is.na(continent)); #Remove pts not assigned to continent
species.sf=select(species.sf,-c(dst)); ##Remove the dst metric from the sp object

hull <- (subset(species.sf,!is.na(continent))) %>%
  group_by(continent) %>%
  summarise( geometry = st_combine( geometry ) ) %>%
  st_convex_hull() #%>% st_shift_longitude()
hull=st_make_valid(hull);

na.sf=subset(species.sf,continent=="North America")
naCent=na.sf %>% summarize(geometry = st_union(geometry)) %>%  st_centroid

mapview::mapview(list(species.sf,hull)) + mapview(naCent,cex=13,col.regions="red")



### ### ###
##Plotting centroids between 70s', both periods', and 2010s' species ###--------------
### ### ###
##Two methods of sample community centroids are used in Fig S18. They are described
##in the sections below:

###___________________###
##1) Centroids created by GBIF entries for species present in 72/73, 2017/2018, or both periods -----------------
###___________________###
##load necessary data:
# GBIF sample set data. Load if not loaded above.
gbif=readr::read_csv("~RawGBIFdata.csv")
# trait data for ESGR samples. Needed to account for persistence across sample periods. Load if not loaded above.
traits=read_excel("~/Grahametal2023_BeeTraits.xlsx",sheet="TraitData",na = "NA")
##recall, the bin column:
## Bin = 1: Species only found in 72/73
## Bin = 2: Species only found in 2017/2018
## Bin = 3: Species found in both 72/73 & 2017/2018

naCents=rep();
for(i in 1:3) {
	##Focus on each period
	periodSP=subset(traits,bin==i)$species
	gbifSP=subset(gbif,species %in% periodSP)
	##Some latitude entries are 0...scrubbing them here.
	Latitudes=gbifSP$decimalLatitude; Latitudes[Latitudes==0] <- NA
	##Some longitude entries are 0...scrubbing them here.
	Longitudes=gbifSP$decimalLongitude; Longitudes[Longitudes==0] <- NA

	species.sf <- as.data.frame(cbind(Longitudes=Longitudes,Latitudes=Latitudes)) %>%
	  st_as_sf(coords = c("Longitudes", "Latitudes"), crs = 4326) #%>%

	##Assign points to continent
	species.sf = species.sf %>% st_join(continent);
	species.sf = subset(species.sf,!is.na(continent)); #Remove pts not assigned to continent

	##Focus on North America...
	na.sf=subset(species.sf,continent=="North America")
	##Get centroid of GBIF instances for all relevant species
	naCents[i]=na.sf %>% summarize(geometry = st_union(geometry)) %>%  st_centroid
}
centCoorsGBIFsp=as.data.frame(matrix(c(unlist(naCents)),3,2,byrow=TRUE)); #rbind(naCents,naCents,naCents);
colnames(centCoorsGBIFsp)=c("X","Y")
centCoorsGBIFsp$time=c("1970s","2010s","Both"); centCoorsGBIFsp=centCoorsGBIFsp[c(1,3,2),];



###___________________###
##2) Centroid of centroids: Each species centroid is WEIGHTED ----------------------
##by Evas/Isaacs sample size for species present in 72/73, 2017/2018, or both periods
###___________________###

esgrNEW=read.csv("~/ESGR_samples_EvansIsaacs_withTraits.csv",sep=",",header=T, na.strings=c("NaN","#NAME?","-Inf",""))
esgrNEW$GenusSpecies=str_replace_all(esgrNEW$GenusSpecies, c("_" = " " ))
##recall, the bin column:
## Bin = 1: Species only found in 72/73
## Bin = 2: Species only found in 2017/2018
## Bin = 3: Species found in both 72/73 & 2017/2018

naCents=rep();
for(i in 1:3) {
	##Some latitude entries are 0...scrubbing them here.
	#Doing it the esgrNEW way is using instances of centroids from each species, will involve repeats
	centLats=esgrNEW$naCentLat[esgrNEW$bin==i]; centLats[centLats==0] <- NA
	##Some longitude entries are 0...scrubbing them here.
	centLongs=esgrNEW$naCentLong[esgrNEW$bin==i]; centLongs[centLongs==0] <- NA

	species.sf <- as.data.frame(cbind(Longitudes=centLongs,Latitudes=centLats)) %>%
	  st_as_sf(coords = c("Longitudes", "Latitudes"), crs = 4326)
	##Assign points to continent
	species.sf = species.sf %>% st_join(continent);
	species.sf = subset(species.sf,!is.na(continent)); #Remove pts not assigned to continent

	##Focus on North America...
	na.sf=subset(species.sf,continent=="North America")
	naCents[i]=na.sf %>% summarize(geometry = st_union(geometry)) %>%  st_centroid
}
centCoorsESGRsp=as.data.frame(matrix(c(unlist(naCents)),3,2,byrow=TRUE)); #rbind(naCents,naCents,naCents);
colnames(centCoorsESGRsp)=c("X","Y")
centCoorsESGRsp$time=c("1970s","2010s","Both"); centCoorsESGRsp=centCoorsESGRsp[c(1,3,2),];


##PLOT CENTROIDS FROM EACH METHOD ---------------------------------
##Use the leaflet package to plot centroids from each method. Centroids per
##method are connected with arrows
leaflet(naCents[1]) %>%
addProviderTiles("Esri.WorldPhysical") %>%
addCircleMarkers(lng = ~centCoorsESGRsp[,1], lat = ~centCoorsESGRsp[,2],color=c("black"),label = ~(c("1970s","Both","2010s")), #) %>%
  labelOptions = labelOptions(noHide = TRUE, direction = "top",offset=c(0,-12), textsize = "18px", textOnly = FALSE)) %>%
addArrowhead(lng= ~ centCoorsESGRsp[,1], lat= ~ centCoorsESGRsp[,2],color="black") %>% addGraticule(interval=2.5) %>%
#leaflet(naCents[1]) %>%
addProviderTiles("Esri.WorldPhysical") %>%
addCircleMarkers(lng = ~centCoorsGBIFsp[,1], lat = ~centCoorsGBIFsp[,2],color=c("blue"),label = ~(c("1970s","Both","2010s")), #) %>%
  labelOptions = labelOptions(noHide = TRUE, direction = "bottom",offset=c(0,12), textsize = "18px", textOnly = FALSE)) %>%
addArrowhead(lng= ~ centCoorsGBIFsp[,1], lat= ~ centCoorsGBIFsp[,2],color="blue") %>% addGraticule(interval=2.5)

