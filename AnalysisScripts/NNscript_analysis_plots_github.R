###______________________________________________________________________###
##Neural Network Analysis: 70s to 2010s sample period extirpation analysis##------
###______________________________________________________________________###

##Code to run neural network analysis used to make Fig 3 and Fig S9-S13
##All necessary data is located in the /data file at https://github.com/prglaum/historicalBeeRecords_ESGR

##   Contents----
# Load packages - Line 17
# Load data - Line 42
# Data preparation - Line 52 
# Neural network analysis on the all the major sample period extirpation/persistence data (1972/1973 and 2017/2018 periods) - Line 95
# Neural network analysis on the dichotomous SUBSETS based on species traits (i.e., polylectic cavity nesting bees) - Line 185
# Recreating Supplementary neural network figures Fig S9-S13 - Line 258
# Partial Dependence Plot [Feature Interaction, Fig S10] - Line 321

##   Load packages ----------------------
##Used for local parallelization
library("future") 
library("future.callr")

##Packages to complete Neural Network analysis:
library(nnet) #Neural Network package
library(iml) #Model-agnostic analysis package, PDPs and ICE plots
library(NeuralNetTools) #For Garson & Olden analysis of feature effects

#plotting/data visualization packages:
library(ggplot2) ##for plotting
library(ggridges)
library(viridis)

#for data wrangling/management
library(dplyr)
library(readxl) ##for loading xlsx files 

library(glmm) #Generalized linear model (mixed?)
library(psych) #I don't remember why this is loaded...
library(mgcv) #For Generalized additive models
library(lme4) #Generalized linear mixed models


##   Load data -----------------------------
##Load trait data
#traits=read.csv("C:/Users/prglaum/Documents/Pollinators-Empirical/CCmanuscript/deathProject/gbifTraits-Combo3.csv")
traits=read_excel("~/Grahametal2023_BeeTraits.xlsx",sheet="TraitData",na = "NA") #input your correct directory 
##We use the "bin" variable to indicate species presence across our two
##major sample periods (72/73 & 2017/2018). 
## Bin = 1: Species only found in 72/73
## Bin = 2: Species only found in 2017/2018
## Bin = 3: Species found in both 72/73 & 2017/2018

##   Data Prep & simple Fisher test-----------------------------
###limit data for NN extinction analysis to species that were present in 72/73 and/or 2017/2018 
wtraits=subset(traits,bin>0)

#Create the dichotomous variables defining combinations of diet and nesting. 
#These are needed to run the neural network
wtraits$oliC=as.numeric(wtraits$Nest_location=="cavity"&wtraits$Dietary_specialization=="oligolectic")
wtraits$oliG=as.numeric(wtraits$Nest_location=="ground"&wtraits$Dietary_specialization=="oligolectic")
wtraits$polyC=as.numeric(wtraits$Nest_location=="cavity"&wtraits$Dietary_specialization=="polylectic")
wtraits$polyG=as.numeric(wtraits$Nest_location=="ground"&wtraits$Dietary_specialization=="polylectic")
wtraits$clepto=as.numeric(wtraits$Sociality=="clepto-parasitic")
wtraits$clepto[is.na(wtraits$clepto)]=0

##Creating dataframe focused solely on species from 72/73 sample period. 
##This is used to run our neural network to predict species extirpation/persistence based on traits. 
##We do this with the bin variable, set it either 1 or 3 to check which species 
##found in in the 70s made it to the 2010s. 
survive=subset(wtraits,bin!=2)

###FISHER EXACT TESTS
##Two sided:
fisher.test(table(survive$Dietary_specialization,survive$extinct))
fisher.test(table(survive$Nest_location,survive$extinct))
##One sided:
fisher.test(table(survive$clepto,survive$extinct),alternative="greater")
fisher.test(table(survive$oliG,survive$extinct),alternative="greater")
fisher.test(table(survive$polyC,survive$extinct),alternative="less")
fisher.test(table(survive$polyG,survive$extinct),alternative="greater") 

##The "use" dataframe creates SCALED variables for neural network analysis.
##Quantitative variables of different scales need to be scaled to best interface
##with neural networks
use=dplyr::select(survive,c(naLatMin,naLatMax,naLongMin,naLongMax,phenoMean,phenoRange,ITD) )
use=as.data.frame(scale(use));
##Add the dichotomous nest/diet variables:
use$oliC=as.factor(survive$oliC); use$oliG=as.factor(survive$oliG); use$polyC=as.factor(survive$polyC); use$polyG=as.factor(survive$polyG); use$clepto=as.factor(survive$clepto);
use=cbind(extinct=survive$extinct,use);

##Remove species from neural network consideration that had missing trait data
use=use[complete.cases(use), ]


###__________________###
##    NNET ON FULL DATA (Fig 3) ##-----------------------------------
###Neural networks on all bees in the 72/73 vs 2017/2018 comparison
###__________________###
##make sure necessary Neural Network packages are loaded. 
## Note, if you don't want to recreate this whole process, skip below to load our 
## NN output data and go straight to plotting. 
## Also note, recreating a new 1000 trained neural networks will result in similar
## but quantitatively different looking plots given the random initial weights. 

##Train a neural network on our full trait dataset for the 72/73 species.
##In this example we're using the nnet package to run a single hidden layer
##neural net with 5 neurons. We use 1500 as the max number of iterations and 
##a decay rate of 5e-3. 
myNNFull=nnet(as.factor(extinct)~oliG+polyC+clepto+polyG+ITD+naLongMax+naLongMin+naLatMax+naLatMin+phenoRange+phenoMean, use, maxit=1500, size=6,decay = 5e-3,trace=FALSE) 

##Test our NN performance by testing its predictions against our data.
test=table(use$extinct, predict(myNNFull, newdata = use[,-1],type = "class"))
sum(diag(test))/sum(test)

##Investigate feature importance of this individual NN instance. These require the NeuralNetTools package
garson(myNNFull)+coord_flip() ##Non-directional importance measurement
olden(myNNFull)+coord_flip() ##This is the technique we use to measure relationships

####
##Now begin out process taking results from 1000 unique training runs. 
####

##Start a dataframe to hold data from NN iterations. 
idata6=olden(myNNFull)$data; #pull data from Olden analysis
#Rescale data between -1 and 1 for relative importance that can be compared 
#across unique trainings. 
rescaledImp=idata6$importance/max(abs(idata6$importance)); 

idata6$RSimportance=rescaledImp; #store rescaled olden importance values (RSO)
idata6$rank=rank(idata6$importance); #store rank importance
idata6$absrank=rank(-abs(idata6$importance)); #store abs(importance) 
idata6$iter=rep(1,11); #label this training iteration
tdf=garson(myNNFull)$data; #include the Garson measure of importance, called "rel_imp" in dataframe
idata6=merge(x = idata6, y = tdf, by = c("x_names"),all.x=TRUE)

##Vector to track accuracy of each training iteration:
test=table(use$extinct, predict(myNNFull, newdata = use[,-1],type = "class"))
accuracy=rep(); accuracy[1]=sum(diag(test))/sum(test);

for(i in 1:1000) { #loop through unique training runs
	##train NN with new starting initial weights
	myNNFull=nnet(as.factor(extinct)~oliG+polyC+clepto+polyG+ITD+naLatMax+naLatMin+naLongMax+naLongMin+phenoRange+phenoMean, use, maxit=1500, size=6,decay = 5e-3,trace=FALSE)
	##get output data from Olden analysis
	tdata=olden(myNNFull)$data;
	##store data as described above
	tdata$RSimportance=tdata$importance/max(abs(tdata$importance));
	tdata$rank=rank(tdata$importance);
	tdata$absrank=rank(-abs(tdata$importance)); 
	tdata$iter=rep(i+1,11);
	tdf=garson(myNNFull)$data;
	tdata=merge(x = tdata, y = tdf, by = c("x_names"),all.x=TRUE)
	idata6=rbind(idata6,tdata);

	test=table(use$extinct, predict(myNNFull, newdata = use[,-1],type = "class"))
	accuracy[i+1]=sum(diag(test))/sum(test)
}  #end of iterations
mean(accuracy)
idata6=idata6[order(idata6$RSimportance),]
idata6$x_names <- factor(idata6$x_names, levels = unique(idata6$x_names))

##Save your own version of the output RSO
#  write.csv(idata6,"naFeatureEffects_allData_6nodesLocal.csv")

##Because each NN starts with random initial weights, to recreate the exact figure
##from our manuscript we need to load in the data from our exact iterations. 
##Data is available here: https://github.com/prglaum/historicalBeeRecords_ESGR/tree/main/data 
idata6=read.csv("~/naOldenEvansIsaacs_allbees_1000iter.csv")

##Plot the RSO importance values from our iterations using box and whisker plot
ggplot(data=idata6,aes(x=factor(x_names),y=RSimportance)) +#,level = c("polyC","maxLong","minLat","phenoMean","polyG","phenoRange","minLong","maxLat","ITD","clepto","oliG")),y=RSimportance)) +
  coord_flip()+ theme(text = element_text(size=15))+
  geom_boxplot(fill="gray") + ylab("Rescaled Olden Imp") + xlab("")+ theme_bw() 

##Plot RSO importance values with RIDGE PLOTS ###
##This creates Figure 3
library(ggridges) #load necessary packages
library(viridis)
ggplot(idata6old, aes(x = `RSimportance`, y = factor(x_names,level = c("polyC1","naLongMax","phenoRange","naLatMin","phenoMean","naLatMax","polyG1","ITD","naLongMin","oliG1","clepto1")) , fill = ..x..)) +
  geom_density_ridges_gradient(scale = 2.1, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Olden\nImportance", option = "C") +
  labs(title = 'Rescaled Importance Values',x="Importance",y="Traits") + theme_bw() +
  scale_y_discrete(labels=c("Poly.\nCavity","Max. Long","Pheno. Range","Min. Lat","Pheno. Mean","Max. Lat","Poly.\nGround","ITD","Min Long","Oligo.\nGround","Clepto."))


####_________________###################_________________###
###   NEURAL NET on SUBSET DATA ###---------------------------
## Neural networks on dichotomous data subsets (e.g., polyC)
####_________________###################_________________#
## Use the subD subset object below to choose one of the 4 dichotomous 
## nest/diet categorical variables used in our analysis: 
#polyC=polylectic cavity nesting bees
#polyG=polylectic ground nesting bees
#oliG=oligolectic ground nesting bees
#clepto=clepto-parasitic bees

subD=subset(use,polyG==1)

##Train neural network. We achieved generall good fits using 3 neurons in our hidden layer here.
myNNsub=nnet(as.factor(extinct)~ITD+naLatMax+naLatMin+naLongMax+naLongMin+phenoRange+phenoMean, subD, maxit=1500, size=3,decay = 5e-3,trace=FALSE)
##test for accuracy of neural network, are we good enoug predictive power to continue? 
test=table(subD$extinct, predict(myNNsub, newdata = subD[,-1],type = "class"))
sum(diag(test))/sum(test)

##View interpretive plots
garson(myNNsub) + coord_flip()
olden(myNNsub) + coord_flip()

####Run ITERATIONS
####Get importance info, make storage dataframe
idataSub=olden(myNNsub)$data
rescaledImp=idataSub$importance/max(abs(idataSub$importance))
idataSub$RSimportance=rescaledImp;
idataSub$rank=rank(idataSub$importance);
idataSub$absrank=rank(-abs(idataSub$importance)); 
idataSub$iter=rep(1,11);
tdf=garson(myNNFull)$data;
idataSub=merge(x = idataSub, y = tdf, by = c("x_names"),all.x=TRUE)

subAccuracy=rep();  subAccuracy=sum(diag(test))/sum(test);

for(i in 1:1000) {
	##train NN with new starting initial weights
	myNNsub=nnet(as.factor(extinct)~ITD+maxLat+minLat+maxLong+minLong+phenoRange+phenoMean, subD, maxit=1500, size=3,decay = 5e-3,trace=FALSE)
	##get output data from Olden analysis
	tdata=olden(myNNsub)$data;
	##store data as described above
	tdata$RSimportance=tdata$importance/max(abs(tdata$importance));
	tdata$rank=rank(tdata$importance);
	tdata$absrank=rank(-abs(tdata$importance)); 
	tdata$iter=rep(i+1,11);
	tdf=garson(myNNsub)$data;
	tdata=merge(x = tdata, y = tdf, by = c("x_names"),all.x=TRUE)

	test=table(subD$extinct, predict(myNNsub, newdata = subD[,-1],type = "class"))
	subAccuracy[i+1]=sum(diag(test))/sum(test)
	idataSub=rbind(idataSub,tdata);
}  #end
mean(subAccuracy)
idataSub=idataSub[order(idataSub$RSimportance),]
idataSub$x_names <- factor(idataSub$x_names, levels = unique(idataSub$x_names))

#Save RSO value output if you recreated your own version of our analysis. Make sure to name appropriately
#write(idataSub,"naFeatureEffects_polyC_3nodesLocal.csv")

##Because each NN starts with random initial weights, to recreate the exact figure
##from our manuscript we need to load in the data from our exact iterations. Here we how the polylectic cavity nesters:
idataSub=read.csv("~/na_polyC_OldenImp_1000rep.csv")

idataSub=idataSub[order(idataSub$RSimportance),]

##Plot the RSO importance values from our iterations using box and whisker plot
##The full SI figure plot code is below. 
ggplot(data=idataSub,aes(x=factor(x_names),y=RSimportance)) +
  coord_flip()+ theme(text = element_text(size=15))+
  geom_boxplot(fill="gray") + ylab("Rescaled Olden Imp") + xlab("")+ theme_bw() 


###____###____###___####____###
##    Recreating supplementary neural network figures (Fig S9 - S13) ------------------------------
##Because each NN starts with random initial weights, to recreate the exact figure
##from our manuscript we need to load in the data from our exact iterations.
##
##Choose which subset of dichotomous bee neural network RSO values to plot. Exact data from our training runs
##are found at https://github.com/prglaum/historicalBeeRecords_ESGR/tree/main/data
## na_polyC_OldenImp_1000rep.csv - polylectic cavity nesters
## na_polyG_OldenImp_1000rep.csv - polylectic ground nesters
## na_oliG_OldenImp_1000rep.csv - oligolectic ground nesters
## na_clepto_OldenImp_1000rep.csv - clepto-parasitic 

##Here we show results from oligolectic ground nesters
idataSub=read.csv("~/na_oliG_OldenImp_1000rep.csv")

## Step 1:
##The following preserve order across layers. Run which ever group you need to
polyClevels=c("phenoRange","phenoMean","naLatMin","naLongMin","naLongMax","naLatMax","ITD")
polyCcols=c("lightblue","lightblue","lightblue","gray","gray","gray","red")
polyCnames=c("Pheno. Range","Pheno. Mean","NA Lat Min","NA Long Min","NA Long Max","NA Lat Max","ITD")

polyGlevels=c("phenoRange","naLongMax","naLatMin","naLatMax","naLongMin","ITD","phenoMean")
polyGcols=c("lightblue","lightblue","gray","gray","gray","gray","red")
polyGnames=c("Pheno. Range","NA Long. Max","NA Lat. Min","NA Lat. Max","NA. Long Min","ITD","Pheno. Mean")

oliGlevels=c("naLongMax","ITD","naLatMin","phenoMean","naLongMin","naLatMax","phenoRange")
oliGcols=c("gray","gray","gray","gray","gray","red","red")
oliGnames=c("NA Long. Max","ITD","NA Lat. Min","Pheno. Mean","NA Long. Min","NA Lat. Max","Pheno. Range")

cleptolevels=c("ITD","naLatMin","phenoMean","naLongMax","naLatMax","phenoRange","naLongMin")
cleptocols=c("lightblue","gray","gray","gray","gray","gray","red")
cleptonames=c("ITD","NA Lat. Min","Pheno. Mean","NA Long. Max","NA Lat. Max","Pheno. Range","NA Long. Min")

## Step 2:
##first plot, make sure to update calls to levels, cols, and names to fit the subset being plotted
firstgraph=ggplot(data=idataSub,aes(x=factor(x_names,levels=oliGlevels),
  y=RSimportance))+ coord_flip()+ theme_bw() +
  theme(text = element_text(size=15)) + geom_boxplot(fill="gray") + 
  ylab("Rescaled Olden Imp") + xlab("Traits") 

## Step 3: extract layer data from firstgraph in order to specifically color the outlying
# points the appropriate color:
ld<-layer_data(firstgraph); #get graph layer data
OLS=ld$outliers; OLSlens=as.numeric(lapply(ld$outliers,length)); # get outliers & length of that vector to create cols length
x_names=rep(oliGlevels,OLSlens)
##change color depending on group plotted:
usecol=rep(c("gray","gray","gray","gray","gray","red","red"),OLSlens)
##new data frame with outliers and their appropriate colors
outs=data.frame(x_names=x_names,RSimportance=unlist(OLS),usecol=usecol )

## Step 4:
##Plot the whole thing. Make sure to do any appropriate updates if plotting different bee group
ggplot(data=idataSub,aes(x=factor(x_names,levels=oliGlevels ),
  y=RSimportance))+ coord_flip()+ theme_bw() +
  theme(text = element_text(size=15)) +  
  geom_hline(aes(yintercept=0.5),linetype = "longdash") + geom_hline(aes(yintercept=-0.5),linetype = "longdash") +
  geom_boxplot(fill=oliGcols,outlier.shape = NA ) +
  geom_point(data = outs,colour=usecol ) +
  scale_x_discrete(labels=oliGnames) +
  xlab("Traits") +  ylab("Rescaled Olden Imp\nPersistace<--------->Extirpation")



###____######____###
###   INTERACTIONS: Partial Dependence Plots (Fig S10) ##------------------------
## Note, the iml package must be loaded to run this: 
## Choose a dichotomous bee group (here we use polylectic cavity nesters)
subD=subset(use,polyC==1)
## Train neural network. We achieved generall good fits using 3 neurons in our hidden layer here.
myNNsub=nnet(as.factor(extinct)~ITD+naLatMax+naLatMin+naLongMax+naLongMin+phenoRange+phenoMean, subD, maxit=1500, size=3,decay = 5e-3,trace=FALSE)

## Create a predictor object with iml package.
predictorSub <- Predictor$new(myNNsub, data = subD[,-1], y = subD[,1])
## Check for general interactions...
interact  <- Interaction$new(predictorSub)
interact$results

## Make 2D partial dependence plot (PDP) for phenoMean and phenoRange: 
## Note, PDPs are made from single neural networks. Recreations made here will look similar
## but not identical to the version published in the Supplementary Results. 
pdp2F=FeatureEffect$new(predictorSub, feature = c("phenoMean","phenoRange"),method="pdp")
pdp2F$plot()+theme_bw()+ theme(text = element_text(size=14)) + scale_fill_continuous(name = "Predicted\nProb. of\nExtinction") +
  ylab("Pheno. Range") + xlab("Pheno. Mean") + labs(title = "PDP: Polylectic cavity nesting")



