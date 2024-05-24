###____________________________________________________________________________________________________###
## Phylogenetic analysis - measuring, mapping, and anlyzing the effects of phylogenetic autocorrelation ##
###____________________________________________________________________________________________________###

# Load packages & data
# Data loading and preparation
# Fig S8 - Visualize extirpation on the phylogeny
# Fig S7 - Create and plot phylogenetic correlogram
# Phylogenetic eigenvector maps (see Caetano et al 2023 for full description of process)
# Accuracy of PEM alone (w/out traits) as predictors of extirpation in neural networks - Line 145
# Effects of PEM on all traits' Olden values from neural networks - Line 202
# Repeating analysis on dichotomous (categorical) bee subsets - Line 300

## Load necessary packages & data -------------------------------------------------------
library(dplyr) # for data wrangling
library(ggplot2) # basic plotting
library(readxl) # loading in excel data
library(stringr) # manipulating strings

library(picante)
library(ape) # basic management of phylogenetic data
library(phylobase) # for attaching trait data to phylogenies
library(phylosignal) # for measuring phylogenetic distance, Moran's I, phylogenetic correlograms (all to determine autocorrelation here)
library(phytools) # did we en up actually using this?
library(ggtree) # for plotting/visualizing phylogenies
library(MPSEM) # for creating PEMs based on traits in our phylogeny

library(nnet) # for creating, training, and using neural networks
library(NeuralNetTools) # for analyzing networks w/ Olden analysis

## Load our species traits file:
traits=read_excel("C:/Users/prglaum/Documents/Pollinators-Empirical/CCmanuscript/pubAttempt2/githubCode/Grahametal2023_BeeTraits.xlsx",sheet="TraitData",na = "NA")

## Load in nwk file of Max. Likelihood tree from http://beetreeoflife.org/ based on our Evans' species list:
beeTree=read.tree("C:/Users/prglaum/Documents/Pollinators-Empirical/CCmanuscript/pubAttempt2/Revisions_R1/ML_nwk_version.txt")
## See http://beetreeoflife.org/ for more on creating species-specific prunedphylogenies from their global phylogeny.


### Data preparation ------------------------------------------------------------------------------------
### prep our phylogeny to account for JUST the 72/73 species....

##Aligning our species list with what came from bee tree of life...
##make naming format the same (use "_" instead of " " for species names)
traits$species2=str_replace_all(traits$species, c(" " = "_" ));

## Create trait data subset of only species that appeared in 1972/1973 AND that exist in our phylogeny...
matchedData=subset(traits,bin>0 & bin!=2 & species2 %in% beeTree$tip.label)
## This covers 127/135 ~ 94% of our species, so we have good coverage here.

## Use keep.tip to remove tree structure not pertaining to our species of interest
matchedTree=keep.tip(beeTree,matchedData$species2 )
length(matchedTree$tip.label) ## Now our tree depicts the 127 species we will analyze

## Align our matchedData trait data subset structure with the matchedTree tip label order...
matchedData=matchedData[order(match(matchedData$species2,matchedTree$tip.label)),];
## Use the species2 label (names w/ "_") to create row labels in the trait data ("matchedData")
## Row labels are used to connect traits to tip labels for phylogenetic trait analysis.
rownames(matchedData)=matchedData$species2
matchedData$species2 == matchedTree$tip.label ##We can check here to make the all align correctly.


##  Visualize extirpation on the phylogeny - Fig S8-------------------------------------------------------
## Note this can be done w/ any other trait. Just change trait name below and
## color pallete if using a continuous trait.

# Create list of species extirpate and those that persisted
exts<-list(extirpated=matchedData$species2[matchedData$extinct==1],
           persisted=matchedData$species2[matchedData$extinct==0])

# now create a grouped tree based on selected OTUs organized by our exts trait list above.
matchedTreeNew <- groupOTU(matchedTree,.node = exts)

## Create Figure S8
# use ggTree to visualize extinction persistence on phylogeny
ggtree(matchedTreeNew, layout='circular') + geom_tree(aes(color=group)) +
  scale_color_manual(values=c("red","blue")) +
  geom_tiplab(size=2.4,aes(color=group))

## Create and plot phylogenetic correlogram - Figure S7 -------------------------------------------
## Use phylobase & phylosignal packages to create phylogenetic correlogram of extinction/persist

# Use phylobase's phylo4d function to create an R object where traits are mapped to the phylogeny.
BaseTreeStats=phylo4d(matchedTree, matchedData) #matchedData) ##via phylobase
#save(BaseTreeStats,file="extCorrelInput.RData") # object connecting traits to phylogeny can be saved

# Using the phylosignal package, create a correlogram of a chosen trait.
# Here we use extirpation/persistence as a "trait" to create Fig S7:
phyCorr=phyloCorrelogram(BaseTreeStats, trait="extinct") ##via phylosignal
plot(phyCorr)
# Check correlation on some other traits, e.g., ITD, which can show significant correlation
phyCorrpm=phyloCorrelogram(BaseTreeStats, trait="ITD") ##via phylosignal
plot(phyCorrpm)

## It is also possible to obtain a general/average Moran's I metric per trait.
stats=phyloSignal(BaseTreeStats[,"extinct"])
stats$stat$I
data.frame(matrix(unlist(stats), nrow=5, byrow=T))


##_______________________________##
## PHYLOGENETIC EIGENVECTOR MAPS ## ----------------------------------------------------------
##_______________________________##
## Here we use the MPSEM package and follow processes outlined in Caetano et al 2023 to measure
## phylogenetic structure w/ PEMs and account for that structure in training our neural networks

library(MPSEM) #load if not already loaded

## making sure all tip labels and species names match correctly
all(matchedTree$tip.label== matchedData$species2)

## First use the phylogeny as a basis to create a directed graph
bee.pgraph <- Phylo2DirectedGraph(matchedTree)

# Use our directed graph and a representative trait to create the phylogenetic
# eigenvector maps (PEM). Note, we use ITD here as it's the most phylogenetically correlated trait.
# PEM contains both eigenvector and corresponding eigenvalue information per vector.
PEM <-
  MPSEM::PEM.fitSimple(
    y = matchedData$ITD, #matchedData$ITD,
    x = NULL,
    w = bee.pgraph,
    d = "distance",
    sp = "species",
    lower = 0.1,
    upper = 1)

## Turn our eigenvectors into a dataframe so we can 1) use them as neural network inputs & 2) combine
## them with our trait variables into a dateframe to use as neural network inputs.
PEM_df <- as.data.frame(PEM)
names(PEM_df) <- paste0("PEM", 1:ncol(PEM_df))

# To avoid a unwieldly number of inputs, we follow protocol from Caetano et al 2023 and
# select only broader scale vectors (higher eigenvalues) at the 0.1 level.
PEM_sel <- PEM_df[,PEM$d > max(PEM$d)*0.1]
# Reorder to match data table
# PEM_sel <- PEM_sel[match(skink_data$tonini_name, matchedTree$tip.label),]

# Save it so you won't have to do it again
#saveRDS(PEM_sel, "PEM_Evans.rds")

# Load neural network packages if not already loaded
library(nnet)
library(NeuralNetTools)

## Accuracy of PEM alone as predictors of extirpation ---------------------------------------------------
## First, we investigate the performance of PEMs as input features in neural networks
# Add out extirpation/persistence column to our selected eigenvectors.
PEM_sel$extinct=matchedData$extinct
# Accuracy tests can be run w/ all eigenvectors above the .1 cutoff or the first 11 (see below).
# We used to the first 11 eigenvectors to balance with the 11 initial trait variables.
# Using the full 50 [1:49] eigenvectors that made the .1 cutoff creates an overloaded input layer for a 6 node neural network.
# To avoid an overloaded input layer and a sparse feature space, we'd need to change our neural neural network structure and obtain more data. This change hinders our
# ability to directly test for phylogenetic structural effects on the exact neural networks used in our analysis.
# This would also interfere with the clarity on trait based inference via edge weights (Olden values). We therefore limited our usage of PEMs to the first 11 eigenvectors.
# However no qualitative changes were found all using 50 eigenvectors in accuracy and original Olden values still generally correlate with Olden values w/ PEM regardless.

## Add aligned "extinction" variable as
PEM_use=data.frame(extinct=PEM_sel$extinct,PEM_sel[,1:11])# [,1:11] for same eigenvectors as inference check & [,1:49] for all eigenvectors above cutoff
accuRec=rep();
for (i in 1:1000) {
  myNNFull=nnet(as.factor(extinct)~., PEM_use, maxit=1500, size=6,decay = 5e-3,trace=FALSE)

  test=table(PEM_use$extinct, predict(myNNFull, newdata = PEM_use[,-1],type = "class"))
  accuRec[i]=sum(diag(test))/sum(test)
}
mean(accuRec) ## Check mean training accuracy
olden(myNNFull)+coord_flip() ## plot olden analysis if interested...

## Prepare data
## Use only the species and species trait data from our species of interest...
traits_use=subset(traits,species2 %in% matchedData$species2)
## Order species and their data correctly to match our
traits_use=traits_use[order(match(traits_use$species2,matchedData$species2)),];
traits_use$species2 == rownames(PEM_use)

traits_use$oliC=as.numeric(traits_use$Nest_location=="cavity"&traits_use$Dietary_specialization=="oligolectic")
traits_use$oliG=as.numeric(traits_use$Nest_location=="ground"&traits_use$Dietary_specialization=="oligolectic")
traits_use$polyC=as.numeric(traits_use$Nest_location=="cavity"&traits_use$Dietary_specialization=="polylectic")
traits_use$polyG=as.numeric(traits_use$Nest_location=="ground"&traits_use$Dietary_specialization=="polylectic")
traits_use$clepto=as.numeric(traits_use$Sociality=="clepto-parasitic")
traits_use$clepto[is.na(traits_use$clepto)]=0

## select our specific traits
use=dplyr::select(traits_use,c(naLatMin,naLatMax,naLongMin,naLongMax,phenoMean,phenoRange,ITD) )
## scale quantitative traits
use=as.data.frame(scale(use));
##Add the dichotomous nest/diet variables:
use$oliC=as.factor(traits_use$oliC); use$oliG=as.factor(traits_use$oliG); use$polyC=as.factor(traits_use$polyC);
use$polyG=as.factor(traits_use$polyG); use$clepto=as.factor(traits_use$clepto);
use=cbind(extinct=traits_use$extinct,use);

## Training iterations on our initial neural network input feature set. Reproduced here only for ease of comparison
normAccuRec=rep();
for (i in 1:1000) {
  myNNFull=nnet(as.factor(extinct)~oliG+polyC+clepto+polyG+ITD+naLongMax+naLongMin+naLatMax+naLatMin+phenoRange+phenoMean, use, maxit=1500, size=6,decay = 5e-3,trace=FALSE)
  ##Test our NN performance by testing its predictions against our data.
  test=table(use$extinct, predict(myNNFull, newdata = use[,-1],type = "class"))
  normAccuRec[i]=sum(diag(test))/sum(test)
}
mean(normAccuRec)

####  Checking trait Olden inference from neural networks using traits & PEM -----------------------------------
## RUNNING OUR FULL 1000 iterations, to see if get back the same inference for our trait variables
## after taking phylo correlation into account
useALL=cbind(use,PEM_use[,-1])
allAccuRec=rep();

myNNFull2=nnet(as.factor(extinct)~., useALL, maxit=1500, size=6,decay = 5e-3,trace=FALSE)
test=table(useALL$extinct, predict(myNNFull2, newdata = useALL[,-1],type = "class"))
allAccuRec[1]=sum(diag(test))/sum(test)

idata6=olden(myNNFull2)$data; #pull data from Olden analysis
#Rescale data between -1 and 1 for relative importance that can be compared
#across unique trainings.
rescaledImp=idata6$importance/max(abs(idata6$importance));
idata6$RSimportance=rescaledImp; #store rescaled olden importance values (RSO)
idata6$rank=rank(idata6$importance); #store rank importance
idata6$absrank=rank(-abs(idata6$importance)); #store abs(importance)
idata6$iter=rep(1,length(23)); #label this training iteration
tdf=garson(myNNFull2)$data; #include the Garson measure of importance, called "rel_imp" in dataframe
idata6=merge(x = idata6, y = tdf, by = c("x_names"),all.x=TRUE)

for (i in 1:500) {
  myNNFull2=nnet(as.factor(extinct)~., useALL, maxit=1500, size=6,decay = 5e-3,trace=FALSE)

  tdata=olden(myNNFull2)$data;
  ##store data as described above
  tdata$RSimportance=tdata$importance/max(abs(tdata$importance));
  tdata$rank=rank(tdata$importance);
  tdata$absrank=rank(-abs(tdata$importance));
  tdata$iter=rep(i+1,61);#23);
  tdf=garson(myNNFull2)$data;
  tdata=merge(x = tdata, y = tdf, by = c("x_names"),all.x=TRUE)
  idata6=rbind(idata6,tdata);

  ##Test our NN performance by testing its predictions against our data.
  test=table(use$extinct, predict(myNNFull, newdata = useALL[,-1],type = "class"))
  allAccuRec[i+1]=sum(diag(test))/sum(test)

}
mean(allAccuRec)
olden(myNNFull2)+coord_flip()

## Save so you don't have to rerun the 1000 iterations...
#write.csv(idata6,"phyloAuto&FullTraits_Olden_1000iter.csv")
idata6=read.csv("~/phyloAuto&FullTraits_Olden_1000iter.csv") ## load saved Olden data if already completed and saved.

## Order and factorize Olden data.
idata6=idata6[order(idata6$RSimportance),]
idata6$x_names <- factor(idata6$x_names, levels = unique(idata6$x_names))

## Basic box and whisker plot of ALL the features' Olden values. Note if using all the eigenvectors, this
## is Not recommended as there are too many features to efficiently display.
ggplot(data=idata6,aes(x=factor(x_names),y=RSimportance)) +#,level = c("polyC","maxLong","minLat","phenoMean","polyG","phenoRange","minLong","maxLat","ITD","clepto","oliG")),y=RSimportance)) +
  coord_flip()+ theme(text = element_text(size=15))+
  geom_boxplot(fill="gray") + ylab("Rescaled Olden Imp") + xlab("")+ theme_bw()

## For the purposes of creating te ridgeplot, create a subset of overall Olden data to just focus on the trait features
## since they are the focus of our analysis for comparison to neural networks trained w/out PEM.
idatatraits=subset(idata6,x_names %in% c("polyC1","naLongMax","phenoRange","naLatMin","phenoMean","naLatMax","polyG1","ITD","naLongMin","oliG1","clepto1") )

library(ggridges) #load necessary packages
library(viridis)
ggplot(idatatraits, aes(x = `RSimportance`, y = factor(x_names,level = c("polyC1","naLongMax","phenoRange","naLatMin","phenoMean","naLatMax","polyG1","ITD","naLongMin","oliG1","clepto1")) , fill = ..x..)) +
  geom_density_ridges_gradient(scale = 2.1, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Olden\nImportance", option = "C") +
  labs(title = 'Rescaled Importance Values w/ PEM',x="Importance",y="Traits") + theme_bw() +
  scale_y_discrete(labels=c("Poly.\nCavity","Max. Long","Pheno. Range","Min. Lat","Pheno. Mean","Max. Lat","Poly.\nGround","ITD","Min Long","Oligo.\nGround","Klepto."))


###Comparing the Olden values from NN training w/ and w/out PEM to see if phylo. autocorrelation inclusion affects the
## their original inference...
## Here again, load necessary data if not running full analysis and instead using existing save file provided in github.
idata6=read.csv("~/phyloAuto&FullTraits_Olden_1000iter.csv") ## load saved Olden data if already completed and saved.
idatatraits=subset(idata6,x_names %in% c("polyC1","naLongMax","phenoRange","naLatMin","phenoMean","naLatMax","polyG1","ITD","naLongMin","oliG1","clepto1") )

newRS=idatatraits %>%
  group_by(x_names) %>%
  summarize(avgRS=mean(RSimportance),sdRS=sd(RSimportance))

## Load original Olden data for traits from neural networks not using PEM:
idata6Original=read.csv("~/na_OldenEvansIsaacs_allbees_1000iter.csv")

OrigRS=idata6Original %>%
  group_by(x_names) %>%
  summarize(avgRS=mean(RSimportance),sdRS=sd(RSimportance))

#match dataframe order:
newRS=newRS[order(match(newRS$x_names,OrigRS$x_names)),];
scatterlabels=c("klepto","ITD","latMax","latMin","longMax","longMin","oliG","pheno\nMean","pheno\nRange","polyC","polyG")

ggplot(data.frame(names=scatterlabels,orig=OrigRS$avgRS,phylo=newRS$avgRS), aes(x=orig,y=phylo,label=names) ) + theme_bw() +
  theme(text = element_text(size=15)) +   geom_text(vjust = 0, nudge_y = 0.05) +
  geom_point(size=2.5) + stat_smooth(method="lm") + xlab("Mean trait Olden values only using traits") + ylab("Mean trait Olden values incorporating PEM")

## Check explanatory power between the two Olden datasets:
summary(lm(newRS$avgRS~OrigRS$avgRS) )

####_________SUBSETS ____________________________________
####SUBSETS!!!!!!!!!!!!!!!! --------SUBSETS-------------------------------------
####_________SUBSETS ____________________________________
## RUNNING OUR FULL 1000 iterations, to see if get back the same inference for our trait variables
## after taking phylo correlation into account
useALLsub=subset(useALL,polyG==1); #cbind(use,PEM_use[,-1])
useALLsub=select(useALLsub,-c(polyG,polyC,oliG,oliC,clepto))
allAccuRec=rep();

myNNFull2=nnet(as.factor(extinct)~., useALLsub, maxit=1500, size=6,decay = 5e-3,trace=FALSE)
test=table(useALLsub$extinct, predict(myNNFull2, newdata = useALLsub[,-1],type = "class"))
allAccuRec[1]=sum(diag(test))/sum(test)

idata6=olden(myNNFull2)$data; #pull data from Olden analysis
#Rescale data between -1 and 1 for relative importance that can be compared
#across unique trainings.
rescaledImp=idata6$importance/max(abs(idata6$importance));
idata6$RSimportance=rescaledImp; #store rescaled olden importance values (RSO)
idata6$rank=rank(idata6$importance); #store rank importance
idata6$absrank=rank(-abs(idata6$importance)); #store abs(importance)
idata6$iter=rep(1,length(18)); #label this training iteration
tdf=garson(myNNFull2)$data; #include the Garson measure of importance, called "rel_imp" in dataframe
idata6=merge(x = idata6, y = tdf, by = c("x_names"),all.x=TRUE)

for (i in 1:1000) {
  myNNFull2=nnet(as.factor(extinct)~., useALLsub, maxit=1500, size=6,decay = 5e-3,trace=FALSE)

  tdata=olden(myNNFull2)$data;
  ##store data as described above
  tdata$RSimportance=tdata$importance/max(abs(tdata$importance));
  tdata$rank=rank(tdata$importance);
  tdata$absrank=rank(-abs(tdata$importance));
  tdata$iter=rep(i+1,18);
  tdf=garson(myNNFull2)$data;
  tdata=merge(x = tdata, y = tdf, by = c("x_names"),all.x=TRUE)
  idata6=rbind(idata6,tdata);

  ##Test our NN performance by testing its predictions against our data.
  test=table(useALLsub$extinct, predict(myNNFull2, newdata = useALLsub[,-1],type = "class"))
  allAccuRec[i+1]=sum(diag(test))/sum(test)

}
mean(allAccuRec)
olden(myNNFull2)+coord_flip()

#write.csv(idata6,"R1_phyloAuto&polyG_Olden_1000iter.csv")

idata6=idata6[order(idata6$RSimportance),]
idata6$x_names <- factor(idata6$x_names, levels = unique(idata6$x_names))

ggplot(data=idata6,aes(x=factor(x_names),y=RSimportance)) +#,level = c("polyC","maxLong","minLat","phenoMean","polyG","phenoRange","minLong","maxLat","ITD","clepto","oliG")),y=RSimportance)) +
  coord_flip()+ theme(text = element_text(size=15))+
  geom_boxplot(fill="gray") + ylab("Rescaled Olden Imp") + xlab("")+ theme_bw()

idatatraits=subset(idata6,x_names %in% c("naLongMax","phenoRange","naLatMin","phenoMean","naLatMax","ITD","naLongMin") )
idatatraits500=subset(idata6,x_names %in% c("polyC1","naLongMax","phenoRange","naLatMin","phenoMean","naLatMax","polyG1","ITD","naLongMin","oliG1","clepto1") )

library(ggridges) #load necessary packages
library(viridis)
ggplot(idatatraits500, aes(x = `RSimportance`, y = factor(x_names,level = c("polyC1","naLongMax","phenoRange","naLatMin","phenoMean","naLatMax","polyG1","ITD","naLongMin","oliG1","clepto1")) , fill = ..x..)) +
  geom_density_ridges_gradient(scale = 2.1, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Olden\nImportance", option = "C") +
  labs(title = 'Rescaled Importance Values w/ PEM',x="Importance",y="Traits") + theme_bw() +
  scale_y_discrete(labels=c("Poly.\nCavity","Max. Long","Pheno. Range","Min. Lat","Pheno. Mean","Max. Lat","Poly.\nGround","ITD","Min Long","Oligo.\nGround","Clepto."))


###Comparing the Olden values from NN training w/ and w/out PEM to see if phylo. autocorrelation inclusion affects the
## their inference...

idata6=read.csv("C:/Users/prglaum/Documents/Pollinators-Empirical/CCmanuscript/pubAttempt2/Revisions_R1/R1_phyloAuto&clepto_Olden_1000iter.csv")
idatatraits=subset(idata6,x_names %in% c("naLongMax","phenoRange","naLatMin","phenoMean","naLatMax","ITD","naLongMin") )

newRS=idatatraits %>%
  group_by(x_names) %>%
  summarize(avgRS=mean(RSimportance),sdRS=sd(RSimportance))

#idata6Original=read.csv("C:/Users/prglaum/Documents/Pollinators-Empirical/CCmanuscript/pubAttempt2/na_OldenEvansIsaacs_allbees_1000iter.csv")
idataOriginalSUB=read.csv("C:/Users/prglaum/Documents/Pollinators-Empirical/CCmanuscript/pubAttempt2/na_clepto_OldenImp_1000rep.csv")

OrigRS=idataOriginalSUB %>%
  group_by(x_names) %>%
  summarize(avgRS=mean(RSimportance),sdRS=sd(RSimportance))

newRS=newRS[order(match(newRS$x_names,OrigRS$x_names)),];
scatterlabelsPOLYC=c("ITD","latMax","latMin","longMax","longMin","pheno\nMean","pheno\nRange")

ggplot(data.frame(names=scatterlabelsPOLYC,orig=OrigRS$avgRS,phylo=newRS$avgRS), aes(x=orig,y=phylo,label=names) ) + theme_bw() +
  theme(text = element_text(size=15)) +   geom_text(vjust = 0, nudge_y = 0.05) + ggtitle("d) Kleptoparasites") +
  geom_point(size=2.5) + stat_smooth(method="lm") + xlab("Mean trait Olden values only using traits") + ylab("Mean trait Olden values incorporating PEM")

summary(lm(newRS$avgRS~OrigRS$avgRS) )
