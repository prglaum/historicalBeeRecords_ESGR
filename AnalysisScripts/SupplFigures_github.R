###_______________________________________________________###
## Supplementary Results figures ##
###_______________________________________________________###

##All necessary data is located in the /data file at https://github.com/prglaum/historicalBeeRecords_ESGR/data &
##  https://agdatacommons.nal.usda.gov/account/articles/25233991
## load necessary packages
library(ggplot2) ##for plotting
library(cowplot) ##for arranging plots in grid
library(dplyr) ##data wrangling
library(FSA) ##Dunn test
library(lme4) ##for LMM
library(mgcv) ##for GAM
library(FSA) ##for Dunn test
library(readxl) ##for loading xlsx files
library(vegan) ##broadly applied package for ecological analysis. Used here for rarefaction
library(mgcv) ##for regressions via generalized additive models (GAM)

### -------------------Supplementary Figures - Contents ---------------------------
# Load Data -
##_________Supplementary Methods Figures________________________
# Fig S1 - Pictorial diagram created in powerpoint
# Fig S2 - Code for Fig S2 is available in geoplotting.R script
# Fig S3 - Code available upon request
# Fig S4 - Simply use the "plot" function on a nnet neural network object, e.g., plot(nnetObject)
# Fig S5 - Code starting on line 49 (for data prep) and line 88 (for plotting)
# Fig S6 - Code starting on line 49 (for data prep) and line 136 (for plotting)
# Fig S7 - Code available in phylogeneticAutocorrelation.R script
# Fig S8 - Code available in phylogeneticAutocorrelation.R script
# Fig S9 - Pictorial diagram created in powerpoint
##_________Supplementary Results Figures________________________
# Fig S10 - Figure made in GraphPad Prism. Contact Kelsey Graham w/ questions.
# Fig S11 - Code for rarefaction starting and for GAM plot starting on line 175
# Fig S12 - Figure made during document revisions. Contact Kelsey Graham w/ questions.
# Fig S13 - Code available in phylogeneticAutocorrelation.R script
# Fig S14 - Code available in commonBees.R script
# Fig S15 - S19 - Neural network figures can be recreated using code in NNscript_analysis_plots.R script
# Fig S20 - Code available in phylogeneticAutocorrelation.R script
# Fig S21 - Code available in commonBees.R script
# Fig S22 - Load data on line 49 and code starting on line 239
# Fig S23 - Load data on line 49 and data prep on line 270 Plot code begins on line 198
# Fig S24 - Load data on line 49 and code starting on line 293
# Fig S25 - Code starting on line 325
# Fig S26 - Code for Fig S18 is available in geoPlotting.R script
# Fig S27 - Figure made in GraphPad Prism. Contact Kelsey Graham w/ questions.
# Fig S28 - Code starting on line 364
# Fig S29 - Figure made in GraphPad Prism. Contact Kelsey Graham w/ questions.

### LOAD NECESSARY DATA --------------------------------------------------------

###Load necessary trait data:
traits=read_excel("~/Grahametal2023_BeeTraits.xlsx",sheet="TraitData",na = "NA")
##We use the "bin" variable to indicate species presence in across our two
##major sample periods (72/73 & 2017/2018).
## Bin = 1: Species only found in 72/73
## Bin = 2: Species only found in 2017/2018
## Bin = 3: Species found in both 72/73 & 2017/2018

###limit data for NN extinction analysis to species that were present in 72/73 and/or 2018/2019
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
##This is used to run our neural network to predict species persistance based on traits.
##We do this with the bin variable, set it either 1 or 3 to check which species
##found in in the 70s made it to the 2010s.
survive=subset(wtraits,bin!=2)

##Load sample data from Evans/Isaacs sample periods w/ trait data attached to each specimen record:
EvansIsaacs=read.csv("C:/Users/prglaum/Documents/Pollinators-Empirical/CCmanuscript/pubAttempt2/githubCode/ESGR_samples_EvansIsaacs_withTraitsMatch.csv",sep=",",header=T, na.strings=c("NaN","#NAME?","-Inf",""))

##We use the vegan package to rarefy our richness data. To use this package,
##the data needs to be in a specific format, which we have made in advance and load below.
esgrV=read.csv("~/AllVettedVeganESGR_Sept14_2022.csv",sep=",",header=T, na.strings=c("NaN","#NAME?","-Inf",""))

##______________________________________________________________
##_________Supplementary Methods Figures________________________
##______________________________________________________________

####Fig S5: DIET, plot differences across dietary range: -----------------------
##Note, letter labels were determined via Dunn Test above (alter code directly below
##to run any specific test of interest).
##Test differences across sample periods for any categorical trait in our data set.
kruskal.test(pheno_q90~Dietary_specialization,data=survive)
#library(FSA)
DT=dunnTest(naLongMax~Dietary_specialization,data=survive,method="bh"); DT
DT=dunnTest(phenoRange~Dietary_specialization,data=survive,method="bh"); DT

##Plot across traits at species level across dietary specialization trait
plot_grid(
ggplot(data=subset(survive,!is.na(phenoRange)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=phenoRange) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("Pheno. Range (days)") + ylim(25,195) + scale_x_discrete(labels=c("klepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","b","c"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(phenoMean)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=phenoMean) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("Pheno. Mean (days)") + ylim(50,330) + scale_x_discrete(labels=c("klepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","ab","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(ITD)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=ITD) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("ITD") + ylim(0.5,8) + scale_x_discrete(labels=c("klepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(na_hull_area)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=log(na_hull_area)) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("log(NA Hull Area)") + ylim(11.5,20) + scale_x_discrete(labels=c("klepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","a","b"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLatMax)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=naLatMax) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("NA Max Lat") + ylim(40,80) + scale_x_discrete(labels=c("klepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","a","b"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLatMin)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=naLatMin) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("NA Min Lat") + ylim(8,55) + scale_x_discrete(labels=c("klepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","a","ab"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLongMax)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=naLongMax) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("NA Max Long") + ylim(-85,-25) + scale_x_discrete(labels=c("klepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","b","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLongMin)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=naLongMin) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("NA Min Long") + ylim(-165,-55) + scale_x_discrete(labels=c("klepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","a","b"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15)),
ncol=4,nrow=2)

##Test group differences for Nest location.
kruskal.test(pheno_q90~Nest_location,data=survive)
DT=dunnTest(naLongMin~Nest_location,data=survive,method="bh"); DT

####Fig S6: NEST, plot differences across nesting strategy: --------------------------
plot_grid(
ggplot(data=subset(survive,!is.na(phenoRange)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=phenoRange) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("Pheno. Range (days)") + ylim(25,195) + scale_x_discrete(labels=c("cavity","klepto","ground") ) +
  stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(phenoMean)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=phenoMean) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("Pheno. Mean (days)") + ylim(50,330) + scale_x_discrete(labels=c("cavity","klepto","ground") ) +
  stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(ITD)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=ITD) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("ITD") + ylim(0.5,8) + scale_x_discrete(labels=c("cavity","klepto","ground") ) +
  stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(na_hull_area)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=log(na_hull_area)) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("log(NA Hull Area)") + ylim(11.5,20) + scale_x_discrete(labels=c("cavity","klepto","ground") ) +
  stat_summary(geom="text", label=c("a","b","c"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLatMax)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=naLatMax) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("NA Max Lat") + ylim(40,80) + scale_x_discrete(labels=c("cavity","klepto","ground") ) +
  stat_summary(geom="text", label=c("a","b","b"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLatMin)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=naLatMin) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("NA Min Lat") + ylim(8,55) + scale_x_discrete(labels=c("cavity","klepto","ground") ) +
  stat_summary(geom="text", label=c("a","b","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLongMax)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=naLongMax) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("NA Max Long") + ylim(-85,-25) + scale_x_discrete(labels=c("cavity","klepto","ground") ) +
  stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLongMin)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=naLongMin) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("NA Min Long") + ylim(-165,-55) + scale_x_discrete(labels=c("cavity","klepto","ground") ) +
  stat_summary(geom="text", label=c("a","b","c"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15)),
ncol=4,nrow=2)

##______________________________________________________________
##_________Supplementary Results Figures________________________
##______________________________________________________________

####Fig S11: Rarefying richness data and plotting via Generalized Additive Models (GAM): --------------------------
##We use the vegan package to rarefy our richness data. To use this package,
##the data needs to be in a specific format, which we have made in advance and load below.
esgrV=read.csv("~/AllVettedVeganESGR_Sept14_2022.csv",sep=",",header=T, na.strings=c("NaN","#NAME?","-Inf",""))

##We use separate cut offs for inclusion into our regressions analysis
#Raw richness
spAbund <- rowSums(esgrV[,2:length(esgrV[1,])])
#Only rarefying years w/ at least 20 species
rarefied=rarefy(esgrV[,2:length(esgrV[1,])],20)
r20=rarefied
#Only rarefying years w/ at least 30 species
rarefied=rarefy(esgrV[,2:length(esgrV[1,])],30)
r30=rarefied
#Only rarefying years w/ at least 40 species
rarefied=rarefy(esgrV[,2:length(esgrV[1,])],40)
r40=rarefied
#Only rarefying years w/ at least 50 species
rarefied=rarefy(esgrV[,2:length(esgrV[1,])],50)
r50=rarefied
#Only rarefying years w/ at least 60 species
rarefied=rarefy(esgrV[,2:length(esgrV[1,])],60)
r60=rarefied
#Only rarefying years w/ at least 70 species
rarefied=rarefy(esgrV[,2:length(esgrV[1,])],70)
r70=rarefied

##Gather them in a dataframe for GAM regression analysis
raredYears=data.frame(Year=esgrV$Year, spAbund, r20, r30, r40, r50, r60, r70)

##	##	##	##
###GAM REGRESSION### library(mgcv) #for gams
##	##	##	##
#To check other cutoffs, change r50 to desired cutoff level.
#Note, the %%1 operator removes any entries turned to one because they were
#below the cutoff level.
GAMmod=gam(r50~s(Year),data=subset(raredYears,(r50%%1)>0) )
summary(GAMmod)

##Predict fit to data, create upper and lower bounds, store them in dataframe.
df_predict2=cbind(subset(raredYears,(r50%%1)>0),predict.gam(GAMmod,subset(raredYears,(r50%%1)>0),type="response",se.fit=TRUE))
df_predict2 <- within(df_predict2, {RR <- fit
LL2 <- fit - 1.96 * se.fit
UL2 <- fit + 1.96 * se.fit
})

#Plot Figure 1
ggplot(data=subset(raredYears,(r50%%1)>0)) +
  theme_bw() + theme(text = element_text(size=20)) + theme(plot.margin = unit(c(.4,.7,.3,.5), "cm")) +
  geom_point( aes(x=esgrV.Year,y=r50) ) +
  geom_ribbon(data=df_predict2,aes(x=esgrV.Year,ymin = LL2, ymax = UL2),inherit.aes = FALSE,alpha = .15) +
  geom_line(data=df_predict2,aes(x=esgrV.Year,y=RR,col="#00FF00"),inherit.aes = FALSE,size = 1.5) +
  geom_point( aes(x=esgrV.Year,y=r50) ) +
  labs(x = "Year", y = "Rarefied Richness") +
  theme(legend.position = "none")


###__________________________________________###
###SOME TRAIT CHANGES FROM 70s to 2010s#########  ##Figures S22 - S24 ------------------
###__________________________________________###
##Load necessary data if not already loaded:
EvansIsaacs=read.csv("C:/Users/prglaum/Documents/Pollinators-Empirical/CCmanuscript/pubAttempt2/githubCode/ESGR_samples_EvansIsaacs_withTraitsMatch.csv",sep=",",header=T, na.strings=c("NaN","#NAME?","-Inf",""))

###__________________________________________###
###PHENOLOGY CHANGES FRMO 70s to 2010s##########  ##Fig S22 -----------------------------
###__________________________________________###
plot_grid(
ggplot(EvansIsaacs,aes(x=Collector,y=phenoRange,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Range (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("a) All samples"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="polylectic"),aes(x=Collector,y=phenoRange,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Range (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("b) Polylectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="oligolectic"),aes(x=Collector,y=phenoRange,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Range (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("c) Oligolectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="clepto-parasitic"),aes(x=Collector,y=phenoRange,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Range (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("d) Klepto."),

ggplot(EvansIsaacs,aes(x=Collector,y=phenoMean,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Mean (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("e) All samples"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="polylectic"),aes(x=Collector,y=phenoMean,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Mean (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("f) Polylectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="oligolectic"),aes(x=Collector,y=phenoMean,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Mean (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("g) Oligolectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="clepto-parasitic"),aes(x=Collector,y=phenoMean,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Mean (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("h) Klepto."),
ncol=4,nrow=2)

kruskal.test(phenoRange~Collector,data=EvansIsaacs)
diets=subset(EvansIsaacs,Dietary_specialization=="polylectic") #clepto-parasitic
kruskal.test(phenoRange~Collector,data=diets)

kruskal.test(phenoMean~Collector,data=EvansIsaacs)
diets=subset(EvansIsaacs,Dietary_specialization=="polylectic") #clepto-parasitic
kruskal.test(phenoMean~Collector,data=diets)

###__________________________________________###
###ITD CHANGES FRMO 70s to 2010s: ##Fig S23 ----------------
###__________________________________________###

justGround=subset(EvansIsaacs,Nest_location=="ground")
justGroundnoBB=subset(EvansIsaacs,Nest_location=="ground"&Genus!="Bombus")
justClepto=subset(EvansIsaacs,Nest_location=="clepto-parasitic")

#Test used to assess difference between sample periods.
kruskal.test(ITD~Collector,data=justGround)
t.test(justGround$phenoRange[diets$Collector=="Evans"],justGround$phenoRange[diets$Collector=="Isaacs"])

allGround=ggplot(justGround,aes(x=Collector,y=ITD,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("ITD(mm)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("Ground nesting")
allGroundnoBB=ggplot(justGroundnoBB,aes(x=Collector,y=ITD,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("ITD(mm)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle(expression("Ground w/ no "~italic(Bombus)~""))
cleptoITD=ggplot(justClepto,aes(x=Collector,y=ITD,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("ITD(mm)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("Kleptoparasitic")

#create_plot()
plot_grid(allGround,allGroundnoBB,cleptoITD,ncol=3)


###__________________________________________###
###LATITUDE CHANGES FROM 70s to 2010s########### ##Figure S24 ---------------
###__________________________________________###
plot_grid(
ggplot(EvansIsaacs,aes(x=Collector,y=naLatMax,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Max") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("a) All samples"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="polylectic"),aes(x=Collector,y=naLatMax,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Max") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("b) Polylectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="oligolectic"),aes(x=Collector,y=naLatMax,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Max") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("c) Oligolectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="clepto-parasitic"),aes(x=Collector,y=naLatMax,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Max") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("d) Klepto."),

ggplot(EvansIsaacs,aes(x=Collector,y=naLatMin,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Min") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("e) All samples"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="polylectic"),aes(x=Collector,y=naLatMin,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Min") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("f) Polylectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="oligolectic"),aes(x=Collector,y=naLatMin,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Min") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("g) Oligolectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="clepto-parasitic"),aes(x=Collector,y=naLatMin,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Min") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("h) Klepto"),
ncol=4,nrow=2)

kruskal.test(naLatMax~Collector,data=EvansIsaacs)
diets=subset(EvansIsaacs,Dietary_specialization=="clepto-parasitic") #clepto-parasitic
kruskal.test(naLatMax~Collector,data=diets)

kruskal.test(naLatMin~Collector,data=EvansIsaacs)
diets=subset(EvansIsaacs,Dietary_specialization=="polylectic") #clepto-parasitic
kruskal.test(naLatMin~Collector,data=diets)


###------------------###
##SPECIES LEVEL CHANGE## Fig S25 ----------------------
#______________________#

## Individual phenology plots
ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=phenoRange,group=bin)) + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") ) + theme_bw()
ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=phenoMean,group=bin)) + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") ) + theme_bw()
ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=pheno_q10,group=bin)) + geom_boxplot(fill="lightblue") + scale_x_discrete(labels=c("70s","Both","2010s") ) + theme_bw()
ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=pheno_q90,group=bin)) + geom_boxplot(fill="lightblue") + scale_x_discrete(labels=c("70s","Both","2010s") ) + theme_bw()

##Fig S25
##Phenological shift showing earlier flying SPECIES
plot_grid(
  ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=pheno_q10,group=bin)) + theme_bw() + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") ) + ylim(80,275) +
    geom_boxplot(fill="light blue") + ylab("Phenology 10th Quantile (days of year)") + xlab("") + stat_summary(geom="text", label=c("a","b","a"),fun=max, vjust=-1, size=10) + theme(text = element_text(size=15)) + ggtitle("a)"),
  ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=pheno_q90,group=bin)) + theme_bw() + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") ) + ylim(140,310) +
    geom_boxplot(fill="light blue") + ylab("Phenology 90th Quantile (day of year)") + xlab("") + stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text = element_text(size=15)) + ggtitle("b)"),
  ncol=2)

library(FSA)
kruskal.test(pheno_q90~bin,data=subset(traits,bin>0))
DT=dunnTest(pheno_q90 ~ as.factor(bin), data=subset(traits,bin>0), method="bh"); DT


##Not in Supplementary Results, but included here for completeness
##Longitudinal shifts
# naLongMax, naLongMin, naLong_q10, naLong_q90, naCentLong
plot_grid(
  ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=naLongMax,group=bin)) + theme_bw() + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") ) + ylim(-90,-40) +
    geom_boxplot(fill="light blue") + ylab("NA Long. Max") + xlab("") + stat_summary(geom="text", label=c("a","a","ab"),fun=max, vjust=-1, size=10) + theme(text = element_text(size=15)) + ggtitle("a)"),
  ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=naLongMin,group=bin)) + theme_bw() + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") ) + ylim(-170,-70) +
    geom_boxplot(fill="light blue") + ylab("NA Long. Min") + xlab("") + stat_summary(geom="text", label=c("a","a","ab"),fun=max, vjust=-1, size=10) + theme(text = element_text(size=15)) + ggtitle("b)"),
  ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=naCentLong,group=bin)) + theme_bw() + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") ) + #ylim(140,310) +
    geom_boxplot(fill="light blue") + ylab("NA Centroidal Longitude ") + xlab("") + stat_summary(geom="text", label=c("a","a","b"),fun=max, vjust=-1, size=10) + theme(text = element_text(size=15)) + ggtitle("c)"),
  ncol=3)

kruskal.test(naLongMin~bin,data=subset(traits,bin>0))
DT=dunnTest(naLongMin ~ as.factor(bin), data=subset(traits,bin>0), method="bh"); DT

###_______________________###
## Graphical representation of phenologica range effects in polylectic bees Fig S28 ------------------
###_______________________###
library(sjPlot) #load package to use plot_model function
##Add the polylectic flag variable, "poly" to the survive data frame
survive$poly=rep("n",length(survive$MinLong))
survive$poly[survive$Dietary_specialization=="polylectic"]="y"

survive$DS=as.factor(survive$Dietary_specialization)

use=dplyr::select(survive,c(naLatMin,naLatMax,naLongMin,naLongMax,pheno_q10,pheno_q90,phenoMean,phenoRange,ITD) )
use=as.data.frame(scale(use));
use$oliC=as.factor(survive$oliC); use$oliG=as.factor(survive$oliG); use$polyC=as.factor(survive$polyC); use$polyG=as.factor(survive$polyG); use$clepto=as.factor(survive$clepto);
use=cbind(extinct=survive$extinct,use,poly=survive$poly,DS=as.factor(survive$Dietary_specialization));

###_______________________###
## Figure 20a - Generalized linear model
mmod=(glmer(as.factor(extinct)~(1+(phenoRange)|DS)+(phenoRange),data=use,family=binomial))
summary(lmod)
library(sjPlot) ## NOTE! this package masks plot_grid from cowplot. Only load when done with cowplot, or detach sjplot when wanting to use cowplot again.

plot_model(lmod,type = "pred", terms = c("phenoRange[all]", "poly")) + theme_bw() +
  ggtitle("")+ylab("Prob. of Local Extirpation") + xlab("Pheno. Range") + guides(color=guide_legend(title="Polylecty")) +
  theme(text = element_text(size = 14)) + theme(legend.position = c(0.2, 0.2))
#plot_model(gmod,type = "pred", terms = c("poly","phenoRange[-1,2]"))

###_______________________###
## Figure 20b - Linear mixed model
## create a factor version of the Dietary_specialization variable, DS for dietary strategy
survive$DS=as.factor(survive$Dietary_specialization)
## load lme4 for LMM
library(lme4)

##Make a scaled version of the variables to limit the effects of singular fits. Name this data frame "use"
use=dplyr::select(survive,c(naLatMin,naLatMax,naLongMin,naLongMax,pheno_q10,pheno_q90,phenoMean,phenoRange,ITD) )
use=as.data.frame(scale(use));
use$oliC=as.factor(survive$oliC); use$oliG=as.factor(survive$oliG); use$polyC=as.factor(survive$polyC); use$polyG=as.factor(survive$polyG); use$clepto=as.factor(survive$clepto);
##Run model, predict, plot
mmod=(glmer(as.factor(extinct)~(1+(phenoRange)|DS)+(phenoRange),data=use,family=binomial))
pred=data.frame(fit=predict(mmod,newdata=use,type="response"))
dfp=data.frame(pR=use$phenoRange,DS=use$DS,extP=pred$fit)
ggplot(dfp,aes(x=pR,y=extP,color=DS)) + geom_line(linewidth=1.5) + theme_bw() +
  xlab("Pheno. Range") + ylab("Prob. of Local Extirpation") +
  scale_color_discrete(name = "Diet") + theme(text = element_text(size = 14))

###_______________________###
## Figure 20C - Generalized additive model
## Again use the scaled version of the trait variables
## model, create confidence intervals, predict, store prediction, plot
use$DS[use$DS=="clepto-parasitic"]="klepto"
gammod=gam(extinct ~ s(phenoRange, DS, bs="fs",k=4), data=use,family=binomial("logit"))
summary(gammod)
fam=family(gammod); ilink <- fam$linkinv;
pred=as.data.frame(predict(gammod,newdata=use,type="link",se.fit=TRUE))

dfp=data.frame(pR=use$phenoRange,DS=use$DS,fit_resp  = ilink(pred$fit),
               right_upr = ilink(pred$fit + (2 * pred$se.fit)),right_lwr = ilink(pred$fit - (2 * pred$se.fit)))

ggplot(dfp,aes(x=pR,y=fit_resp,color=DS,fill=DS)) + geom_line(linewidth=1.5) + theme_bw() +
  geom_ribbon(aes(ymin = right_lwr, ymax = right_upr), alpha = 0.1) +
  guides(fill="none",color=guide_legend(title="Diet")) + theme(text = element_text(size = 14)) +
  xlab("Pheno. Range") + ylab("Prob. of Local Extirpation") + scale_color_discrete(name="Diet",labels=c("klepto.","oligoectic","polylectic"))


