###_______________________________________________________###
## Supplementary Results figures ##
###_______________________________________________________###

##All necessary data is located in the /data file at https://github.com/prglaum/historicalBeeRecords_ESGR

## load necessary packages
library(ggplot2) ##for plotting
library(cowplot) ##for arranging plots in grid
library(dplyr) ##data wrangling
library(FSA) ##Dunn test
library(lme4) ##for LMM
library(sjPlot) ##for plotting
library(mgcv) ##for GAM
library(readxl) ##for loading xlsx files

### Supplementary Figures - Contents ---------------------------
# Fig S1 - Pictorial diagram created in powerpoint
# Fig S2 - Code for Fig S2 is available in geoplotting.R script
# Fig S3 - Code available upon request
# Fig S4 - Simply use the "plot" function on a nnet neural network object, e.g., plot(nnetObject)
# Fig S5 - Code starting on line 36 (for data prep) and line 74 (for plotting)
# Fig S6 - Code starting on line 36 (for data prep) and line 113 (for plotting)
# Fig S7 - Pictorial diagram created in powerpoint
# Fig S8 - Figure made in GraphPad Prism. Contact Kelsey Graham w/ questions.
# Fig S9 - S13 - Neural network figures can be recreated using code in NNscript_analysis_plots.R script
# Fig S14 - Load data on line 150 and code starting on line 156
# Fig S15 - Load data on line 150 and data prep on line 187. Plot code begins on line 198
# Fig S16 - Load data on line 150 and code starting on line 210
# Fig S17 - Code starting on line 242
# Fig S18 - Code for Fig S18 is available in geoPlotting.R script
# Fig S19 - Figure made in GraphPad Prism. Contact Kelsey Graham w/ questions.
# Fig S20 - Code starting on line 255


###Quantitative trait differences across categorical traits (sp. level): Fig S5-Fig S6 -----------------

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

##Test differences across sample periods for any categorical trait in our data set.
kruskal.test(pheno_q90~Dietary_specialization,data=survive)
#library(FSA)
DT=dunnTest(naLongMax~Dietary_specialization,data=survive,method="bh"); DT
DT=dunnTest(phenoRange~Dietary_specialization,data=survive,method="bh"); DT

##Plot across traits at species level.
##Note, letter labels were determined via Dunn Test above (alter code directly above
##to run any specific test of interest).

####Fig S5: DIET, plot differences across dietary range: -----------------------
plot_grid(
ggplot(data=subset(survive,!is.na(phenoRange)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=phenoRange) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("Pheno. Range (days)") + ylim(25,195) + scale_x_discrete(labels=c("clepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","b","c"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(phenoMean)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=phenoMean) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("Pheno. Mean (days)") + ylim(50,330) + scale_x_discrete(labels=c("clepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","ab","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(ITD)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=ITD) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("ITD") + ylim(0.5,8) + scale_x_discrete(labels=c("clepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(na_hull_area)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=log(na_hull_area)) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("log(NA Hull Area)") + ylim(11.5,20) + scale_x_discrete(labels=c("clepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","a","b"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLatMax)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=naLatMax) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("NA Max Lat") + ylim(40,80) + scale_x_discrete(labels=c("clepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","a","b"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLatMin)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=naLatMin) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("NA Min Lat") + ylim(8,55) + scale_x_discrete(labels=c("clepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","a","ab"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLongMax)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=naLongMax) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("NA Max Long") + ylim(-85,-25) + scale_x_discrete(labels=c("clepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","b","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLongMin)&!is.na(Dietary_specialization)),aes(x=as.factor(Dietary_specialization),y=naLongMin) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Diet") + ylab("NA Min Long") + ylim(-165,-55) + scale_x_discrete(labels=c("clepto","oligo","poly") ) +
  stat_summary(geom="text", label=c("a","a","b"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15)),
ncol=4,nrow=2)

##Test group differences for Nest location.
kruskal.test(pheno_q90~Nest_location,data=survive)
DT=dunnTest(naLongMin~Nest_location,data=survive,method="bh"); DT

####Fig S6: NEST, plot differences across nesting strategy: --------------------------
plot_grid(
ggplot(data=subset(survive,!is.na(phenoRange)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=phenoRange) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("Pheno. Range (days)") + ylim(25,195) + scale_x_discrete(labels=c("cavity","clepto","ground") ) +
  stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(phenoMean)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=phenoMean) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("Pheno. Mean (days)") + ylim(50,330) + scale_x_discrete(labels=c("cavity","clepto","ground") ) +
  stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(ITD)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=ITD) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("ITD") + ylim(0.5,8) + scale_x_discrete(labels=c("cavity","clepto","ground") ) +
  stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(na_hull_area)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=log(na_hull_area)) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("log(NA Hull Area)") + ylim(11.5,20) + scale_x_discrete(labels=c("cavity","clepto","ground") ) +
  stat_summary(geom="text", label=c("a","b","c"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLatMax)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=naLatMax) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("NA Max Lat") + ylim(40,80) + scale_x_discrete(labels=c("cavity","clepto","ground") ) +
  stat_summary(geom="text", label=c("a","b","b"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLatMin)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=naLatMin) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("NA Min Lat") + ylim(8,55) + scale_x_discrete(labels=c("cavity","clepto","ground") ) +
  stat_summary(geom="text", label=c("a","b","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLongMax)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=naLongMax) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("NA Max Long") + ylim(-85,-25) + scale_x_discrete(labels=c("cavity","clepto","ground") ) +
  stat_summary(geom="text", label=c("a","a","a"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15))
,
ggplot(data=subset(survive,!is.na(naLongMin)&!is.na(Nest_location)),aes(x=as.factor(Nest_location),y=naLongMin) ) + theme_bw() +
  geom_boxplot(fill="light blue") + xlab("Nest") + ylab("NA Min Long") + ylim(-165,-55) + scale_x_discrete(labels=c("cavity","clepto","ground") ) +
  stat_summary(geom="text", label=c("a","b","c"),fun=max, vjust=-1, size=10) + theme(text=element_text(size=15)),
ncol=4,nrow=2)


###__________________________________________###
###SOME TRAIT CHANGES FROM 70s to 2010s#########  ##Figures S14 - S17 ------------------
###__________________________________________###
##Load necessary data:
EvansIsaacs=read.csv("C:/Users/prglaum/Documents/Pollinators-Empirical/CCmanuscript/pubAttempt2/ESGR_samples_EvansIsaacs_withTraits.csv",sep=",",header=T, na.strings=c("NaN","#NAME?","-Inf",""))

###__________________________________________###
###PHENOLOGY CHANGES FRMO 70s to 2010s##########  ##Fig 14 -----------------------------
###__________________________________________###
plot_grid(
ggplot(EvansIsaacs,aes(x=Collector,y=phenoRange,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Range (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("a) All samples"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="polylectic"),aes(x=Collector,y=phenoRange,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Range (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("b) Polylectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="oligolectic"),aes(x=Collector,y=phenoRange,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Range (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("c) Oligolectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="clepto-parasitic"),aes(x=Collector,y=phenoRange,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Range (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("d) Clepto"),

ggplot(EvansIsaacs,aes(x=Collector,y=phenoMean,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Mean (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("e) All samples"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="polylectic"),aes(x=Collector,y=phenoMean,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Mean (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("f) Polylectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="oligolectic"),aes(x=Collector,y=phenoMean,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Mean (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("g) Oligolectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="clepto-parasitic"),aes(x=Collector,y=phenoMean,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("Pheno. Mean (Jul. Days)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("h) Clepto"),
ncol=4,nrow=2)

kruskal.test(phenoRange~Collector,data=EvansIsaacs)
diets=subset(EvansIsaacs,Dietary_specialization=="polylectic") #clepto-parasitic
kruskal.test(phenoRange~Collector,data=diets)

kruskal.test(phenoMean~Collector,data=EvansIsaacs)
diets=subset(EvansIsaacs,Dietary_specialization=="clepto-parasitic") #clepto-parasitic
kruskal.test(phenoMean~Collector,data=diets)

###__________________________________________###
###ITD CHANGES FRMO 70s to 2010s: ##Fig 15 ----------------
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
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("Ground w/ no Bombus")
cleptoITD=ggplot(justClepto,aes(x=Collector,y=ITD,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("ITD(mm)") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("Clepto-parasitic")

create_plot()
plot_grid(allGround,allGroundnoBB,cleptoITD,ncol=3)


###__________________________________________###
###LATITUDE CHANGES FRMO 70s to 2010s########### ##Figure S16 ---------------
###__________________________________________###
plot_grid(
ggplot(EvansIsaacs,aes(x=Collector,y=naLatMax,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Max") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("a) All samples"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="polylectic"),aes(x=Collector,y=naLatMax,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Max") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("b) Polylectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="oligolectic"),aes(x=Collector,y=naLatMax,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Max") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("c) Oligolectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="clepto-parasitic"),aes(x=Collector,y=naLatMax,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Max") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("d) Clepto"),

ggplot(EvansIsaacs,aes(x=Collector,y=naLatMin,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Min") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("e) All samples"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="polylectic"),aes(x=Collector,y=naLatMin,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Min") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("f) Polylectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="oligolectic"),aes(x=Collector,y=naLatMin,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Min") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("g) Oligolectic"),
ggplot(subset(EvansIsaacs,Dietary_specialization=="clepto-parasitic"),aes(x=Collector,y=naLatMin,group=Collector)) + theme_bw() + geom_boxplot(fill="light blue") + ylab("NA Lat. Min") +
  theme(text = element_text(size=15)) + scale_x_discrete(labels=c("'72&'73","'17&'18") ) + ggtitle("h) Clepto"),
ncol=4,nrow=2)

kruskal.test(naLatMax~Collector,data=EvansIsaacs)
diets=subset(esgrNEW,Dietary_specialization=="clepto-parasitic") #clepto-parasitic
kruskal.test(naLatMax~Collector,data=diets)

kruskal.test(naLatMin~Collector,data=EvansIsaacs)
diets=subset(EvansIsaacs,Dietary_specialization=="polylectic") #clepto-parasitic
kruskal.test(naLatMin~Collector,data=diets)


###------------------###
##SPECIES LEVEL CHANGE## Fig S17 ----------------------
#______________________#

##Fig S17
##Phenological shift showing earlier flying SPECIES
ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=phenoRange,group=bin)) + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") )
ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=phenoMean,group=bin)) + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") )
ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=pheno_q10,group=bin)) + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") )
ggplot(subset(traits,bin>0),aes(x=factor(bin,levels=c(1,3,2)),y=pheno_q90,group=bin)) + geom_boxplot() + scale_x_discrete(labels=c("70s","Both","2010s") )



###_______________________###
## Graphical representation of phenologica range effects in polylectic bees Fig S20 ------------------
###_______________________###
library(sjPlot) #load package to use plot_model function
##Add the polylectic flag variable, "poly" to the survive data frame
survive$poly=rep("n",length(survive$MinLong))
survive$poly[survive$Dietary_specialization=="polylectic"]="y"

###_______________________###
## Figure 20a - Generalized linear model
mmod=(glmer(as.factor(extinct)~(1+(phenoRange)|DS)+(phenoRange),data=use,family=binomial))
summary(lmod)
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
gammod=gam(extinct ~ s(phenoRange, DS, bs="fs",k=4), data=use,family=binomial("logit"))
summary(gammod)
fam=family(gammod); ilink <- fam$linkinv;
pred=as.data.frame(predict(qw,newdata=use,type="link",se.fit=TRUE))

dfp=data.frame(pR=use$phenoRange,DS=use$DS,fit_resp  = ilink(pred$fit),
               right_upr = ilink(pred$fit + (2 * pred$se.fit)),right_lwr = ilink(pred$fit - (2 * pred$se.fit)))

ggplot(dfp,aes(x=pR,y=fit_resp,color=DS,fill=DS)) + geom_line(linewidth=1.5) + theme_bw() +
  geom_ribbon(aes(ymin = right_lwr, ymax = right_upr), alpha = 0.1) +
  guides(fill="none",color=guide_legend(title="Diet")) + theme(text = element_text(size = 14)) +
  xlab("Pheno. Range") + ylab("Prob. of Local Extirpation")


