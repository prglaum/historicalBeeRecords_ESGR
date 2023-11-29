###_______________________________________________________###
## Main text figures focused on ##
## Simple trait analysis & richness rarefaction ##
###_______________________________________________________###

##All necessary data is located in the /data file at https://github.com/prglaum/historicalBeeRecords_ESGR

##Load libraries -----------------------
library(vegan) ##for rarefaction 
library(stringr) ##for modifying strings
library(ggplot2) ##for plotting
library(modEvA) #for Dsquared
library(mgcv) #for gams
library(dplyr) ##for data wrangling/manipulation
library(readxl) ##for loading xlsx files

## Contents -------
## Figure 1 - Rarefying richness data and plotting via Generalized Linear Models (GAM). Line 25
## Figure 2b - Multivariate analysis and NMDS plot. Line 90
## Figure 3 - See NNscript_analysis_plots_github.R script. 
## Figure 4 - Categorical trait changes between major sample periods. Line 159


###___________###
###RAREFACTION### ## Figure 1----------------------------
###___________###
library(vegan)

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

#	 ,-.		 ,-.
#	 \ /	 	 \ /
#	{|||)<	{|||)<
#	 / \  	 / \ hjw
#	 `-^		 `-^		
###__________________###
###TRAIT BASED NMDS  ###  ##Figure 2b ---------------
###__________________###

##Load necessary data 
EvansIsaacs=read.csv("~/ESGR_samples_EvansIsaacs_withTraits.csv",sep=",",header=T, na.strings=c("NaN","#NAME?","-Inf",""))

##Create data structure for the NMDS
##We remove max/min latitudes & longitudes due to negative numbers. Just use
##latitude and longitue ranges instead
tryer=structure(list(time=as.factor(EvansIsaacs$Collector),
  month=EvansIsaacs$Month,
  nest=EvansIsaacs$Nest_location,
  soc=EvansIsaacs$Sociality,
  diet=EvansIsaacs$Dietary_specialization,
  ITD=as.numeric(EvansIsaacs$ITD),
  latRange=as.numeric(EvansIsaacs$latRange),longRange=as.numeric(EvansIsaacs$longRange),
  phenoMean=as.numeric(EvansIsaacs$phenoMean),phenoRange=as.numeric(EvansIsaacs$phenoRange)),
  class = "data.frame", row.names = seq(1,length(EvansIsaacs$Collector)))

##Remove rows with NAs to remove species with incomplete trait data
tryer=tryer[complete.cases(tryer), ]

##separate by months
df=tryer %>%   
  group_by(time,month) %>%
  summarize(grPerc=sum(nest=='ground')/n(),cavPerc=sum(nest=='cavity')/n(),clepPerc=sum(nest=="clepto-parasitic")/n(),
  polyPerc=sum(diet=="polylectic")/n(),oligoPerc=sum(diet=="oligolectic")/n(),
  avgITD=mean(ITD),avgLatR=mean(latRange),avgPhM=mean(phenoMean),avgPhR=mean(phenoRange) )
df$month[1:10]=c("May","June","July","Aug","Sept","May","June","July","Aug","Sept")

##Check for similar dispersion
dispTest=df[1:10,3:11]
dis2 <- vegdist(dispTest,method="bray") ##test for differences in dispersion between the two sets of data, if sig. you cannot run a PERMANOVA
mod2 <- betadisper(dis2, df$time[1:10])
anova(mod2) #variance is NOT sig different.

##permutational multivariate analysis of variance using distance matrices
anova3 <- adonis(df[1:10,3:11] ~df$time[1:10], method="bray",permutations=999)
summary(anova3)

##Make NMDS
znmds2 <- metaMDS(df[1:10,3:11], distance = "bray", k=2,na.rm=TRUE,autotransform = FALSE)
##Plot NMDS
plot(znmds2,type='n',xlim=c(-.2,.3),ylim=c(-.2,.2),cex=0.8)
points(znmds2, display = "sites",labels = (df$month[1:10]),pch=20,col = c("blue", "green") [df$time])
ordiellipse(znmds2, groups = df$time[1:10], draw = "polygon", lty = 1,border=c("blue","green"),kind="sd", conf=0.95, label=FALSE, cex=0.75, show.groups=c("Evans","Isaacs")) #col = "grey90"
##include significant axes
fitter <- envfit(znmds2, df[1:10,3:11], permutations = 999)
plot(fitter, p.max = 0.005, col = "black", cex = 0.7)
##Add legend and text labels on data points
legend(.15,.2, title=NULL, pch=c(19,19), ncol=1, text.width=1,bty = "n",
       col=c("blue", "green"),
       cex=1, legend=c( "Historical", "Contemporary"))
text(znmds2$points[1:10], znmds2$points[11:20], labels=df$month[1:10], cex= 0.7, 
 pos=c(1,3,3,3,4, 3,3,3,1,3))
#1=bottom, 2=left ,3=top, 4=right


#                .' '.            __
#       .        .   .           (__\_
#        .         .         . -{{_(|8)
#jgs       ' .  . ' ' .  . '     (__/
#                .' '.            __
#       .        .   .           (__\_
#        .         .         . -{{_(|8)
#          ' .  . ' ' .  . '     (__/


###____________________________________###
####		SIMPLE CATEGORICAL COMPARISONS	###### Fig 4----------------------
###____________________________________###

##Load data 
traits=read_excel("~/Grahametal2023_BeeTraits.xlsx",sheet="TraitData",na = "NA")
##We use the "bin" variable to indicate species presence in across our two
##major sample periods (72/73 & 2017/2018). 
## Bin = 1: Species only found in 72/73
## Bin = 2: Species only found in 2017/2018
## Bin = 3: Species found in both 72/73 & 2017/2018

## Prepare data:
## Create dataframe only including species from major sample periods, 
## 72&73 & 2017/2018.
wtraits=subset(traits,bin>0)

###Creating few additional variables and caregorical variableswtraits$oliC=as.numeric(wtraits$Nest_location=="cavity"&wtraits$Dietary_specialization=="oligolectic")
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


### MAKING FIGURE PANELS ###
### NESTING LOCATION
## Figure 4a
NestSurvive=c(length(subset(survive,Nest_location=="ground"&survive$extinct==0)$species)/length(subset(survive,Nest_location=="ground")$species),
length(subset(survive,Nest_location=="cavity"&survive$extinct==0)$species)/length(subset(survive,Nest_location=="cavity")$species),
length(subset(survive,Nest_location=="clepto-parasitic"&survive$extinct==0)$species)/length(subset(survive,Nest_location=="clepto-parasitic")$species) )
barplot(NestSurvive,names=c("ground","cavity","clepto"),ylim=c(0,1),ylab="Percent of persistent species")

NestSurvive=data.frame(nest=c("ground","cavity","clepto"), persist=NestSurvive );
NestSurvive$nest=factor(NestSurvive$nest)
ggplot(data=NestSurvive,aes(y=persist,x=factor(nest,level=c("ground","cavity","clepto") ) )) + 
  geom_bar(stat="identity",position='dodge') + theme_bw() + theme(text = element_text(size=20)) +
  labs(x = "", y = "Percent of persistent species") + ylim(0,1)

## Figure 4d
dfNest=wtraits %>%
  group_by(bin) %>%
  count(Nest_location) %>%
  mutate(Freq = n/sum(n))
dfNest=dfNest[complete.cases(dfNest), ]
dfNest$bin=as.factor(dfNest$bin); dfNest$Nest_location=as.factor(dfNest$Nest_location); 

ggplot(data=dfNest,aes(y=Freq,fill=factor(bin,level=c(1,3,2)),x=factor(Nest_location,level=c('ground','cavity','clepto-parasitic') ) )) +
  theme_bw() + theme(text = element_text(size=20)) +
  geom_bar(stat="identity",position='dodge') +
  labs(x = "Nest Type", y = "% of Raw Species Richness") + labs(fill = "Years") + 
  scale_fill_discrete(labels=c("'72&'73","both","'17&'18")) +
  scale_x_discrete(labels=c("ground","cavity","clepto")) +
  theme(legend.position = c(.7,.75))

###DIETARY SPECIALIZATION
## Fig 4b
DietSurvive=c(length(subset(survive,Dietary_specialization=="polylectic"&survive$extinct==0)$species)/length(subset(survive,Dietary_specialization=="polylectic")$species),
length(subset(survive,Dietary_specialization=="oligolectic"&survive$extinct==0)$species)/length(subset(survive,Dietary_specialization=="oligolectic")$species),
length(subset(survive,Dietary_specialization=="clepto-parasitic"&survive$extinct==0)$species)/length(subset(survive,Dietary_specialization=="clepto-parasitic")$species) )
barplot(DietSurvive,names=c("polylectic","oligolectic","clepto"),ylim=c(0,1),ylab="Percent of persistent species")

DietSurvive=data.frame(diet=c("polylectic","oligolectic","clepto"), persist=DietSurvive );
DietSurvive$diet=factor(DietSurvive$diet)
ggplot(data=DietSurvive,aes(y=persist,x=factor(diet,level=c("polylectic","oligolectic","clepto") ) )) + 
  geom_bar(stat="identity",position='dodge') + theme_bw() + theme(text = element_text(size=20)) +
  labs(x = "", y = "Percent of persistent species") + ylim(0,1)

## Fig 4e
dfDiet=wtraits %>%
  group_by(bin) %>%
  count(Dietary_specialization) %>%
  mutate(Freq = n/sum(n))
dfDiet=dfDiet[complete.cases(dfDiet), ]
dfDiet$bin=as.factor(dfDiet$bin); dfDiet$Dietary_specialization=as.factor(dfDiet$Dietary_specialization); 

ggplot(data=dfDiet,aes(y=Freq,fill=factor(bin,level=c(1,3,2)),x=factor(Dietary_specialization,level=c("polylectic","oligolectic","clepto-parasitic") ) )) +
  theme_bw() + theme(text = element_text(size=20)) +
  geom_bar(stat="identity",position='dodge') +
  labs(x = "Diet Type", y = "% of Raw Species Richness") + labs(fill = "Years") + 
  scale_fill_discrete(labels=c("'72&'73","both","'17&'18")) +
  scale_x_discrete(labels=c("polylectic","oligolectic","clepto")) +
  theme(legend.position = c(.7,.75))

###SOCAILITY
## Fig 4c
SocSurvive=c(length(subset(survive,Sociality=="solitary"&survive$extinct==0)$species)/length(subset(survive,Sociality=="solitary")$species),
length(subset(survive,Sociality=='Eusocial'&survive$extinct==0)$species)/length(subset(survive,Sociality=='Eusocial')$species),
length(subset(survive,Sociality=="clepto-parasitic"&survive$extinct==0)$species)/length(subset(survive,Sociality=="clepto-parasitic")$species),
length(subset(survive,Sociality=="Facultative"&survive$extinct==0)$species)/length(subset(survive,Sociality=="Facultative")$species) )
barplot(SocSurvive,names=c("solitary","eusocial","clepto","facultative"),ylim=c(0,1),ylab="Percent of persistent species")

SocSurvive=data.frame(soc=c("solitary","eusocial","clepto","facultative"), persist=SocSurvive );
SocSurvive$soc=factor(SocSurvive$soc)
ggplot(data=SocSurvive,aes(y=persist,x=factor(soc,level=c("solitary","eusocial","clepto","facultative") ) )) + 
  geom_bar(stat="identity",position='dodge') + theme_bw() + theme(text = element_text(size=20)) +
  labs(x = "", y = "Percent of persistent species") + ylim(0,1)

## Fig 4e
dfSoc=wtraits %>%
  group_by(bin) %>%
  count(Sociality) %>%
  mutate(Freq = n/sum(n))
dfSoc=dfSoc[complete.cases(dfSoc), ]
dfSoc$bin=as.factor(dfSoc$bin); dfSoc$Sociality=as.factor(dfSoc$Sociality); 

ggplot(data=dfSoc,aes(y=Freq,fill=factor(bin,level=c(1,3,2)),x=factor(Sociality,level=c("solitary","Eusocial","clepto-parasitic","Facultative") ) )) +
  theme_bw() + theme(text = element_text(size=20)) +
  geom_bar(stat="identity",position='dodge') +
  labs(x = "Sociality", y = "% of Raw Species Richness") + labs(fill = "Years") + 
  scale_fill_discrete(labels=c("'72&'73","both","'17&'18")) +
  scale_x_discrete(labels=c("solitary","eusocial","clepto","facultative")) +
  theme(legend.position = c(.7,.75))