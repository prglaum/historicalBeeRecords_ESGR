###_______________________________________________________###
## Main text figures ##
## Focus species composition, diversity metrics, & trait analysis ##
###_______________________________________________________###

##All necessary data is located in the /data file at https://github.com/prglaum/historicalBeeRecords_ESGR/data &
##  https://agdatacommons.nal.usda.gov/account/articles/25233991

##Load libraries -----------------------
library(stringr) ##for modifying strings
library(ggplot2) ##for plotting
library(modEvA) #for Dsquared
library(mgcv) #for gams
library(dplyr) ##for data wrangling/manipulation
library(readxl) ##for loading xlsx files

## Contents -------
##  _________________  Species Composition Analysis ___________________
## Calculating rarefied species richness across 7 bins spanning 1921-2018. Line 29
## Species analysis between historical (1972/1973) and contemporary (2018/2019) major sampling periods. Line 60
## Figure 1a - Multivariate analysis and NMDS plot focused on species composition. Line 107

##  _________________  Trait Composition Analysis _____________________
## Figure 1b - Multivariate analysis and NMDS plot focused on trait composition. Line 152
## Figure 2 - See NNscript_analysis_plots_github.R script.
## Figure 3 - Categorical trait changes between major sample periods. Line 225


##Calculating rarefied species richness across bins ----------------------------

WBs <- read.csv(file.choose(),header=T) #bring in the file "ESGRbees_9.9.22_AllBees_UndetRemoved_ApisRemoved_RepeatsRemoved_Updated.csv"
head(WBs)

##Manipulate dataframe
colnames(WBs)
WBs_con<-WBs[c(1,7,10,11,21)] #condense dataframe to necessary columns only
head(WBs_con)
library(reshape2)
bees_group <- dcast(WBs_con, YearGroup ~ Updated_Genus_species, length) #change from long format to wide format, with species in columns, grouped by YearGroup (rows)
head(bees_group)
colnames(bees_group)
bees_g <- bees_group[2:222] #dataframe with only the species data

##Calculate species richness for each Year Grouping
library(vegan)
specrich <- specnumber(bees_g) #get number of species total for each grouping
specrich
info2 <- bees_group[1:1] #row information (year groupings)
info2
specrich2 <- cbind(info2, specrich) #connect the richness data to the row information
specrich2

#rarefy by YearGroup
head(bees_g)
rarefy(bees_g, se = TRUE) ##determine the smallest site maximum = 320
rare <- rarefy(bees_g,320, se=TRUE) #calculate rarefied species richness and standard error
rare


#####Historical (1972/1973) versus contemporary (2017/2018) analyses-----------------------------------------------

evans <- read.csv(file.choose(),header=T) #bring in the file "EvansIsaacs_ESGR_Apisremoved_Undetremoved_MayOctb.csv"
head(evans)
colnames(evans)
##calculate species richness for each group (historical=Evans or contemporary=Isaacs)
egroup2 <- dcast(evans, Collector ~ GenusSpecies, length)
egroup2$dummy <- rep(1,nrow(egroup2)) #include a dummy species for NMDS analysis
head(egroup2)
colnames(egroup2)
bees_e2 <- egroup2[2:168]
specrich22 <- specnumber(bees_e2) #get number of species total for each grouping
specrich22
info32 <- egroup2[1:1] #group info
head(info32)
specrich32 <- cbind(info32, specrich22) #connect the richness numbers with the group info
specrich32

##rarefy by collector grouping
head(bees_e2)
rarefy(bees_e2,se=TRUE) #get the smallest site maximum = 1272
rare22 <- rarefy(bees_e2,1272,se=TRUE)
rare22 #rarefied species richness and standard error


##calculate shannon diversity with Hill numbers
library(hilldiv)
beese_transpose <- t(bees_e2)
head(beese_transpose)
hill_div(beese_transpose,1) #q = 1 is for shannon diversity, group 1 = historical, group 2 = contemporary

#use Hutcheson t test to test for sig. difference between hill numbers
library(ecolTest)
evans2 <- as.numeric(bees_e2[1,])
isaacs2 <- as.numeric(bees_e2[2,])
# two-sided test
Hutcheson_t_test(x=evans2,
                 y=isaacs2)

#Calculate Pielou's Evenness for each group
H <- diversity(bees_e2) #calculate diversity
# Observed Richness
richness <- specnumber(bees_e2)  #calculate richness
evenness <- H/log(richness) # Pielou's Evenness
evenness

######NMDS ordination of species composition ----------------------------------------------------------

library(reshape2)
egroup <- dcast(evans, Collector+Month ~ GenusSpecies, length)
egroup$dummy <- rep(1,nrow(egroup))
head(egroup)
colnames(egroup)
bees_e <- egroup[3:169]
info3 <- egroup[1:2]
info3

library(vegan)
head(bees_e)
groups3 <- factor(c(rep("Evans",5),rep("Isaacs",5))) #designate groupings

dis2 <- vegdist(bees_e,method="bray") ##test for differences in dispersion between the two groups, if sig. you cannot run a PERMANOVA
mod2 <- betadisper(dis2, groups3)
anova(mod2) #variance is NOT sig different, proceed with the NMDS and PERMANOVA

ord2<-metaMDS(bees_e, autotransform=FALSE) #run NMDS
ord2
#plot NMDS
plot(ord2, disp='site', type="n", xlim=c(-2.5,2.5), ylim=c(-2.5,2.5), cex=0.8)
points(ord2, display="site",label = info3$Month, select=which(info3$Collector=="Evans"), pch=19, col="blue", cex=0.75)
points(ord2, display="site", select=which(info3$Collector=="Isaacs"), pch=19, col="green", cex=0.75)
levels(info3$Collector)=c("Evans", "Isaacs")
ordiellipse(ord2, info3$Collector, draw="polygon", border="blue", kind="sd", conf=0.95, label=FALSE, cex=0.75, show.groups="Evans")
ordiellipse(ord2, info3$Collector, draw="polygon", border="green", kind="sd", conf=0.95, label=FALSE, cex=0.75, show.groups="Isaacs")
legend(-3.5,2.5, title=NULL, pch=c(19,19), ncol=1, text.width=1,bty = "n",
       col=c("blue", "green"),
       cex=1, legend=c( "Historical", "Contemporary"))

ef <- envfit(ord2,bees_e, permu = 999) #calculate which species are most driving the ordination
ef
plot(ef, p.max = 0.01, cex = .5, col = "black") #fit vectors corresponding to species composition #sig set at p<0.01

#this is the Permanova test to see if the two communities are significantly different from each other
anova2 <- adonis2(bees_e ~info3$Collector, method="bray",permutations=999)
anova2

#	 ,-.		 ,-.
#	 \ /	 	 \ /
#	{|||)<	{|||)<        ________TRAIT BASED ANALYSES_________
#	 / \  	 / \ hjw
#	 `-^		 `-^
###__________________###
###TRAIT BASED NMDS  ###  ##Figure 1b ---------------
###__________________###

##Load necessary data
EvansIsaacs=read.csv("C:/Users/prglaum/Documents/Pollinators-Empirical/CCmanuscript/pubAttempt2/githubCode/ESGR_samples_EvansIsaacs_withTraitsMATCH.csv",sep=",",header=T, na.strings=c("NaN","#NAME?","-Inf",""))
  #read.csv("~/ESGR_samples_EvansIsaacs_withTraits.csv",sep=",",header=T, na.strings=c("NaN","#NAME?","-Inf",""))

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
anova3

##Make NMDS
znmds2 <- metaMDS(df[1:10,3:11], distance = "bray", k=2,na.rm=TRUE,autotransform = FALSE)
##Plot NMDS

#png("fig1b3.png",units="in", width=6.5, height=5.25,res=450)
plot(znmds2,type='n',xlim=c(-.2,.3),ylim=c(-.2,.2),cex=0.8)
points(znmds2, display = "sites",pch=20,col = c("blue", "green") [df$time]) #labels = (df$month[1:10])
ordiellipse(znmds2, groups = df$time[1:10], draw = "polygon", lty = 1,border=c("blue","green"),kind="sd", conf=0.95, label=FALSE, cex=0.75, show.groups=c("Evans","Isaacs")) #col = "grey90"
##include significant axes
fitter <- envfit(znmds2, df[1:10,3:11], permutations = 999)
plot(fitter, p.max = 0.01, col = "black", cex = 0.7) #p.max=0.005
##Add legend and text labels on data points
legend(.12,.22, title=NULL, pch=c(19,19), ncol=1, text.width=1,bty = "n",
       col=c("blue", "green"),
       cex=1, legend=c( "Historical", "Contemporary"))
text(znmds2$points[1:10], znmds2$points[11:20], labels=df$month[1:10], cex= 0.8,
 pos=c(1,3,3,1,4, 3,3,3,1,3))
#1=bottom, 2=left ,3=top, 4=right
#dev.off()

#                .' '.            __
#       .        .   .           (__\_
#        .         .         . -{{_(|8)
#jgs       ' .  . ' ' .  . '     (__/
#                .' '.            __
#       .        .   .           (__\_
#        .         .         . -{{_(|8)
#          ' .  . ' ' .  . '     (__/


###____________________________________###
####		SIMPLE CATEGORICAL COMPARISONS	###### Fig 3----------------------
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
## Figure 3a
NestSurvive=c(length(subset(survive,Nest_location=="ground"&survive$extinct==0)$species)/length(subset(survive,Nest_location=="ground")$species),
length(subset(survive,Nest_location=="cavity"&survive$extinct==0)$species)/length(subset(survive,Nest_location=="cavity")$species),
length(subset(survive,Nest_location=="clepto-parasitic"&survive$extinct==0)$species)/length(subset(survive,Nest_location=="clepto-parasitic")$species) )
barplot(NestSurvive,names=c("ground","cavity","clepto"),ylim=c(0,1),ylab="Percent of persistent species")

NestSurvive=data.frame(nest=c("ground","cavity","clepto"), persist=NestSurvive );
NestSurvive$nest=factor(NestSurvive$nest)
Fig3a=ggplot(data=NestSurvive,aes(y=persist,x=factor(nest,level=c("ground","cavity","clepto") ) )) +
  geom_bar(stat="identity",position='dodge') + theme_bw() + theme(text = element_text(size=20)) +
  labs(x = "", y = "Percent of persistent species") + ylim(0,1) + scale_x_discrete(labels=c("ground","cavity","klepto"))

## Figure 3d
dfNest=wtraits %>%
  group_by(bin) %>%
  count(Nest_location) %>%
  mutate(Freq = n/sum(n))
dfNest=dfNest[complete.cases(dfNest), ]
dfNest$bin=as.factor(dfNest$bin); dfNest$Nest_location=as.factor(dfNest$Nest_location);

Fig3d=ggplot(data=dfNest,aes(y=Freq,fill=factor(bin,level=c(1,3,2)),x=factor(Nest_location,level=c('ground','cavity','clepto-parasitic') ) )) +
  theme_bw() + theme(text = element_text(size=20)) +
  geom_bar(stat="identity",position='dodge') +
  labs(x = "Nest Type", y = "% of Raw Species Richness") + labs(fill = "Years") +
  scale_fill_discrete(labels=c("'72&'73","both","'17&'18")) +
  scale_x_discrete(labels=c("ground","cavity","klepto")) +
  theme(legend.position = c(.7,.75))

###DIETARY SPECIALIZATION
## Fig 3b
DietSurvive=c(length(subset(survive,Dietary_specialization=="polylectic"&survive$extinct==0)$species)/length(subset(survive,Dietary_specialization=="polylectic")$species),
length(subset(survive,Dietary_specialization=="oligolectic"&survive$extinct==0)$species)/length(subset(survive,Dietary_specialization=="oligolectic")$species),
length(subset(survive,Dietary_specialization=="clepto-parasitic"&survive$extinct==0)$species)/length(subset(survive,Dietary_specialization=="clepto-parasitic")$species) )
barplot(DietSurvive,names=c("polylectic","oligolectic","clepto"),ylim=c(0,1),ylab="Percent of persistent species")

DietSurvive=data.frame(diet=c("polylectic","oligolectic","clepto"), persist=DietSurvive );
DietSurvive$diet=factor(DietSurvive$diet)
Fig3b=ggplot(data=DietSurvive,aes(y=persist,x=factor(diet,level=c("polylectic","oligolectic","klepto") ) )) +
  geom_bar(stat="identity",position='dodge') + theme_bw() + theme(text = element_text(size=20)) +
  labs(x = "", y = "Percent of persistent species") + ylim(0,1) + scale_x_discrete(labels=c("polylectic","oligolectic","klepto"))

## Fig 3e
dfDiet=wtraits %>%
  group_by(bin) %>%
  count(Dietary_specialization) %>%
  mutate(Freq = n/sum(n))
dfDiet=dfDiet[complete.cases(dfDiet), ]
dfDiet$bin=as.factor(dfDiet$bin); dfDiet$Dietary_specialization=as.factor(dfDiet$Dietary_specialization);

Fig3e=ggplot(data=dfDiet,aes(y=Freq,fill=factor(bin,level=c(1,3,2)),x=factor(Dietary_specialization,level=c("polylectic","oligolectic","clepto-parasitic") ) )) +
  theme_bw() + theme(text = element_text(size=20)) +
  geom_bar(stat="identity",position='dodge') +
  labs(x = "Diet Type", y = "% of Raw Species Richness") + labs(fill = "Years") +
  scale_fill_discrete(labels=c("'72&'73","both","'17&'18")) +
  scale_x_discrete(labels=c("polylectic","oligolectic","klepto")) +
  theme(legend.position = c(.7,.75))

###SOCAILITY
## Fig 3c
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

SocSurvive=c(length(subset(survive,Sociality=="solitary"&survive$extinct==0)$species)/length(subset(survive,Sociality=="solitary")$species),
             length(subset(survive,Sociality=='Eusocial'&survive$extinct==0)$species)/length(subset(survive,Sociality=='Eusocial')$species),
             length(subset(survive,Sociality=="clepto-parasitic"&survive$extinct==0)$species)/length(subset(survive,Sociality=="clepto-parasitic")$species)
)
SocSurvive=data.frame(soc=c("solitary","eusocial","clepto"), persist=SocSurvive );
SocSurvive$soc=factor(SocSurvive$soc)
Fig3c=ggplot(data=SocSurvive,aes(y=persist,x=factor(soc,level=c("solitary","eusocial","clepto") ) )) +
  geom_bar(stat="identity",position='dodge') + theme_bw() + theme(text = element_text(size=20)) +
  labs(x = "", y = "Percent of persistent species") + ylim(0,1) + scale_x_discrete(labels=c("solitary","eusocial","klepto"))

## Fig 3f
dfSoc=wtraits %>% #subset(wtraits,Sociality!= "Facultative") %>%
  group_by(bin) %>%
  count(Sociality) %>%
  mutate(Freq = n/sum(n))
dfSoc=dfSoc[complete.cases(dfSoc), ]
dfSoc$bin=as.factor(dfSoc$bin); dfSoc$Sociality=as.factor(dfSoc$Sociality); dfSoc=subset(dfSoc,Sociality!="Facultative")

Fig3f=ggplot(data=dfSoc,aes(y=Freq,fill=factor(bin,level=c(1,3,2)),x=factor(Sociality,level=c("solitary","Eusocial","clepto-parasitic") ) )) +
  theme_bw() + theme(text = element_text(size=20)) +
  geom_bar(stat="identity",position='dodge') +
  labs(x = "Sociality", y = "% of Raw Species Richness") + labs(fill = "Years") +
  scale_fill_discrete(labels=c("'72&'73","both","'17&'18")) +
  scale_x_discrete(labels=c("solitary","eusocial","klepto")) +
  theme(legend.position = c(.7,.75))

FIG3=plot_grid(Fig3a, Fig3b, Fig3c,
          Fig3d, Fig3e, Fig3f,
          ncol=3,nrow=2,labels=c("a)","b)","c)","d)","e)","f)"),
          label_size=18,label_x = 0.9,label_y = 0.95)

ggsave(FIG3,file="Figure3.tiff", width=16, height=10)
