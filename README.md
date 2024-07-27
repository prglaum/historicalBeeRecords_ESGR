# historicalBeeRecords_ESGR
Analysis of historical bee community sampling in the E.S. George Reserve: an ecological preserve. Published as Graham et al 2024. The contents here are split into two folders:

1 - R scripts used in the analysis of Graham et al 2024

2 - Data used in the analysis of Graham et al 2024


# R scripts used in the analysis of Graham et al 2024:
All R scripts are commented with detailed instructions. 
We recommend reviewing NNscript_analysis_plots_github.R before investigating code/analysis presented in the
phylogeneticAutocorrelation.R script. This is because this later script utilizes analysis techniques introduced in the NNscript_analysis_plots_github.R script. 

AnalyticalScript_plots:
R script for all the figures from the main text.

NNscripts_analysis_plots:
R script detailing the neural network analysis used in the main text and supplementary material. 

SupplFigures:
R script with all the necessary code for the (non-neural network) figures in the Supplementary Results. 

geoPlotting:
R script detailing the mapping and GIS work done in Graham et al 2024. Also includes code to recreate Figure S26. 

phyloAutCorr:
R scripts detailing the use of Phylogenetic Eigenvector Mapping and other techniques used in analyzing phylogenetic autocorrelation. 

# Data used in the analysis of Graham et al 2024:
AllVettedVeganESGR.csv (Fig S11)
Sampled specimen data from the Evans historical (1972/1973) and Isaacs contemporary (2017/2018) sampling periods in the 
necessary format to rarify species richness using the vegan package in R.

ESGR_samples_EvansIsaacs_withTraits.csv (Fig S22-S24)
Sampled specimen data from the Evans historical (1972/1973) and Isaacs contemporary (2017/2018) sampling periods with 
species trait data attached to each individual specimen. 

Grahametal2023_BeeTraits.xlsx (all figures describing species traits)
Species taxonomic information, natural/life history traits, physiological traits, geographic distribution, etc. 
The first tab describes each trait, the second tab houses trait data per species. Bee Trait data is also availale on the Ag Data Commons: https://agdatacommons.nal.usda.gov/articles/dataset/Data_from_A_century_of_wild_bee_sampling_historical_data_and_neural_network_analysis_reveal_ecological_traits_associated_with_species_loss_/25233991

ML_nwk_version.txt:
Bee phylogeny used in analysis of phylogenetic autocorrelation. Phylogeny can be sourced:
Henríquez-Piskulich P, Hugall AF, Stuart-Fox D. 2024 A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila). 
Mol. Phylogenet. Evol. 190, 107963. (doi:10.1016/J.YMPEV.2023.107963) 
OR at: http://beetreeoflife.org/

Neural Network Olden values (Fig 2 & Figs S15,S17,S18,S19):
Necessary Re-scaled Olden (RSO) importance values from the unique training iterations of neural networks trained on North
American trait variables and extirpation/persistence data from the historical (1972/1973) and contemporary (2017/2018) 
sampling periods. 

na_OldenEvansIsaacs_allbees_1000iter.csv - Used to make Figure 2
na_polyC_OldenImp_1000rep.csv
na_polyG_OldenImp_1000rep.csv
na_oliG_OldenImp_1000rep.csv
na_clepto_OldenImp_1000rep.csv
Used to make figures S15,S17,S18,S19 in the Supplementary Results. 
Code to use the Olden data is located in the NNscript_analysis_plots_github.R script. 

RawGBIFdata_final.csv - File containing the raw dated georefenced data from GBIF. See citation below:
GBIF Occurrence Download https://doi.org/10.15468/dl.ta4zxp 
Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-04-26
Data was used to create phenology and geospatial profiles per bees species. Due to the large size, the actual file is available on Ag Data Commons: https://agdatacommons.nal.usda.gov/articles/dataset/Data_from_A_century_of_wild_bee_sampling_historical_data_and_neural_network_analysis_reveal_ecological_traits_associated_with_species_loss_/25233991
