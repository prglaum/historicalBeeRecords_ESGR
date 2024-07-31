# historicalBeeRecords_ESGR
#### Analysis of historical bee community sampling in the E.S. George Reserve: an ecological preserve. Published as Graham et al 2024. The contents here are split into two folders:

1 - R scripts used in the analysis of Graham et al 2024

2 - Data used in the analysis of Graham et al 2024

Further info at available at https://github.com/prglaum/historicalBeeRecords_ESGR & 
https://agdatacommons.nal.usda.gov/articles/dataset/Data_from_A_century_of_wild_bee_sampling_historical_data_and_neural_network_analysis_reveal_ecological_traits_associated_with_species_loss_/25233991


# R scripts used in the analysis of Graham et al 2024:
#### All analysis code is presented as R scripts from the R statistical language. All R scripts are commented with detailed instructions. All R scripts can be opened in base R (https://cran.r-project.org/), RStudio (https://posit.co/download/rstudio-desktop/), or any other IDE. <br>    We recommend reviewing NNscript_analysis_plots_github.R before investigating code/analysis presented in the phylogeneticAutocorrelation.R script. This is because this later script utilizes analysis techniques introduced in the NNscript_analysis_plots script. 

#### AnalyticalScript_plots:
R script for all the figures from the main text.

#### NNscripts_analysis_plots:
R script detailing the neural network analysis used in the main text and supplementary material. 

#### SupplFigures:
R script with all the necessary code for the (non-neural network) figures in the Supplementary Results. 

#### geoPlotting:
R script detailing the mapping and GIS work done in Graham et al 2024. Also includes code to recreate Figure S26. 

#### phyloAutCorr:
R scripts detailing the use of Phylogenetic Eigenvector Mapping and other techniques used in analyzing phylogenetic autocorrelation. 

# Data used in the analysis of Graham et al 2024:
#### Most data is presented in .csv or .xlsx formats. These are common file types which can be opened in Excel or via languages like R and Python. The bee phylogeny is a .txt file written in Newick format and should be opened via the phyloAuto script for our purposes. Each data file is followed by the figures it was used in creating and the R script(s) which call load that data. 

#### AllVettedVeganESGR.csv (Fig S11; SupplFigures script):
Sampled specimen data from the Evans historical (1972/1973) and Isaacs contemporary (2017/2018) sampling periods in the 
necessary format to rarify species richness using the vegan package in R. Each column heading represents a unique species and each row indicates a sample year with values indicating the number of that species sampled per year. This file is the input data required to run the rarefaction analysis using the vegan package in R (SupplFigures script). 

#### ESGR_samples_EvansIsaacs_withTraits.csv (Fig S22-S24; SupplFigures script):
Sampled specimen data from the Evans historical (1972/1973) and Isaacs contemporary (2017/2018) sampling periods with 
species trait data attached to each individual specimen. Each row represents a unique sampling of an individual bee. The species traits of each bee are described in every row it appears. See Grahametal2023_BeeTraits.xlsx's ReadMe sheet for a full description of the trait variables.  
This file is also availale on the Ag Data Commons: https://agdatacommons.nal.usda.gov/articles/dataset/Data_from_A_century_of_wild_bee_sampling_historical_data_and_neural_network_analysis_reveal_ecological_traits_associated_with_species_loss_/25233991

#### Grahametal2023_BeeTraits.xlsx (all figures describing species traits; called in every script):
Species taxonomic information, natural/life history traits, physiological traits, geographic distribution, etc. 
The .xlsx file contains two sheets. TraitData: all trait data per species. ReadMe: a mini "ReadMe" with all trait data columns explained. 
Bee Trait data is also availale as Table S2 as Supplementary Material and on the Ag Data Commons: https://agdatacommons.nal.usda.gov/articles/dataset/Data_from_A_century_of_wild_bee_sampling_historical_data_and_neural_network_analysis_reveal_ecological_traits_associated_with_species_loss_/25233991

#### ML_nwk_version.txt (Fig S7,S8, S13, S20; phyloAutCorr script):
Bee phylogeny used in analysis of phylogenetic autocorrelation in Newick formatting. Load into a session of R using the phyloAuroCorr script. Phylogeny can be sourced:
Henríquez-Piskulich P, Hugall AF, Stuart-Fox D. 2024 A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila). 
Mol. Phylogenet. Evol. 190, 107963. (doi:10.1016/J.YMPEV.2023.107963) 
OR at: http://beetreeoflife.org/

#### Neural Network Olden values (Fig 2 & Figs S15,S17,S18,S19; NNscript_analysis_plots script):
Necessary Re-scaled Olden (RSO) importance values from the unique training iterations of neural networks trained on North
American trait variables and extirpation/persistence data from the historical (1972/1973) and contemporary (2017/2018) 
sampling periods. <br>
**All Olden value data files have the following variables:** <br>
- x_names: The names of the input feature (trait variable in our case) evaluated for RSO Importance.<br>
- importance: The raw Olden importance value from the weights of one of the individual neural network out of the 1000 trained.<br>
- RSimportance: The Olden importance rescaled to be between -1 and 1 to be able to better compare across individually trained networks. _This is the variable used in the RSO Importance analysis (e.g., Fig 2)._ <br>
- order: The relative order of the trait variable's RSO Importance variable from each unique training instance out of the 1000 uniquely neural networks. There are 7 input trait features, so each trained network ordered traits w/in a range from 1 (most associated with extirpation) to 7 (most associated with persistence).  <br>
- rel_imp: Value for the Garson importance, an importance metric without direction only weight of predictive power. 

**Neural Network Olden Files:** <br>
- na_OldenEvansIsaacs_allbees_1000iter.csv - Olden data from neural networks trained on all bee data in extirpation analysis between 1972/1973 & 2017/2018. Used in NNscripts_analysis_plots script to make Figure 2. <br>
- na_polyC_OldenImp_1000rep.csv - Olden data from neural networks trained on polylectic cavity nesting bee data in extirpation analysis between 1972/1973 & 2017/2018. Used in NNscript_analysis_plots script to make Fig S15. <br>
- na_polyG_OldenImp_1000rep.csv - Olden data from neural networks trained on polylectic ground nesting bee data in extirpation analysis between 1972/1973 & 2017/2018. Used in NNscript_analysis_plots script to make Fig S16. <br>
- na_oliG_OldenImp_1000rep.csv - Olden data from neural networks trained on oligolectic ground nesting bee data in extirpation analysis between 1972/1973 & 2017/2018. Used in NNscript_analysis_plots script to make Fig S18. <br>
- na_clepto_OldenImp_1000rep.csv - Olden data from neural networks trained on kleptoparasitic bee data in extirpation analysis between 1972/1973 & 2017/2018. Used in NNscript_analysis_plots script to make Fig S19. <br>

#### RawGBIFdata_final.csv - File containing the raw dated georefenced data from GBIF. See citation below:
GBIF Occurrence Download https://doi.org/10.15468/dl.ta4zxp 
Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-04-26
Data was used to create phenology and geospatial profiles per bees species. See [GBIF](https://www.gbif.org/) for more on their data formatting. Due to the large size and licensing issues between GBIF and Dryad, this file is NOT uploaded here. The actual file is available on Ag Data Commons: https://agdatacommons.nal.usda.gov/articles/dataset/Data_from_A_century_of_wild_bee_sampling_historical_data_and_neural_network_analysis_reveal_ecological_traits_associated_with_species_loss_/25233991
