# pal-int-revised
  
This is the repository for the manuscript **Extinction risk controlled by interaction of long-term and short-term climate change** accepted for publication at *Nature Ecology and Evolution*.  
  
# Disclaimer  
  
Please note that we currently develop a R package to facilitate the calculation of palaeoclimate interactions. This package will make the code shown here in this repository more user-friendly and accessible. Please contact me via my [E-mail account](mailto:gregor.mathes@uni-bayreuth.de) if you have further questions or want early access to the package.  
  
# Scripts  

## Overview  
  
To reproduce our results, you either need to set the working directory for each R script to your local path or use the Rstudio project environment.  
We use the `here` package for reasons of workflow maintenance, so reading in data and saving results either requires the same folder structure as shown in this repository or local paths with each call.  
Other packages mainly used in our analysis are the `tidyverse` and `divDyn`.
  
### autocorrelation.R  
  
This script calculates the Durbin-Watson statistic for autocorrelation of each final GLMM of each fossil clade and produces Extended Data Figure 9b.  
  
  
### extended-data.R  
  
This script produces all table outputs directly within R using the `flextable` and the `officer` packages and exports them to a word file. The majority of Extended Data Figures and Supplemental Figures are produced with this script as well.  
  
### foram-preparation.R  
  
This script used the raw range-through data for foraminifera, cleans it and then processes it to produce the input for the final model in *glmm-analysis.R*. 
  
### glmm-analysis.R  
  
This script uses the processed extinction signal and climate proxy data for each fossil clade and calculates 10 GLMMs based on trends with varying lengths. 
The final model for each clade is then selected based on *AIC* and the effect intensity (*change in extinction risk*) is calculated based on this final model. 
The final output (model + effect) is then saved. Note that this script requires high computational power and might run for quite a while.  
  
### glmm-evaluation.R  
  
This script analyses each final GLMM for each fossil clade created in *glmm-analysis.R*. First, it calculates AIC and BIC values for models based on palaeoclimate interactions and for models based on short-term change only. 
These are then used for the model comparison and to produce Figure 2.  
Second, the effect intensity is extracted from the model output and combined with data from the null model simulation from *simulations-null-model.R* to produce Figure 3.  
Finally, the effect size is plotted against the duration of each fossil clade to produce Figure 4.  
  
### isotope-binning.R  
  
This script processes the raw climate proxy data from Veizer and Prokoph 2015. The raw data will then be further processed for merging with the extinction signal in *temperature-trends.R*.   
  
### mass-extinctions.R  
  
This script tests whether mass extinction bias our results. It runs a GLMM on all stages and one without the stages covering mass extinctions and then compares model quality via pseudo-R<sup>2</sup>.  
  
### now-preparation.R  
  
This script uses raw data downloaded from the NOW database for mammals, cleans it and then processes it to produce the input for the final model in *glmm-analysis.R*.  
  
### nr_species.R  
  
This script simply calculates the total number of distinct species per groups and sums them up.  
  
### pbdb-preparation.R  
  
This script shows how to download data from the Paleobiology database (PBDB) within R using the new API. It then cleans the data and processes it to produce the input for the final model in *glmm-analysis.R*.  
  
### simulations-autocorrelation.R  
  
This script conducts test whether autocorrelation between the extinction signal and climate proxy data could bias our results, using a simulation approach.
Note that this script requires high computational power and might run for quite a while.    
  
### simulations-null-model.R  
  
This script produces the null model to compare empirical results from *glmm-evaluation.R* with, using a simulation framework. 
Note that this script requires high computational power and might run for quite a while.   
  
### subsample-robustness-test.R  
  
This script conducts another robustness test, checking whether subsampling via classical rarefaction or shareholder quorum might change our results.  
  
### temperature-trends.R  
  
This script uses the processed climate proxy data and calculates short-term temperature change and long-term trends.  

# Folders  
  
## data  
  
This folder contains the processed data sets ready for GLMM calculaten. These are the files that generally end with *_trends.csv*. 
It further contains the file with geologic stage information *gradstein.csv* and files containing climate proxy data *Prokoph & Veizer 2015.csv*, *TimeSeriesUsed.csv*, and *isotemp_trends.csv*. 
This folder also contains subfolders:  
  
### model-output  
  
This subfolder contains the final model output (lists of the GLMM model and calculated effect). Unfortunately, size restrictions won't allow an upload of these. 
The models however can be reproduced using the R script *glmm-analysis.R*.  
  
### raw-fossil-data
  
This subfolder contains all raw and unprocessed data used in our analysis, for all 8 fossil clades.  
  
### results  
  
This subfolder contains all results from e.g. GLMM evaluation, simulations etc. saved in csv format.  
  
### simulation-results   
  
This subfolder contains the inbetween data from the simulations. Unfortunately, size restrictions won't allow an upload of these. 
They however can be reproduced using the R script *simulations-null-model.R* and *simulation-autocorrelation.R*.  
  
## figures  
  
This folder contains all figures shown in our manuscript. 
  
## tables  
  
This folder contains all tables for supplemental material and extended data as produced in *extended-data.R*. All tables are savedto word format and are exported to a single word file called *extended-tables.docx*.  
  