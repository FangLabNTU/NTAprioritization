# **NTAprioritization**
## **Summary**
NTAprioritization is an R package for analyzing raw candidate list obtained from non-target analysis. The raw candidate list can be obtained by pre-processing in Waters Porgenesis QI after acquisition by Waters LC-MS-QTOF. NTA prioritization enables a complete workflow from raw feature lists to a final prioritized identification list with multiple filters. The workflow starts with three parameters: spectral library searching scores, predicted retention time, and toxicity data. It includes 2 functions:
1. getRTMSMS():return candidate list with MSMS macthing score(fragmentation score) > 0 and RT-MS/MS level 
2. getTierlevel(): return data table in csv format with toxicity level and final priority in 5 tiers and another plot with highlighted tier 1 compounds.

## **Installation**
1. `NTAprioritization` requires R version equal or newer than 4.1.1. The following script allows to install rexposome dependencies:

`if (!requireNamespace("BiocManager", quietly=TRUE)){
		install.packages("BiocManager")}`
    
`packages = c('rcdk','rcdklibs','randomForest','leaps','caret','corrplot','tidyverse',
			'mlr','dplyr','Metrics','ggpubr','ggplot2','miceadds','rio','openxlsx','ggrepel')`
               
`for( pkg in packages ) {
		if( !pkg %in% rownames( installed.packages() ) ) {
				message( "Installing ", pkg )
					BiocManager::install( pkg )}
					}`


2. To install the package using: 

(1)Download the packages to local folder (path_to_file)

`install.packages(path_to_file, repos = NULL, type="source")`

(2) The package can be installed using the R package devtools. devtools can be installed win the following code:

`install.packages("devtools")`

`devtools::install_github("FangLabNTU/NTAprioritization")`

## **Data processing**
1. Get predicted retention time and RT-MS/MS level by typing getRTMSMS() in Rstudio console:

`getRTMSMS()`

interface pops up and direct to the source raw data in csv format

2. Get Toxicity level and Compound priority in 5 tiers by typing getRTMSMS() in Rstudio console:

`getTierlevel()`

refer to the folder where the raw data is stored and find the result table and plot.
