## Cran package
install.packages(qqman)
install.packages(ggplot2)
install.packages(plotly)
install.packages('gmodels')
install.packages(scales)
install.packages(reshape2)
install.packages(viridis)
install.packages(dplyr)
install.packages(purrr)
install.packages(stringr)
install.package(data.table)
install.packages(matrixStats)
install.packages(pwr)
install.packages(ggplot2)
install.packages(c("WriteXLS","MASS","dplR"),repos="http://cran.r-project.org")

## Biocmanager package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c( 
  "ChAMP", 
  "DMRcate",
  "minfi",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylationEPICmanifest",
  "wateRmelon",
  "readxl",
  "RPMM",
  "ChAMPdata",
  "DMRcatedata",
  "ChIPpeakAnno",
  "minfiData",
  "FlowSorted.Blood.450k",
  "FlowSorted.Blood.EPIC",
  "FDb.InfiniumMethylation.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "geneLenDataBase",
  "illuminaio",
  "minfiDataEPIC",
  "preprocessCore",
  "limma",
  "bioDist"),update = TRUE, ask = FALSE, quiet = TRUE)


## from github sources

remotes::install_github("MAnalytics/akmedoids") #needed for LRTools
devtools::install_github("LionelRohner/LRTools")
remotes::install_github("gabrielodom/DMRcomparePaper")
devtools::install_github('tamartsi/Aclust', dependencies = TRUE)
devtools::install_github('sjczheng/EpiDISH', dependencies = TRUE)

## For the AClust_import files based on DMRCompare -> github
install.packages("IMA",repos=c("http://rforge.net")) 


## Instead, the region-level annotation library for the 450k 
##microarray could be downloaded from here : https://rforge.net/IMA/fullannotInd.rda

##Then users can load the regional-level annotation library by 
##issuing the following command at the R prompt:
  
load("./fullannotInd.rda") 

