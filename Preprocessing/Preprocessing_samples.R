
library(ChAMP)
library(LRTools)
library(ChAMPdata)
library(minfi)
library(akmedoids)
library(matrixStats)
library(pwr)
library(ggplot2)
library(illuminaio)
library(dplyr)

simu <-  "/Path/to/your/idat_simu"

myLoad_Simulation <- champ.load_extended(directory = simu,
                                    method = "minfi",
                                    methValue = "B",
                                    autoimpute = TRUE,
                                    filterDetP = TRUE,
                                    ProbeCutoff = 0,
                                    SampleCutoff = 0.1,
                                    detPcut = 0.01,
                                    filterBeads = TRUE,
                                    beadCutoff = 0.05,
                                    filterNoCG = TRUE,
                                    filterSNPs = TRUE,
                                    population = NULL,
                                    filterMultiHit = TRUE,
                                    filterXY = TRUE,
                                    force = FALSE,
                                    arraytype = "450K",
                                    sampleSheet.csv = "simulation.csv",
                                    preproc = "Noob",
                                    dyeMethod = "single")


colnames(myLoad_Simulation$mset) <- paste(myLoad_Simulation$mset$Sample_ID,myLoad_Simulation$mset$Group)

## re ordering the colums for later.
Beta_simu <- as.data.frame(getBeta(myLoad_Simulation$mset))
treat <- grep("Treatment", colnames(Beta_simu))
control <- grep("Control", colnames(Beta_simu))
new_order <- c(treat, control)

Beta_simu_reordered <- Beta_simu[,new_order]
colnames(Beta_simu_reordered)
rownames(Beta_simu_reordered)

champ.QC(beta = myLoad_Simulation$beta, pheno = myLoad_Simulation$pd$Group ,
         mdsPlot = TRUE, densityPlot = TRUE,
         dendrogram = TRUE, PDFplot = TRUE, 
         Rplot = TRUE, Feature.sel = "None",
         resultsDir = "/Path/to/your/idat_simu/QC_res")


write.csv(Beta_simu_reordered, file ="/Path/to/your/beta_simu.csv" )



champ.QC(beta = myLoad_Simulation$beta, pheno = myLoad_Simulation$pd$Group ,
         mdsPlot = TRUE, densityPlot = TRUE,
         dendrogram = TRUE, PDFplot = TRUE, 
         Rplot = TRUE, Feature.sel = "None",
         resultsDir = "/Path/to/your/idat_simu/QC_res")

kabuki <-  "/Path/to/your/idat_kabuki"

myload_Kabuki <- champ.load_extended(directory = kabuki ,
                                     method = "minfi",
                                     methValue = "B",
                                     autoimpute = TRUE,
                                     filterDetP = TRUE,
                                     ProbeCutoff = 0,
                                     SampleCutoff = 0.1,
                                     detPcut = 0.01,
                                     filterBeads = TRUE,
                                     beadCutoff = 0.05,
                                     filterNoCG = TRUE,
                                     filterSNPs = TRUE,
                                     population = NULL,
                                     filterMultiHit = TRUE,
                                     filterXY = TRUE,
                                     force = FALSE,
                                     arraytype = "450K",
                                     sampleSheet.csv = "KMT2D.csv",
                                     preproc = "Noob",
                                     dyeMethod = "single")

champ.QC(beta = myload_Kabuki$beta, pheno = myload_Kabuki$pd$disease.state.ch1 ,
         mdsPlot = TRUE, densityPlot = TRUE,
         dendrogram = TRUE, PDFplot = TRUE, 
         Rplot = TRUE, Feature.sel = "None",
         resultsDir = "/Path/to/your/idat_kabuki/QC_res")

