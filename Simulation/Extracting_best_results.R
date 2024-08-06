library(ChAMP)
library(ChAMPdata)
library(minfi)
library("ChIPpeakAnno")
library("limma")
library(DMRcate)
library(idDMR)
library(bumphunter)
library("DMRcatedata")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(dplyr)
library(plyr)
library(PRROC)


##### Data needed #### 
data.dir <-"/path/to/your/files/Simulation/"
Aclusters_df <- read.csv(paste0(data.dir, "A-clust-results.csv"),
                         header = TRUE)
rownames(Aclusters_df) <- Aclusters_df$cpg
beta_mat  <- read.csv(paste0(data.dir, "Beta_for_simu.csv"),
                      row.names = 1, header = TRUE)
CPGs_df <-  read.csv(paste0(data.dir, "cpgLocation_df.csv"),
                     row.names = 1,header = TRUE)

### Make sure that row names are the cpgnames 
rownames(Aclusters_df) 
rownames(beta_mat)


ResultProbelasso  <- read.csv( "/path/to/your/files/Simulation/Probelasso/clean_DMR_PB_Simulation.csv",
                      row.names = 1, header = TRUE) ### Beta_mat filtered. in A-Clust-results
ResultDMRcate  <- read.csv( "/path/to/your/files/Simulation/DMRCate/clean_DMR_cate_Simulation.csv",
                               row.names = 1, header = TRUE) ### Beta_mat filtered. in A-Clust-results
ResultidDMR  <- read.csv( "/path/to/your/files/Simulation/idDMR/clean_DMR_idDMR_Simulation.csv",
                               row.names = 1, header = TRUE) ### Beta_mat filtered. in A-Clust-results
ResultBumphunter  <- read.csv( "/path/to/your/files/Simulation/Bumphunter/clean_DMR_Bumphunter_Simulation.csv",
                               row.names = 1, header = TRUE) ### Beta_mat filtered. in A-Clust-results

bestResultsDir <- "/path/to/your/files/Simulation/Best_results/"

# Combine into a list

datasets <- list(ResultProbelasso, ResultDMRcate, ResultidDMR,ResultBumphunter)#ResultBumphunter

# Extract the row with the maximum AUPR value from each dataset

best_aupr_rows <- lapply(datasets, function(df) {
  max_index <- which.max(df$AuPR)
  df[max_index, ]
})

best_aupr_rows

rm(datasets,ResultProbelasso, ResultDMRcate, ResultidDMR,ResultBumphunter,best_aupr_rows)
####### best results for Probelasso 

source("/path/to/your/files/WriteProbeLassofunctions.R")
pVal_num_best <- 0.01
aveLassoRad_best <- 1000
minDmrSep_best <- 500

WriteProbeLassoResults(beta_mat,           
                         CPGs_df,Aclusters_df,
                       deltas_num = c(0.10, 0.25, 0.40),
                       seeds_int = c(110, 220, 340, 460, 690),
                         pVals_num =pVal_num_best,
                         aveLassoRad_int = aveLassoRad_best,
                         minDmrSep_int = minDmrSep_best,
                         resultsDir = bestResultsDir, 
                         verbose=TRUE)

out_PB <- ProcessProbeLassoResults(bestResultsDir,
                                   beta_mat,
                                   Aclusters_df,
                                   verbose = TRUE)


write.csv(out_PB, file = paste0(bestResultsDir, "best_PB_results.csv"))

rm(WriteProbeLassoResults,ProcessProbeLassoResults,RunProbeLasso,
   CleanResults, MergeDMRsWithCPGs, ProcessidDMRResults, 
   StandardizeOutput, pVal_num_best, aveLassoRad_best, minDmrSep_best, out_PB)

source("/path/to/your/files/WriteDMRCatefunctions.R")
lambda_best <- 1000
C_best <- 2 

WriteDMRcateResults(beta_mat,
                    CPGs_df,
                    Aclusters_df,
                    deltas_num = c(0.10,0.25, 0.40),
                    seeds_int = c(110, 220, 340, 460, 690),
                    lambdas_num = lambda_best,
                    Cs_int = C_best,
                    resultsDir = bestResultsDir,
                    verbose = TRUE)


out_DMRcate <- ProcessDMRcateResults(resultsDir = bestResultsDir,
                                beta_mat,
                                Aclusters_df,
                                verbose = TRUE)

write.csv(out_DMRcate, file = paste0(bestResultsDir, "best_dmrcate_results.csv"))

rm(WriteDMRcateResults,ProcessDMRcateResults, RunDMRcate,
   CleanResults, MergeDMRsWithCPGs, ProcessidDMRResults, 
    C_best, lambda_best, out_DMRcate)

source("/path/to/your/files/WriteBumphunterfunctions.R")

cutoffQ_best = 0.95
maxGap_best = 750

WriteBumphunterResults(beta_mat,
                       CPGs_df,
                       Aclusters_df,
                       deltas_num = c(0.10,0.25, 0.40),
                       seeds_int = c(110, 220, 340, 460, 690),
                       cutoffQ_num = cutoffQ_best,
                       maxGap_int = maxGap_best,
                       resultsDir = bestResultsDir,
                       verbose = TRUE)

out_Bumphunter <- ProcessBumphunterResults(bestResultsDir,
                                            beta_mat,
                                            Aclus_df,
                                            verbose = TRUE)
write.csv(out_Bumphunter, file = paste0(bestResultsDir, "best_Bumphunter_results.csv"))


rm(WriteBumphunterResults, ProcessBumphunterResults,
   CleanResults, MergeDMRsWithCPGs,  
   StandardizeOutput)

source("/path/to/your/files/WriteIdDMRfunctions.R")

G_best <- 1000
FDR_anno_best <- 0.05 

WriteidDMRResults(beta_mat,
                  CPGs_df,
                  Aclusters_df,
                  deltas_num = c(0.10,0.25, 0.40),
                  seeds_int = c(110, 220, 340, 460, 690),
                  G_int = G_best,
                  FDR_anno = FDR_anno_best,
                  resultsDir = bestResultsDir,
                  verbose = TRUE)

out_idDMR <- ProcessidDMRResults(bestResultsDir,
                                  beta_mat,
                                  Aclusters_df,
                                  verbose = TRUE)

write.csv(out_idDMR, file = paste0(bestResultsDir, "best_idDMR_results.csv"))
