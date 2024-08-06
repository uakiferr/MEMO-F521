library(ChAMP)
library(ChAMPdata)
library(minfi)
library("ChIPpeakAnno")
library("limma")
library(DMRcate)
library(idDMR)
library("DMRcatedata")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(dplyr)
library(plyr)
library(PRROC)


################ load all data needed for the analysis #######
data.dir <-"path/to/your/files/Kabuki/"

Kabuki_dmr_df <- read.csv(paste0(data.dir, "Kabuki_DMR_to_find_true.csv"), header = TRUE, row.names = 1)

Kabuki_dmr_df$actual <- c(rep("positive", length(nrow(Kabuki_dmr_df))))

rownames(Kabuki_dmr_df) <- Kabuki_dmr_df$cpg # making sure the row are the cgnames 

beta_kabuki  <- read.csv(paste0(data.dir, "beta_kabuki.csv"),
                         row.names = 1, header = TRUE) ### Beta_mat filtered. in A-Clust-results
CPGs_df <- read.csv("path/to/your/files/Simulation/cpgLocation_df.csv",
                    row.names = 1, header = TRUE) ## based on the DMRcompare, see Aclust-results. 

resultsDir <- "path/to/your/files/Kabuki/Best_results_kabuki/"

ResultProbelasso  <- read.csv( "path/to/your/files/Kabuki/Probelasso/clean_DMR_PB_Kabuki.csv",
                               row.names = 1, header = TRUE) ### Beta_mat filtered. in A-Clust-results
ResultDMRcate  <- read.csv( "path/to/your/files/Kabuki/DMRCate/clean_DMR_cate_Kabuki.csv",
                            row.names = 1, header = TRUE) ### Beta_mat filtered. in A-Clust-results
ResultidDMR  <- read.csv( "path/to/your/files/Kabuki/idDMR/clean_DMR_idDMR_Kabuki.csv",
                          row.names = 1, header = TRUE) ### Beta_mat filtered. in A-Clust-results

source("path/to/your/files/Kabuki/Function_Needed_for_best_results.R")

datasets <- list(ResultProbelasso, ResultDMRcate, ResultidDMR)

# Extract the row with the maximum AUPR value from each dataset

best_F1_rows <- lapply(datasets, function(df) {
  filtered_df <- df[df$size == 10, ]  # Filter rows where size == 3, 5, 7, 10 
  if (nrow(filtered_df) > 0) {  # Check if there are any rows left after filtering
    max_index <- which.max(filtered_df$F1)
    filtered_df[max_index, ]
  } else {
    NULL  # Return NULL if no rows match the condition
  }
})
best_F1_rows <- (Filter(Negate(is.null), best_F1_rows))
best_F1_rows

best_aupr_rows <- lapply(datasets, function(df) {
  filtered_df <- df[df$size == 10, ]  # Filter rows where size == 3, 5, 7, 10 
  if (nrow(filtered_df) > 0) {  # Check if there are any rows left after filtering
    max_index <- which.max(filtered_df$AuPR)
    filtered_df[max_index, ]
  } else {
    NULL  # Return NULL if no rows match the condition
  }
})
best_aupr_rows <- Filter(Negate(is.null), best_aupr_rows)
best_aupr_rows

best_TP_rows <- lapply(datasets, function(df) {
  filtered_df <- df[df$size == 10, ]  # Filter rows where size == 3, 5, 7, 10 
  if (nrow(filtered_df) > 0) {  # Check if there are any rows left after filtering
    max_index <- which.max(filtered_df$TP)
    filtered_df[max_index, ]
  } else {
    NULL  # Return NULL if no rows match the condition
  }
})

best_TP_rows <- Filter(Negate(is.null), best_TP_rows)
best_TP_rows

rm(datasets,ResultProbelasso, ResultDMRcate, ResultidDMR,ResultBumphunter,best_aupr_rows)


### Probelasso 
adjPval_best <- c(0.01, 0.05)
mLassoRad_best <- c(700, 1000) 
minDmrSep_best <- c(200, 750) 


WriteProbeLassoResults(beta_kabuki,           
                       CPGs_df,
                       sizes_num = 11,
                       seeds_int = 50,
                       pVals_num = adjPval_best,
                       aveLassoRad_int = mLassoRad_best,
                       minDmrSep_int =minDmrSep_best,
                       resultsDir = resultsDir, 
                       verbose=TRUE)


### idDMR 

G <- c(200, 250, 1000)
FDR <- c(0.001, 0.01, 0.1)

WriteidDMRResults(beta_kabuki,
                  CPGs_df,
                  sizes_num = 11 ,
                  seeds_int =50,
                  G_int =  G,
                  FDR_anno = FDR,
                  resultsDir = resultsDir,
                  verbose = TRUE)


### DMRcate 

lambda <- 1000
C_int <- c(1,2,5) 


WriteDMRcateResults(beta_kabuki,
                    CPGs_df,
                    sizes_num = 11,
                    seeds_int = c(50, 100, 210, 330, 450),
                    lambdas_num = lambda,
                    Cs_int = C_int,
                    resultsDir = resultsDir ,
                    verbose = TRUE)

################# Getting the results overlapps and PR curve ##############


DMR_cate_Kabuki <- ProcessDMRcateResults(resultsDir= resultsDir,
                                       beta_kabuki,
                                       Kabuki_dmr_df,
                                       verbose = TRUE)

PbL_Kabuki <- ProcessProbeLassoResults(resultsDir,
                                beta_kabuki,
                                Kabuki_dmr_df,
                                verbose = TRUE)

idDMR_Kabuki <- ProcessidDMRResults(resultsDir= resultsDir,
                                    beta_kabuki,
                                    Kabuki_dmr_df,
                                    verbose = TRUE)




res_best <- list(PbL_Kabuki, idDMR_Kabuki, DMR_cate_Kabuki)
####### Best results ########

best_F1 <- lapply(res_best , function(df) {
  max_index <- which.max(df$F1)
  df[max_index, ]
})

best_F1

best_TP <- lapply(res_best , function(df) {
  max_index <- which.max(df$TP)
  df[max_index, ]
})

best_TP

best_aupr <- lapply(res_best , function(df) {
  max_index <- which.max(df$AuPR)
  df[max_index, ]
})

ranges_df <- idDMRResults_samplesize11_seed50_G1000_FDRcpg0.001[[1]]
idDMR <-ranges_df[ranges_df$dmr.n.cpgs > 2 & ranges_df$dmr.pval < 0.05 & abs(ranges_df$maxdiff)> 0.1, ]
ranges_df <- ProbeLassoResults_samplesize11_seed50_adjPvalProbe0.05_meanLassoRd700_minDmrSep200[[1]]
PB <- ranges_df[ranges_df$dmr.n.cpgs > 2 & ranges_df$dmr.pval < 0.05 & abs((ranges_df$betaAv_Control)-(ranges_df$betaAv_Kabuki)) > 0.1,  ]
ranges_df <- DMRcateResults_samplesize11_seed50_lambda1000_C5[[1]]
DMRc <-ranges_df[ranges_df$dmr.n.cpgs > 2 & ranges_df$dmr.pval < 0.05 & abs(ranges_df$maxdiff)> 0.1, ]

total <- sum(idDMR$no.cpgs) + sum(PB$dmr.n.cpgs) +sum(DMRc$no.cpgs)

write.csv(idDMR, file = paste0(BestResultsDir, "res_idDMR_best_DMRs.csv")) 
write.csv(DMRc, file = paste0(BestResultsDir, "res_DMRcate_best_DMRs.csv"))
write.csv(PB, file = paste0(BestResultsDir, "res_ProbeLasso_best_DMRs.csv"))
