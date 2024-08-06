##### All functions needed for probelasso in order to generate the results #####

SimulateData <- function(beta_mat, ## the normalized beta_mat that had no probes to much over 0.95 ect.. 
                         Aclusters_df, ## Clusters defined in Aclust_import
                         delta_num, ## treatment added
                         seed_int, ## choose random seed from specified vector
                         betaCols_idx = 9:49, ## from 9 to 49 is beta_value of samples
                         numEx_int = 21, numClusters_int = 500){
  
  ###  Setup  ###
  # set the seed value from the 'rp'-th element from the vector 'seed_values'
  set.seed(seed_int)
  
  # randomly pick numClusters_int clusters
  clusts_int <- Aclusters_df$Clusternumber
  randClusts_int <- sample(max(clusts_int), numClusters_int)
  
  # full data info of randomly picked clusters
  inflateClust_df <- Aclusters_df[which(clusts_int %in% randClusts_int), ]
  bval_mat <- as.matrix(inflateClust_df[, betaCols_idx])
  
  
  ### Check average b-val/m-val group-wise in each cpg, and add delta value ###
  if (delta_num > 0) {
    
    # for all cpgs of numClusters_int random picked clusters
    for (i in 1:nrow(bval_mat)) {
      
      avgbetagr1 <- mean(bval_mat[i, 1:numEx_int])
      avgbetagr2 <- mean(bval_mat[i, (numEx_int + 1):ncol(bval_mat)])
      
      # Switch statement to ensure that adding delta only increases the signal
      if (avgbetagr1 > avgbetagr2) {
        bval_mat[i, 1:numEx_int] <- bval_mat[i, 1:numEx_int] + delta_num
      } else {
        
        bval_mat[i, (numEx_int + 1):ncol(bval_mat)] <-
          bval_mat[i, (numEx_int + 1):ncol(bval_mat)] + delta_num
        
      }
      
    } # END for()
    
    # If any entries of the resultant bval_mat, after adding delta_num, are
    #   greater than 1, replace them by 1 (the max beta score)
    bval_mat[bval_mat >= 1] <- 0.999
    
  } # END if()
  bval_df <- as.data.frame(bval_mat)
  
  inflateClust_df$actual <- "positive"
  
  
  ###  The other clusters  ###
  # Full data info of remaining CPGs, i.e. except randomly picked clusters
  nonInflateClust_df <-
    Aclusters_df[!(rownames(Aclusters_df) %in% inflateClust_df$cpg), ]
  nonInflateClust_df$actual <- "negative"
  
  # Full data info (beta values) of remaining CPGs, i.e. except randomly picked
  #   clusters
  nonInflateBvals_mat <-
    beta_mat[!(rownames(beta_mat) %in% inflateClust_df$cpg), ]
  
  
  ###  Combine Treated and Untreated Clusters  ###
  # Combine the CPGs (belonging to selected clusters) and remianing cpgs
  #   (belonging to non-selected clusters)
  treatedAclusters_df <- rbind(inflateClust_df, nonInflateClust_df)
  # Combine the treated (delta-value added) CPGs (belonging to the clusters
  #   selected at random) and untreated (delta-value not added) CPGs (belonging
  #   to the non-selected clusters)
  betaResults_df <- rbind(bval_df, nonInflateBvals_mat)
  
  
  ###  Order Columns and Rows  ###
  # Reorder the 'actual' column to previous of beta values
  treatedAclusters_df <- as.data.frame(treatedAclusters_df)
  impVars_df <- treatedAclusters_df[c("Clusternumber", "cpg", "CHR", "MAPINFO",
                                      "start_position", "end_position",
                                      "coordinate_37", "chromosome", "actual")]
  otherVars_df <- treatedAclusters_df[setdiff(names(treatedAclusters_df),
                                              names(impVars_df))]
  treatedAclusters_df <- cbind(impVars_df, otherVars_df)
  
  # Reorder the rows (CPGs) according to clusternumber
  treatedAclustCPGordered_df <-
    treatedAclusters_df[order(treatedAclusters_df$Clusternumber), ]
  
  if (delta_num == 0) {
    treatedAclustCPGordered_df$actual <- "negative"
  }
  
  
  ###  Return  ###
  list(simBetaVals_df = betaResults_df,
       simAclusters_df = treatedAclustCPGordered_df)
  
}


MergeDMRsWithCPGs <- function(DMRs_df, CPGs_df, alpha = 0.05){
  
  ###  1. Make ranges of DMR info  ###
  DMRs_df <- DMRs_df[DMRs_df$dmr.pval < alpha , ]
  sig.ranges <- IRanges(DMRs_df$dmr.start, DMRs_df$dmr.end)
  dmr.ranges <- GRanges(seqnames = DMRs_df$dmr.chr, ranges = sig.ranges)
  
  
  ###  2. Make ranges of CPG info  ###
  temp.cpg.ranges <- IRanges(CPGs_df$MAPINFO, CPGs_df$MAPINFO)
  cpg.ranges <- GRanges(seqnames = CPGs_df$chr, ranges = temp.cpg.ranges)
  
  
  ###  3. Find overlaps  ###
  overlaps_df <- as.data.frame(
    findOverlaps(query = dmr.ranges, subject = cpg.ranges, type = "any")
  )
  
  overlaps_df$dmr.order <- as.numeric(overlaps_df$queryHits)
  overlaps_df$cpg.order <- as.numeric(overlaps_df$subjectHits)
  
  
  ###  4. Merge with DMR info  ###
  DMRs_df$dmr.row <- 1:nrow(DMRs_df)
  DMRsInfo_df <- merge(x = overlaps_df, y = DMRs_df,
                       by.x = "dmr.order", by.y = "dmr.row")
  
  
  ###  5. Merge with CPGs info  ###
  CPGs_df$row <- 1:nrow(CPGs_df)
  overlapInfo_df <- merge(DMRsInfo_df, CPGs_df,
                          by.x = "cpg.order", by.y = "row")
  overlapInfo_df <- overlapInfo_df[order(overlapInfo_df$dmr.order), ]
  
  
  ###  6. Add number of CPGs and Return  ###
  numCPGs <- as.data.frame(table(overlapInfo_df$dmr.order))
  
  merge(overlapInfo_df, numCPGs,
        by.x = "dmr.order", by.y = "Var1")
  
}


StandardizeOutput <- function(methodOut_df,
                              method = c("DMRcate",
                                         "ProbeLasso",
                                         "Bumphunter"),
                              cpgLocation_df,
                              dmr.sig.threshold = 0.05,
                              min.cpgs = 5){
  
  method <- match.arg(method)
  
  ###  Extract Output  ###
  switch(method,
         
         DMRcate = {
           
           methodOut_df$dmr.pval  <- methodOut_df$Stouffer
           methodOut_df$dmr.chr   <- methodOut_df$seqnames
           methodOut_df$dmr.start <- methodOut_df$start
           methodOut_df$dmr.end   <- methodOut_df$end
           
         },
         ProbeLasso = {
           
           methodOut_df$dmr.chr   <- methodOut_df$seqnames
           methodOut_df$dmr.start <- methodOut_df$start
           methodOut_df$dmr.end   <- methodOut_df$end
           methodOut_df$dmr.pval  <- methodOut_df$dmrP
           
         },
         Bumphunter = {
           
           methodOut_df$dmr.pval  <- methodOut_df$p.valueArea
           methodOut_df$dmr.chr   <- paste0("chr", as.character(methodOut_df$chr))
           methodOut_df$dmr.start <- methodOut_df$start
           methodOut_df$dmr.end   <- methodOut_df$end
           
         }
         
  )
  
  
  ###  Transform  ###
  temp <- MergeDMRsWithCPGs(DMRs_df = methodOut_df,
                            CPGs_df = cpgLocation_df,
                            alpha = dmr.sig.threshold)
  temp.ncpgs <- unique(
    subset(temp, select = c("dmr.chr", "dmr.start", "dmr.end", "Freq"))
  )
  
  clean_df <-
    merge(methodOut_df, temp.ncpgs, by = c("dmr.chr", "dmr.start", "dmr.end"))
  clean_df$dmr.n.cpgs <- clean_df$Freq
  clean_df$Freq <- NULL
  clean_df <- clean_df[clean_df$dmr.n.cpgs >= min.cpgs, ]
  
  if(nrow(clean_df) == 0){
    clean_df <- NULL
  }
  
  
  ###  Return  ###
  clean_df
  
}



RunProbeLasso <- function(betaVals_mat,
                          labels_fct = factor(c(rep("Treatment", 21),
                                                rep("control", 20))),
                          cpgLocation_df,
                          adjPvalProbe_num,
                          meanLassoRadius_int,
                          minDmrSep_int,
                          dmr.sig.threshold = 0.05,
                          min.cpgs = 5) {
  
  ###  Calculate ProbeLasso Results  ###
  ptm <- proc.time()
  # Requires ChAMPdata::probe.features
  probeLasso_out <- tryCatch(
    
    suppressMessages(champ.DMR(beta = betaVals_mat,
                               pheno = labels_fct,
                               method = "ProbeLasso",
                               minProbes = 3,
                               adjPvalDmr = 0.05,
                               meanLassoRadius = meanLassoRadius_int,
                               minDmrSep = minDmrSep_int,
                               adjPvalProbe = adjPvalProbe_num,
                               PDFplot = FALSE,
                               Rplot = FALSE) ),
    error = function(e1){ NULL }
    
  )
  
  
  
  elapsedtime <- proc.time() - ptm
  
  ####Called in StandardizeOutput 
  MergeDMRsWithCPGs <- function(DMRs_df, CPGs_df, alpha = 0.05){
    
    ###  1. Make ranges of DMR info  ###
    DMRs_df <- DMRs_df[DMRs_df$dmr.pval < alpha , ]
    sig.ranges <- IRanges(DMRs_df$dmr.start, DMRs_df$dmr.end)
    dmr.ranges <- GRanges(seqnames = DMRs_df$dmr.chr, ranges = sig.ranges)
    
    
    ###  2. Make ranges of CPG info  ###
    temp.cpg.ranges <- IRanges(CPGs_df$MAPINFO, CPGs_df$MAPINFO)
    cpg.ranges <- GRanges(seqnames = CPGs_df$chr, ranges = temp.cpg.ranges)
    
    
    ###  3. Find overlaps  ###
    overlaps_df <- as.data.frame(
      findOverlaps(query = dmr.ranges, subject = cpg.ranges, type = "any")
    )
    
    overlaps_df$dmr.order <- as.numeric(overlaps_df$queryHits)
    overlaps_df$cpg.order <- as.numeric(overlaps_df$subjectHits)
    
    
    ###  4. Merge with DMR info  ###
    DMRs_df$dmr.row <- 1:nrow(DMRs_df)
    DMRsInfo_df <- merge(x = overlaps_df, y = DMRs_df,
                         by.x = "dmr.order", by.y = "dmr.row")
    
    
    ###  5. Merge with CPGs info  ###
    CPGs_df$row <- 1:nrow(CPGs_df)
    overlapInfo_df <- merge(DMRsInfo_df, CPGs_df,
                            by.x = "cpg.order", by.y = "row")
    overlapInfo_df <- overlapInfo_df[order(overlapInfo_df$dmr.order), ]
    
    
    ###  6. Add number of CPGs and Return  ###
    numCPGs <- as.data.frame(table(overlapInfo_df$dmr.order))
    
    merge(overlapInfo_df, numCPGs,
          by.x = "dmr.order", by.y = "Var1")
    
  }
  #### Standardize output 
  StandardizeOutput <- function(methodOut_df,
                                method = c("DMRcate",
                                           "ProbeLasso",
                                           "Bumphunter",
                                           "Comb_p"),
                                cpgLocation_df,
                                dmr.sig.threshold = 0.05,
                                min.cpgs = 5){
    
    method <- match.arg(method)
    
    ###  Extract Output  ###
    switch(method,
           
           DMRcate = {
             
             methodOut_df$dmr.pval  <- methodOut_df$Stouffer
             methodOut_df$dmr.chr   <- methodOut_df$seqnames
             methodOut_df$dmr.start <- methodOut_df$start
             methodOut_df$dmr.end   <- methodOut_df$end
             
           },
           ProbeLasso = {
             
             methodOut_df$dmr.chr   <- methodOut_df$seqnames
             methodOut_df$dmr.start <- methodOut_df$start
             methodOut_df$dmr.end   <- methodOut_df$end
             methodOut_df$dmr.pval  <- methodOut_df$dmrP
             
           },
           Bumphunter = {
             
             methodOut_df$dmr.pval  <- methodOut_df$p.valueArea
             methodOut_df$dmr.chr   <- paste0("chr", as.character(methodOut_df$chr))
             methodOut_df$dmr.start <- methodOut_df$start
             methodOut_df$dmr.end   <- methodOut_df$end
             
           }
           
    )
    
    
    ###  Transform  ###
    temp <- MergeDMRsWithCPGs(DMRs_df = methodOut_df,
                              CPGs_df = cpgLocation_df,
                              alpha = dmr.sig.threshold)
    temp.ncpgs <- unique(
      subset(temp, select = c("dmr.chr", "dmr.start", "dmr.end", "Freq"))
    )
    
    clean_df <-
      merge(methodOut_df, temp.ncpgs, by = c("dmr.chr", "dmr.start", "dmr.end"))
    clean_df$dmr.n.cpgs <- clean_df$Freq
    clean_df$Freq <- NULL
    clean_df <- clean_df[clean_df$dmr.n.cpgs >= min.cpgs, ]
    
    if(nrow(clean_df) == 0){
      clean_df <- NULL
    }
    
    
    ###  Return  ###
    clean_df
    
  }
  
  ###  Extract Results  ###
  # extract results if any predicted cluster is identified from ProbeLasso
  if(!is.null(probeLasso_out)){
    
    probeLassoOut_df <- probeLasso_out$ProbeLassoDMR
    
    results_df <- StandardizeOutput(
      methodOut_df = probeLassoOut_df,
      method = "ProbeLasso",
      cpgLocation_df = cpgLocation_df,
      dmr.sig.threshold = dmr.sig.threshold,
      min.cpgs = min.cpgs
    )
    
  } else {
    results_df <- NULL
  }
  
  ###  Return  ###
  list(results_df, elapsedtime[3])
  
}


WriteProbeLassoResults <- function(beta_mat,
                                   CPGs_df,
                                   Aclusters_df,
                                   deltas_num = c(0, 0.10, 0.20, 0.30, 0.40),
                                   seeds_int = c(100, 210, 330, 450, 680),
                                   pVals_num = c(0.001, 0.01, 0.05, 0.1),
                                   aveLassoRad_int = c(375, 700, 1000),
                                   minDmrSep_int = c(200, 250, 500, 750, 1000),
                                   resultsDir = "ProbeLasso_compare/",
                                   verbose = TRUE){

  dir.create(resultsDir, showWarnings = TRUE)
  
  ###  Data Simulation Outer Loop  ###
  designPts_mat <- expand.grid(deltas_num, seeds_int)
  paramsGrid_mat <- expand.grid(pVals_num, aveLassoRad_int, minDmrSep_int)
  
  for(i in 1:nrow(designPts_mat)){
    
    ###  Generate Data Set  ###
    delta <- designPts_mat[i, 1]
    seed  <- designPts_mat[i, 2]
    
    treatment_ls <- SimulateData(beta_mat = beta_mat,
                                 Aclusters_df = Aclusters_df,
                                 delta_num = delta,
                                 seed_int = seed)
    betas_df <- treatment_ls$simBetaVals_df
    
    
    for(j in 1:nrow(paramsGrid_mat)){
      
      ###  Calculate Method Output  ###
      adjPval   <- paramsGrid_mat[j, 1]
      mLassoRad <- paramsGrid_mat[j, 2]
      minDmrSep <- paramsGrid_mat[j, 3]
      
      res_ls <- RunProbeLasso(betaVals_mat = betas_df,
                              cpgLocation_df = CPGs_df,
                              adjPvalProbe_num = adjPval,
                              meanLassoRadius_int = mLassoRad,
                              minDmrSep_int = minDmrSep)
      
      ###  Define NULL Data  ###
      if(is.null(res_ls[[1]])){
        res_ls[[1]] <- data.frame(
          dmr.chr      = NA,
          dmr.start    = NA_integer_,
          dmr.end      = NA_integer_,
          seqnames     = NA,
          start        = NA_integer_,
          end          = NA_integer_,
          width        = NA_integer_,
          strand       = NA,
          dmrNo        = NA_integer_,
          dmrP         = NA_real_,
          dmrpRank     = NA_integer_,
          dmrChrom     = NA_character_,
          dmrStart     = NA_integer_,
          dmrEnd       = NA_integer_,
          dmrSize      = NA_integer_,
          dmrCoreStart = NA_integer_,
          dmrCoreEnd   = NA_integer_,
          dmrCoreSize  = NA_integer_,
          ensemblID    = NA_character_,
          geneSymbol   = NA_character_,
          betaAv_Normal= NA_real_,
          betaAV_Tumor = NA_real_,
          dmr.pval     = NA_real_,
          dmr.n.cpgs   = NA_integer_
        )
      }
      
      ###  Save Results  ###
      file_char <- paste0(
        resultsDir, "ProbeLassoResults_delta", delta,
        "_seed", seed,
        "_adjPvalProbe", adjPval,
        "_meanLassoRd", mLassoRad,
        "_minDmrSep", minDmrSep,
        ".RDS"
      )
      
      if(verbose){
        message("Saving results to file ", file_char, "\n")
      }
      
      saveRDS(res_ls, file = file_char)
      
    } # END for(j)
    
  } # END for(i)  
  
}
ProcessProbeLassoResults <- function(resultsDir,
                                     beta_mat,
                                     AclustCPG_df,
                                     verbose = TRUE){
  # browser()
  
  ###  Vector of Appropriate File Names  ###
  files_char <- list.files(path = resultsDir, pattern = "^ProbeLassoResults_.*\\.RDS$")
  files_char <- gsub(pattern = ".RDS", replacement = "", files_char)
  splitFileNames_ls <- strsplit(files_char, split = "_")
  
  
  ###  Initialize Design Points  ###
  ExtractParamLevels <- function(splitNames_ls, param){
    
    posIdx <- grep(param, splitNames_ls[[1]])
    point_char <- sapply(splitNames_ls, `[`, posIdx)
    point_num <- as.numeric(
      gsub(pattern = param, replacement = "", point_char)
    )
    
    sort(unique(point_num))
    
  }
  
  deltas_num <- ExtractParamLevels(splitFileNames_ls, param = "delta")
  seeds_int <- ExtractParamLevels(splitFileNames_ls, param = "seed")
  pVals_num <- ExtractParamLevels(splitFileNames_ls, param = "adjPvalProbe")
  aveLassoRad_int <- ExtractParamLevels(splitFileNames_ls, param = "meanLassoRd")
  minDmrSep_int <- ExtractParamLevels(splitFileNames_ls, param = "minDmrSep")
  
  # This ensures that the delta values stay together, not the seed values.
  designPts_mat <- expand.grid(seeds_int, deltas_num)
  paramsGrid_mat <- expand.grid(pVals_num, aveLassoRad_int, minDmrSep_int)
  
  
  ###  Simulate Gold Standard Outer Loop  ###
  out_ls <- vector(mode = "list", length = nrow(designPts_mat))
  for(i in 1:nrow(designPts_mat)){
    # browser()
    
    ###  Generate Data Set  ###
    seed  <- designPts_mat[i, 1]
    delta <- designPts_mat[i, 2]
    
    treatment_ls <- SimulateData(beta_mat = beta_mat,
                                 Aclusters_df = AclustCPG_df,
                                 delta_num = delta,
                                 seed_int = seed)
    trueClusters_df <- treatment_ls$simAclusters_df
    
    
    ###  Inner Results Comparison  ###
    innerOut_ls <- vector(mode = "list", length = nrow(paramsGrid_mat))
    for(j in 1:nrow(paramsGrid_mat)){
      # browser()
      
      ###  Calculate Method Output  ###
      adjPval   <- paramsGrid_mat[j, 1]
      mLassoRad <- paramsGrid_mat[j, 2]
      minDmrSep <- paramsGrid_mat[j, 3]
      
      
      ###  Load Results  ###
      resFileName_char <- paste0(resultsDir, "ProbeLassoResults",
                                 "_delta", delta, "_seed", seed,
                                 "_adjPvalProbe", adjPval,
                                 "_meanLassoRd", mLassoRad,
                                 "_minDmrSep", minDmrSep,
                                 ".RDS")
      res_ls <- readRDS(resFileName_char)
      
      
      ###  Clean and Summarize Results  ###
      all_df <- CleanResults(dmrResults_ls = res_ls,
                             Aclusters_df = trueClusters_df)
      innerOut_df <- SummarizeResults(cleanDMR_df = all_df,
                                      time_num = res_ls[[2]])
      innerMeta_df <- data.frame(method = "ProbeLasso",
                                 delta = delta,
                                 seed = seed,
                                 adjPval = adjPval,
                                 mLassoRad = mLassoRad,
                                 minDmrSep = minDmrSep,
                                 stringsAsFactors = FALSE)
      innerOut_ls[[j]] <- cbind(innerMeta_df, innerOut_df)
      
    } # END for(j)
    
    
    ###  Bind Inner Results List  ###
    out_ls[[i]] <- do.call(rbind, innerOut_ls)
    
    if(verbose){
      message("Completed summary for delta = ", delta,
              " and seed = ", seed, ".")
    }
    
  } # END for(i)
  
  ###  Bind Outer Results List and Return  ###
  out_df <- do.call(rbind, out_ls)
  # This would return the results to the console, need to save it 
  out_df
}



CleanResults <- function(dmrResults_ls, Aclusters_df) {
  
  ranges_df <- dmrResults_ls[[1]]
  # Remove rows with all NAs
  ranges_df <- ranges_df[rowSums(is.na(ranges_df)) != ncol(ranges_df), ]
  elapsedtime <- dmrResults_ls[[2]]
  
  clusters_df <- unique(
    Aclusters_df[, c("Clusternumber",
                     "chromosome",
                     "start_position",
                     "end_position",
                     "actual")]
  )
  
  if(nrow(ranges_df) > 0){
    
    ###  Create GRanges  ###
    # query = significant DMRs; need to limit to min.cpgs > 4 and pval < 0.05
    signifRanges_df <-
      ranges_df[ranges_df$dmr.n.cpgs > 4 & ranges_df$dmr.pval < 0.05, ]
    query_GR <- GRanges(seqnames = signifRanges_df$dmr.chr,
                        ranges = IRanges(signifRanges_df$dmr.start,
                                         signifRanges_df$dmr.end))
    
    # subject = Aclusters
    subject_GR <- GRanges(seqnames = clusters_df$chromosome,
                          ranges = IRanges(clusters_df$start_position,
                                           clusters_df$end_position))
    
    # subject-query overlap
    overlap_df <- as.data.frame(
      findOverlaps(query_GR, subject_GR, type = "any", select = "all")
    )
    overlap_df$dmr.order <- as.numeric(overlap_df$queryHits)
    overlap_df$aclust.order <- as.numeric(overlap_df$subjectHits)
    
    
    ###  Merge Results  ###
    # merge with dmr info
    signifRanges_df$dmr.row <- 1:nrow(signifRanges_df)
    signifRanges_df$predicted <- "positive"
    overlapDMRs_df <- merge(x = overlap_df, y = signifRanges_df,
                            by.x = "dmr.order", by.y = "dmr.row", all = TRUE)
    
    # merge with aclust info
    clusters_df$aclust.row <- 1:nrow(clusters_df)
    
    # merge all
    all_df <- merge(x = overlapDMRs_df, y = clusters_df,
                    by.x = "aclust.order", by.y = "aclust.row", all = TRUE)
    
  } else {
    
    all_df <- clusters_df
    all_df$predicted <- "negative"
    
  }
  
  ###  Add Status Column  ###
  all_df$predicted[is.na(all_df$predicted)] <- "negative"
  all_df$actual[is.na(all_df$actual)] <- "negative"
  
  # True Positive
  all_df$status[
    all_df$actual == "positive" & all_df$predicted == "positive"
  ] <- "TP"
  
  # False Negative
  all_df$status[
    all_df$actual == "positive" & all_df$predicted == "negative"
  ] <- "FN"
  
  # False Positive
  all_df$status[
    all_df$actual == "negative" & all_df$predicted == "positive"
  ] <- "FP"
  
  # True Negative
  all_df$status[
    all_df$actual == "negative" & all_df$predicted == "negative"
  ] <- "TN"
  
  ###  Return  ###
  all_df
  
}

