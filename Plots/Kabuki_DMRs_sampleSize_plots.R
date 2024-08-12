###### Resulting table from Probelasso, idDMR and DMRcate#####

######## Functions ########

CleanResults_PB <- function(dmrResults_ls, Kabuki_dmr_df) {
  
  ranges_df <- dmrResults_ls[[1]]
  # Remove rows with all NAs
  ranges_df <- ranges_df[rowSums(is.na(ranges_df)) != ncol(ranges_df), ]
  elapsedtime <- dmrResults_ls[[2]]
  
  clusters_df <- unique(
    Kabuki_dmr_df[, c("Clusternumber",
                      "chromosome",
                      "start_position",
                      "end_position", 
                      "actual")]
  )
  
  if(nrow(ranges_df) > 0){
    
    ###  Create GRanges  ###
    # query = significant DMRs; need to limit to min.cpgs > 2 and pval < 0.05
    signifRanges_df <-ranges_df[ranges_df$dmr.n.cpgs > 2 & ranges_df$dmr.pval < 0.05 & abs((ranges_df$betaAv_Control)-(ranges_df$betaAv_Kabuki)) > 0.1,  ]
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
  
  
  ###  Return  ###
  all_df
  
}


ProcessProbeLassoResults <- function(resultsDir,
                                     beta_kabuki,
                                     Kabuki_dmr_df,
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
  
  sizes_num <- ExtractParamLevels(splitFileNames_ls, param = "samplesize")
  seeds_int <- ExtractParamLevels(splitFileNames_ls, param = "seed")
  pVals_num <- ExtractParamLevels(splitFileNames_ls, param = "adjPvalProbe")
  aveLassoRad_int <- ExtractParamLevels(splitFileNames_ls, param = "meanLassoRd")
  minDmrSep_int <- ExtractParamLevels(splitFileNames_ls, param = "minDmrSep")
  
  # This ensures that the size values stay together, not the seed values.
  designPts_mat <- expand.grid(seeds_int, sizes_num)
  paramsGrid_mat <- expand.grid(pVals_num, aveLassoRad_int, minDmrSep_int)
  
  
  ###  Simulate Gold Standard Outer Loop  ###
  out_ls <- vector(mode = "list", length = nrow(designPts_mat))
  for(i in 1:nrow(designPts_mat)){
    # browser()
    
    ###  Generate Data Set  ###
    seed  <- designPts_mat[i, 1]
    size <- designPts_mat[i, 2]
    
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
                                 "_samplesize", size, "_seed", seed,
                                 "_adjPvalProbe", adjPval,
                                 "_meanLassoRd", mLassoRad,
                                 "_minDmrSep", minDmrSep,
                                 ".RDS")
      res_ls <- readRDS(resFileName_char)
      
      
      ###  Clean and Summarize Results  ###
      all_df <- CleanResults_PB(dmrResults_ls = res_ls,
                                Kabuki_dmr_df = Kabuki_dmr_df)
      innerOut_df <- SummarizeResults(cleanDMR_df = all_df,
                                      time_num = res_ls[[2]])
      innerMeta_df <- data.frame(method = "ProbeLasso",
                                 size = size,
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
      message("Completed summary for size = ", size,
              " and seed = ", seed, ".")
    }
    
  } # END for(i)
  
  ###  Bind Outer Results List and Return  ###
  out_df <- do.call(rbind, out_ls)
  # This would return the results to the console, need to save it 
  out_df
}

CleanResults <- function(dmrResults_ls, Kabuki_dmr_df) {
  
  ranges_df <- dmrResults_ls[[1]]
  # Remove rows with all NAs
  ranges_df <- ranges_df[rowSums(is.na(ranges_df)) != ncol(ranges_df), ]
  elapsedtime <- dmrResults_ls[[2]]
  
  clusters_df <- unique(
    Kabuki_dmr_df[, c("Clusternumber",
                      "chromosome",
                      "start_position",
                      "end_position", 
                      "actual")]
  )
  
  if(nrow(ranges_df) > 0){
    
    ###  Create GRanges  ###
    # query = significant DMRs; need to limit to min.cpgs > 2 and pval < 0.05
    signifRanges_df <-ranges_df[ranges_df$dmr.n.cpgs > 2 & ranges_df$dmr.pval < 0.05 & abs(ranges_df$meandiff)> 0.1,  ]
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
  
  
  ###  Return  ###
  all_df
  
}


######## Functions for the plots ###########
BuildOverlaps <- function(BestResultsDir,
                          size = c(3,5,7,10),
                          seed = c(50, 100, 210, 330, 450),
                          CPGs_df = cpgLocation_df,
                          min.cpgs = 3) {
  ### List Results Files  ###
  fileNames_char <- list.files(BestResultsDir)
  targetNames_char <- paste0("samplesize", size, "_seed", seed)
  correctFiles_idx <- grep(targetNames_char, fileNames_char)
  correctNames_char <- fileNames_char[correctFiles_idx]
  names(correctNames_char) <- sapply(
    strsplit(correctNames_char, "_"),
    function(x){
      gsub(pattern = "Results", replacement = "", x[1])
    }
  )
  
  ### Load and Clean Appropriate Files  ###
  allRes_ls <- lapply(correctNames_char, function(x){
    
    res_ls <- readRDS(paste0(BestResultsDir, x))
    
    if(grepl("ProbeLasso", x)){
      
      ranges_df <- res_ls[[1]]
      results_df<- ranges_df[ranges_df$dmr.n.cpgs > 2 & ranges_df$dmr.pval < 0.05 & abs((ranges_df$betaAv_Control)-(ranges_df$betaAv_Kabuki)) > 0.1,  ]
    } else {
      ranges_df <- res_ls[[1]]
      res_ls[[1]] <-ranges_df[ranges_df$dmr.n.cpgs > 2 & ranges_df$dmr.pval < 0.05 & abs(ranges_df$maxdiff)> 0.1, ]
    }
    
    results_df <- res_ls[[1]]
    keepRows <- rowSums(is.na(results_df)) != ncol(results_df)
    results_df <- results_df[keepRows, ]
    
    
    if(nrow(results_df) > 0){
      
      results_IRanges <- IRanges(
        start = results_df$dmr.start,
        end = results_df$dmr.end
      )
      GRanges(
        seqnames = results_df$dmr.chr,
        ranges = results_IRanges
      )
      
    }
    
  })
  
  ###  Return  ###
  attr(allRes_ls, "size") <- size
  attr(allRes_ls, "repl")  <- which(c(50, 100, 210, 330, 450) %in% seed)
  allRes_ls
  
}



PlotOverlaps <- function(BestResultsDir,
                         figFileName,
                         device = pdf,
                         plotTitle = "default",
                         sizes_num = c(3, 5, 7, 10),
                         seeds_int = c(50, 100, 210, 330, 450),
                         CPGs_df = cpgLocation_df,
                         min.cpgs = 3,
                         totalDMR = 4569,
                         ...) {
  
  # Group by seeds
  design_mat <- expand.grid(seeds_int, sizes_num)
  cpg_df <- CPGs_df
  minCPGs <- min.cpgs
  
  overlapsByMethods_ls <- lapply(1:nrow(design_mat), function(i) {
    
    overlaps_ls <- BuildOverlaps(
      BestResultsDir = BestResultsDir,
      size = design_mat[i, 2], seed = design_mat[i, 1],
      CPGs_df = cpg_df,
      min.cpgs = minCPGs
    )
    
    overlapAttr <- attributes(overlaps_ls)
    
    null_idx <- sapply(overlaps_ls, is.null)
    overlaps_ls <- overlaps_ls[!null_idx]
    
    attr(overlaps_ls, "samplesize") <- overlapAttr$size
    attr(overlaps_ls, "repl")  <- overlapAttr$repl
    overlaps_ls
    
  })
  
  CreateHue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  # Loop through each combination of size and seed
  for (i in 1:nrow(design_mat)) {
    size <- design_mat[i, 2]
    seed <- design_mat[i, 1]
    
    # Construct file name for each plot
    currentFigFileName <- file.path(figFileName, paste0("plot_samplesize", size, "_seed", seed, ".pdf"))
    
    device(file = currentFigFileName, ...)
    
    overlaps_ls <- overlapsByMethods_ls[[i]]
    
    size_num <- attr(overlaps_ls, "size")
    repl_int  <- attr(overlaps_ls, "repl")
    
    currentPlotTitle <- if (plotTitle == "default") {
      paste0("Venn Diagram for size = ", size_num, ", rep = ", repl_int)
    } else {
      plotTitle
    }
    
    makeVennDiagram(
      overlaps_ls, NameOfPeaks = names(overlaps_ls), totalTest = totalDMR,
      by = "region", fill = CreateHue(length(overlaps_ls)),
      main = currentPlotTitle
    )
    
    dev.off()
  }
}
######### Build a PR curve list for each function best results #############
BuildPRcurve <- function(BestResultsDir,
                         size = c(3, 5, 7, 10),
                         seed = c(50, 100, 210, 330, 450),
                         beta_mat = beta_kabuki,
                         Kabuki_dmr_df,
                         CPGs_df = cpgLocation_df,
                         min.cpgs = 3){
  
  ### List Results Files  ###
  fileNames_char <- list.files(BestResultsDir)
  targetNames_char <- paste0("samplesize", size, "_seed", seed)
  correctFiles_idx <- grep(targetNames_char, fileNames_char)
  correctNames_char <- fileNames_char[correctFiles_idx]
  names(correctNames_char) <- sapply(
    strsplit(correctNames_char, "_"),
    function(x){
      gsub(pattern = "Results", replacement = "", x[1])
    }
  )
  
  ### Load and Clean Appropriate Files  ###
  allRes_ls <- lapply(correctNames_char, function(x){
    
    res_ls <- readRDS(paste0(BestResultsDir, x))
    
    if(grepl("ProbeLasso", x)){
      
      CleanResults_PB(dmrResults_ls = res_ls,
                      Kabuki_dmr_df)
      
    } else {CleanResults(dmrResults_ls = res_ls,
                         Kabuki_dmr_df)}
  })
  
  ###  Build PR Curves  ###
  allPRs_ls <- lapply(allRes_ls, function(cleanDMR_df){
    
    if(!is.null(cleanDMR_df$dmr.n.cpgs)){
      
      x_df <- cleanDMR_df[, c("aclust.order", "dmr.pval", "actual")]
      x_df$status <- ifelse(x_df$actual == "positive", 1, 0)
      x_df$dmr.pval[is.na(x_df$dmr.pval)] <- 1
      
      # PR curve
      pr.curve(scores.class0 = 1 - x_df$dmr.pval,
               weights.class0 = x_df$status,
               curve = TRUE)
      
    } else {
      NULL
    }
    
  })
  
  ###  Return  ###
  attr(allPRs_ls, "size") <- size
  attr(allPRs_ls, "repl")  <- which(c(50, 100, 210, 330, 450) %in% seed)
  allPRs_ls
  
}


PlotPRCurve <- function(prCurves_ls,
                        plotTitle = "default",
                        new = TRUE,
                        lineWidth = 1,
                        colours = NULL){
  
  # Extract Meta
  size <- attr(prCurves_ls, "samplesize")
  repl  <- attr(prCurves_ls, "repl")
  
  # Colours
  if(is.null(colours)){
    
    CreateHue <- function(n){
      
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
      
    }
    
    colours <- CreateHue(length(prCurves_ls))
    
  }
  
  # Foundation Plot
  if(new){
    
    if(plotTitle == "default"){
      plotTitle <- paste0(
        "Precision-recall curve: size = ", size, ", rep = ", repl
      )
    }
    
    # Main
    plot(
      prCurves_ls[[1]],
      color = colours[1],
      auc.main = FALSE, legend = TRUE,
      main = plotTitle,
      cex.main = 1, lwd = lineWidth
    )
    
    # Legend
    legend(x = 0.75, y = 1,
           legend = names(prCurves_ls),
           col = colours, lty = 1, cex = 0.75, box.lty = 0)
    
  } else {
    
    plot(
      prCurves_ls[[1]],
      color = colours[1],
      add = TRUE, legend = TRUE,
      lwd = lineWidth
    )
    
  }
  
  
  # Subsequent
  for(i in 2:length(prCurves_ls)){
    
    if(!is.null(prCurves_ls[[i]])){
      plot(
        prCurves_ls[[i]],
        color = colours[i],
        add = TRUE, legend = TRUE,
        lwd = lineWidth
      )
    }
    
  }
  
}

######### Plots ######
BestResultsDir <- "path/to/your/files/Kabuki/Best_per_size_F1/"


BuildPRC <- BuildPRcurve(BestResultsDir,
                         size = 10, ##3, 5, 7, 10
                         seed = 450, ##50, 100, 210, 330, 450
                         beta_mat = beta_kabuki,
                         Kabuki_dmr_df,
                         CPGs_df = cpgLocation_df,
                         min.cpgs = 3) 


PlotPRCurve(BuildPRC,
            plotTitle = "size 10, rep 5" ,
            new = TRUE,
            lineWidth = 1,
            colours = NULL)


figFileName <- paste0(BestResultsDir, "Fig/")

PlotOverlaps(BestResultsDir,
             figFileName,
             device = pdf,
             plotTitle = "default",
             size = c( 5, 7, 10),
             seed = c(50, 100, 210, 330, 450),
             CPGs_df = CPGs_df ,
             min.cpgs = 3)



