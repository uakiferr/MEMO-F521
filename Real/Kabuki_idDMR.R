library(idDMR)
library(ChAMP)
library(ChAMPdata)
library(minfi)
library("ChIPpeakAnno")
library("limma")
library(DMRcate)
library("DMRcatedata")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(dplyr)
library(plyr)
library(PRROC)


################ load all data needed for the analysis #######
data.dir <-"/Users/Asus/Documents/GSEs/Kabuki/"

Kabuki_dmr_df <- read.csv(paste0(data.dir, "Kabuki_DMR_to_find_true.csv"), header = TRUE)

Kabuki_dmr_df$actual <- c(rep("positive", length(nrow(Kabuki_dmr_df))))

rownames(Kabuki_dmr_df) <- Kabuki_dmr_df$cpg # making sure the row are the cgnames 

beta_kabuki  <- read.csv(paste0(data.dir, "beta_kabuki.csv"),
                      row.names = 1, header = TRUE) ### Beta_mat filtered. in A-Clust-results
CPGs_df <- read.csv(paste0(data.dir,"cpgLocation_df.csv"),
                    row.names = 1, header = TRUE) ## based on the DMRcompare, see Aclust-results. 

resultsDir <- "/Users/Asus/Documents/GSEs/Kabuki/idDMR/"

### Make sure that row names are the cpgnames 
rownames(Kabuki_dmr_df) 
rownames(CPGs_df) 

### Fonction to extract ranges from iDDMR 
id_extractRanges <- function(dmroutput,
                             genome = c("hg19", "hg38", "mm10")) {
  genome <- match.arg(genome)
  if (!is(dmroutput, "DMResults")) {
    stop("Error: dmroutput is not a DMResults object. Please create one with aaDMR().")
  }
  coords <- extractCoords(dmroutput@coord)
  coords <- cbind(
    coords,
    dmroutput@no.cpgs,
    dmroutput@min_smoothed_fdr,
    dmroutput@Stouffer,
    dmroutput@Fisher,
    dmroutput@maxdiff,
    dmroutput@meandiff
  )
  coords$chromStart <- as.integer(as.character(coords$chromStart))
  coords$chromEnd <- as.integer(as.character(coords$chromEnd))
  ranges <-
    makeGRangesFromDataFrame(coords, keep.extra.columns = TRUE)
  eh <- ExperimentHub()
  switch(genome,
         hg19 = {
           grt <- eh[["EH3132"]]
         },
         hg38 = {
           grt <- eh[["EH3134"]]
         },
         mm10 = {
           grt <- eh[["EH3136"]]
         }
  )
  genesidx <- as.data.frame(findOverlaps(ranges, grt))
  genesover <-
    tapply(genesidx$subjectHits, genesidx$queryHits, function(x) {
      grt$symbol[x]
    })
  op.A <- sapply(genesover, function(l) {
    paste(l, collapse = ", ")
  })
  name.A <- names(genesover)
  m.A <- as.numeric(name.A)
  M <- length(ranges)
  overlapping.genes <- rep(NA_character_, M)
  overlapping.genes[m.A] <- op.A
  ranges$overlapping.genes <- overlapping.genes
  colnames(values(ranges)) <-
    sub("dmroutput@", "", colnames(values(ranges)))
  ranges
}

### Merging DMRs with CpG 
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

### The change FDR function of the idDMR which tell you to change the FDR 
changeFDR <- function (annot, FDR)
{
  if(!is(annot, "CpGannotated")){
    stop("Error: annot is not a CpGsiteAnnotated object. Please create one with cpgsite.annotate()")
  }
  if(FDR <=0 | FDR >=1){
    stop("Error: please enter an appropriate FDR value, 0 < FDR < 1.")
  }
  annot@ranges$is.sig <- annot@ranges$ind.fdr < FDR
  cat(paste0("Threshold is now set at FDR=", FDR, ", resulting in ",
             sum(annot@ranges$is.sig), " significantly differential CpGs."))
  annot
  
}

#### Give a CPGAnnotate object for the aadmr function of the package. 
### A bit modify for the data. 
cpgsite.annotate <-function(datatype = c("array", "sequencing"),
                            object,
                            what = c("Beta", "M"),
                            arraytype = c("EPIC", "450K"),
                            analysis.type = "differential",
                            design,
                            contrasts = FALSE,
                            cont.matrix = NULL,
                            fdr = 0.05,
                            coef,
                            cpgLocation_df= cpgLocation_df,## in order to have the location
                            ...) {
  analysis.type <- match.arg(analysis.type)
  what <- match.arg(what)
  arraytype <- match.arg(arraytype)
  if (datatype == "array") {
    stopifnot(class(object)[1] %in% c("matrix", "GenomicRatioSet"))
    if (is(object, "matrix")) {
      if (arraytype == "450K") {
        grset <- makeGenomicRatioSetFromMatrix(
          mat = object,
          array = "IlluminaHumanMethylation450k",
          annotation = "ilmn12.hg19",
          what = what
        )
      }
      if (arraytype == "EPIC") {
        grset <- makeGenomicRatioSetFromMatrix(
          mat = object,
          array = "IlluminaHumanMethylationEPIC",
          annotation = "ilm10b4.hg19",
          what = what
        )
      }
    } else {
      grset <- object
    }
    object <- getM(grset)
    switch(analysis.type,
           differential = {
             stopifnot(is.matrix(design))
             if (!contrasts) {
               stopifnot(colnames(design)[1] == "(Intercept)")
             } else {
               stopifnot(!is.null(cont.matrix))
             }
             fit <- lmFit(object, design, ...)
             if (contrasts) {
               stopifnot(coef %in% colnames(cont.matrix))
               fit <- contrasts.fit(fit, cont.matrix)
             }
             fit <- eBayes(fit)
             tt <- topTable(fit, coef = coef, number = nrow(object))
             nsig <- sum(tt$adj.P.Val < fdr)
             if (nsig == 0) {
               message(
                 "Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in aaDMR() to return DMRs, but be warned there is an increased risk of Type I errors."
               )
             }
             if (nsig > 0 & nsig <= 100) {
               message(
                 paste(
                   "Your contrast returned",
                   nsig,
                   "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."
                 )
               )
             }
             if (nsig > 100) {
               message(
                 paste(
                   "Your contrast returned",
                   nsig,
                   "individually significant probes. We recommend the default setting of pcutoff in aaDMR()."
                 )
               )
             }
             betafit <- lmFit(ilogit2(object), design, ...)
             if (contrasts) {
               betafit <- contrasts.fit(betafit, cont.matrix)
             }
             betafit <- eBayes(betafit)
             betatt <-
               topTable(betafit, coef = coef, number = nrow(object))
             m <- match(rownames(tt), rownames(betatt))
             tt$diff <- betatt$logFC[m]
             m <- match(rownames(object), rownames(tt))
             tt <- tt[m,]
             common_rows <-intersect(rownames(object), rownames(cpgLocation_df))
             anno <- cpgLocation_df[common_rows, ] ## this was change in order for it to work
             stat <- tt$t
             annotated <-
               GRanges(
                 as.character(anno$chr),
                 anno$MAPINFO,
                 stat = stat,
                 diff = tt$diff,
                 ind.fdr = tt$adj.P.Val,
                 is.sig = tt$adj.P.Val < fdr
               )
             names(annotated) <- rownames(tt)
           })
    annotated <- sort(annotated)
    return(new("CpGsiteAnnotated", ranges = annotated))
  }
  if (datatype == "sequencing") {
    stop("Sequencing mode is not functional at the moment.")
  } else {
    message("Error: datatype must be one of 'array' or 'sequencing'")
  }
}

### Standardize the output for later comparison and processing of the values. 
StandardizeOutput <- function(methodOut_df,
                              method = c("DMRcate",
                                         "ProbeLasso",
                                         "Bumphunter",
                                         "idDMR"),
                              cpgLocation_df = CPGs_df,
                              dmr.sig.threshold = 0.05,
                              min.cpgs = 3){
  
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
           methodOut_df$dmr.chr   <- ethodOut_df$seqnames
           methodOut_df$dmr.start <- methodOut_df$start
           methodOut_df$dmr.end   <- methodOut_df$end
           
         },
         idDMR = {
           methodOut_df$dmr.pval  <- methodOut_df$Stouffer
           methodOut_df$dmr.chr   <- methodOut_df$seqnames
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

######## Creating a function to samples in groups ##########

sampling <- function(df, patient_size , control_size , seed) {
  
  # Indices for patients and controls
  patient_indices <- 1:11
  control_indices <- 12:22
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Sample indices for patients and controls
  sampled_patient_indices <- sample(patient_indices, patient_size, replace = FALSE)
  sampled_control_indices <- sample(control_indices, control_size, replace = FALSE)
  
  # Select columns based on sampled indices
  sampled_df <- df[, c(sampled_patient_indices, sampled_control_indices)]
  
  return(sampled_df)
}

### Function to run aadmr and then, collect the standardize data. 
RunidDMR <- function(betaVals_mat,
                     labels_fct,
                     cpgLocation_df, G_int, FDR_anno,
                     dmr.sig.threshold = 0.05,
                     min.cpgs = 3, genome = "hg19"){
  
  
  ###  Calculate iddmr Results  ###
  ptm <- proc.time()
  design_mat <- model.matrix(~labels_fct)
  
  # This takes 53.22585 sec
  myannotation <- cpgsite.annotate("array",
                                   as.matrix(betaVals_mat),
                                   what = "Beta",
                                   arraytype = c("EPIC", "450K"),
                                   analysis.type = c("differential"),
                                   design_mat,
                                   contrasts = FALSE,
                                   cont.matrix = NULL,
                                   fdr = FDR_anno, 
                                   coef = 2, 
                                   cpgLocation_df = cpgLocation_df)
  
  # This takes 8.167 sec: 13% of computing time.
  
  # Parallel no supported on Windows
  aadmr_out <- tryCatch(aaDMR(myannotation, g = G_int,  min.cpgs = 2, pcutoff = "fdr"),
                        
                        error = function(e1){ NULL }
                        
  )
  
  elapsedtime <- proc.time() - ptm
  
  
  ###  Extract Results  ###
  
  # If any predicted cluster is identified from DMRcate
  if(!is.null(aadmr_out)){
    
    # Requires DMRcatedata::dmrcatedata
    aadmrOut_df <- data.frame(id_extractRanges(aadmr_out, genome = "hg19"))
    
    results_df <- StandardizeOutput(
      methodOut_df = aadmrOut_df,
      method = "idDMR",
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

### This write the results of the runiddmr 
WriteidDMRResults <- function(beta_mat,
                              CPGs_df,
                              sizes_num = c(3, 5, 7, 10),
                              seeds_int = c(50, 100, 210, 330, 450),
                              G_int = c(200, 250, 500, 750, 1000),
                              resultsDir = "idDMR_compare/",
                              FDR_anno = c(0.001, 0.01, 0.05, 0.1),
                              verbose = TRUE){
  
  
  
  dir.create(resultsDir, showWarnings = FALSE, recursive = TRUE)
  
  ###  Data Simulation Outer Loop  ###
  designPts_mat <- expand.grid(sizes_num, seeds_int)
  paramsGrid_mat <- expand.grid(G_int, FDR_anno)
  
  for(i in 1:nrow(designPts_mat)){
    
    ###  Generate Data Set  ###
    size <- designPts_mat[i, 1]
    seed  <- designPts_mat[i, 2]
    
    sample_beta <- sampling(beta_mat, size , size , seed)
    row.names(sample_beta) <- rownames(beta_mat)
    
    for(j in 1:nrow(paramsGrid_mat)){
      
      ###  Calculate Method Output  ###
      G <- paramsGrid_mat[j, 1]
      FDR  <- paramsGrid_mat[j, 2]
      
      suppressMessages(
        res_ls <- RunidDMR(sample_beta,
                           labels_fct = factor(c(rep("Kabuki", size),
                                                 rep("Control", size))),
                           CPGs_df, G , FDR,
                           dmr.sig.threshold = 0.05,
                           min.cpgs = 3, genome = "hg19")
      )
      
      ###  Define NULL Data  ###
      if(is.null(res_ls[[1]])){
        res_ls[[1]] <- data.frame(
          dmr.chr    = NA_character_,
          dmr.start  = NA_integer_,
          dmr.end    = NA_integer_,
          seqnames   = NA,
          start      = NA_integer_,
          end        = NA_integer_,
          width      = NA_integer_,
          strand     = NA,
          no.cpgs    = NA_integer_,
          min_smoothed_fdr = NA_real_,
          Stouffer   = NA_real_,
          Fisher   = NA_real_,
          maxdiff  = NA_real_,
          meandiff = NA_real_,
          overlapping.genes = NA_character_,
          dmr.pval   = NA_real_,
          dmr.n.cpgs = NA_integer_
          
        )
      }
      
      ###  Save Results  ###
      file_char <- paste0(
        resultsDir, "idDMRResults_samplesize", size, "_seed", seed,
        "_G", G ,"_FDRcpg",FDR , ".RDS"
      )
      
      if(verbose){
        message("Saving results to file ", file_char, "\n")
      }
      
      saveRDS(res_ls, file = file_char)
      
    } # END for(j)
    
    
  } # END for(i)
  
}



### Make sure that row names are the cpgnames 


WriteidDMRResults(beta_kabuki,
                  CPGs_df,
                  sizes_num = c(3, 5, 7, 10),
                  seeds_int = c(50, 100, 210, 330, 450),
                  G_int =  c(200, 250, 500, 750, 1000),
                  FDR_anno = c(0.001, 0.01, 0.05, 0.1),
                  resultsDir = resultsDir,
                  verbose = TRUE)

###################### to process the data, clean and summarize results ###########
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
    signifRanges_df <-ranges_df[ranges_df$dmr.n.cpgs > 2 & ranges_df$dmr.pval < 0.05 & abs(ranges_df$maxdiff)> 0.1,  ]
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
SummarizeResults <- function(cleanDMR_df, time_num){
  # browser()
  
  ###  Table Power Results  ###
  # Frequency count of each status for power
  #   calculation
  
  powerSummary_tbl <- table(
    factor(cleanDMR_df$status,
           levels = c("FN", "FP", "TN", "TP"))
  )
  powerSummary_df <- data.frame(
    matrix(as.numeric(powerSummary_tbl), ncol = 4)
  )
  colnames(powerSummary_df) <- names(powerSummary_tbl)
  
  
  ###  Add Power Calculation Summary ###
  powerSummary_df$power <-
    powerSummary_df$TP / (powerSummary_df$TP + powerSummary_df$FN)
  powerSummary_df$time <- time_num
  powerSummary_df$nPower <- powerSummary_df$TP + powerSummary_df$FN
  powerSummary_df <-
    powerSummary_df[, c("time", "FN", "FP", "TN", "TP", "power", "nPower")]
  
  ###  Table Precision Results  ###
  # Frequency count of each status - based DMRs, for precision
  #   calculation
  if(!is.null(cleanDMR_df$dmr.order)){
    
    statusDMR_df <- cleanDMR_df[, c("dmr.order", "status")]
    
    precisSummary_tbl <- table(
      factor(statusDMR_df$status,
             levels = c("FN", "FP", "TN", "TP"))
    )
    precisSummary_df <- data.frame(
      matrix(as.numeric(precisSummary_tbl), ncol = 4)
    )
    colnames(precisSummary_df) <- paste0(names(precisSummary_tbl), "precis")
    precisSummary_df$FNprecis <- precisSummary_df$TNprecis <- NULL
    
    
    ###  Add Precision Calculation Summary ###
    precisSummary_df$precision <-
      precisSummary_df$TP / (precisSummary_df$TP + precisSummary_df$FP)
    precisSummary_df$nPrecis <- precisSummary_df$TP + precisSummary_df$FP
    
    
    powerSummary_df <- cbind(powerSummary_df, precisSummary_df)
    
  } else {
    
    powerSummary_df$FPprecis  <- NA_integer_
    powerSummary_df$TPprecis  <- NA_integer_
    powerSummary_df$precision <- NA_real_
    powerSummary_df$nPrecis   <- NA_integer_
    
  }
  
  ###  Area under Precision-Recall Curve  ###
  if(!is.null(cleanDMR_df$dmr.pval)){
    
    x_df <- cleanDMR_df[, c("aclust.order", "dmr.pval", "actual")]
    x_df$status <- ifelse(x_df$actual == "positive", 1, 0)
    x_df$dmr.pval[is.na(x_df$dmr.pval)] <- 1
    
    # PR curve
    prCurve_ls <- pr.curve(scores.class0 = 1 - x_df$dmr.pval,
                           weights.class0 = x_df$status)
    powerSummary_df$AuPR <- prCurve_ls$auc.integral
    rm(x_df, prCurve_ls)
    
  } else {
    powerSummary_df$AuPR <- NA_real_
  }
  
  # F1-score
  CalcF1 <- function(tp, fp, fn){ (2 * tp) / (2 * tp + fp + fn) }
  powerSummary_df$F1  <- CalcF1(tp = powerSummary_df$TPprecis,
                                fp = powerSummary_df$FPprecis,
                                fn = powerSummary_df$FN)
  
  
  ###  Summary of numCPGs  ###
  # I've tried to fit the Poisson and Gamma distributions to the number of CPGs,
  #   and a Gamma to the log of the number of CPGs. Nothing fits well. The
  #   counts are minima-inflated. We will report the three quartiles.
  nCPGsSummary_num <- summary(cleanDMR_df$dmr.n.cpgs)
  nCPGs_df <- data.frame(
    matrix(nCPGsSummary_num[c(2, 3, 5)], nrow = 1),
    stringsAsFactors = FALSE
  )
  colnames(nCPGs_df) <- paste0("nCPG_", c("q1", "med", "q3"))
  
  powerSummary_df <- cbind(powerSummary_df, nCPGs_df)
  
  
  ###  Return  ###
  powerSummary_df
  
}

ProcessidDMRResults <- function(resultsDir,
                                Kabuki_dmr_df,
                                verbose = TRUE){
  # browser()
  
  ###  Vector of Appropriate File Names  ###
  files_char <- list.files(path = resultsDir, pattern = "^idDMRResults_.*\\.RDS$")
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
  G_num <- ExtractParamLevels(splitFileNames_ls, param = "G")
  FDR_int <- ExtractParamLevels(splitFileNames_ls, param = "FDRcpg")
  
  # This ensures that the size values stay together, not the seed values.
  designPts_mat <- expand.grid(seeds_int, sizes_num)
  paramsGrid_mat <- expand.grid(G_num , FDR_int)
  
  
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
      
      ###  Calculate Method Output  ###
      G <- paramsGrid_mat[j, 1]
      FDR  <- paramsGrid_mat[j, 2]
      
      
      ###  Load Results  ###
      resFileName_char <- paste0(
        resultsDir, "idDMRResults_samplesize", size, "_seed", seed,
        "_G", G ,"_FDRcpg",FDR , ".RDS"
      )
      res_ls <- readRDS(resFileName_char)
      
      
      ###  Clean and Summarize Results  ###
      all_df <- CleanResults(dmrResults_ls = res_ls,
                             Kabuki_dmr_df = Kabuki_dmr_df)
      innerOut_df <- SummarizeResults(cleanDMR_df = all_df,
                                      time_num = res_ls[[2]])
      innerMeta_df <- data.frame(method = "idDMR",
                                 size = size,
                                 seed = seed,
                                 G = G,
                                 FDR_anno = FDR,
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
  # This would return the results to the console, but we should probably save
  #   the results data frame instead?
  out_df
  
}

idDMR_kabuki <- ProcessidDMRResults(resultsDir= resultsDir,
                                  Kabuki_dmr_df,
                                  verbose = TRUE)

write.csv(idDMR_kabuki , file = paste0(resultsDir,"clean_DMR_idDMR_Kabuki.csv"))
