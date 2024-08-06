######  Library  ##############################################################
data.dir <- "/path/to/Simulation/"


library(Aclust)
library(DMRcate)
library(data.table)
library(DMRcompare)



######  1. Aclust - prepare annotation file  ##################################

# data downloaded from https://rforge.net/IMA/ , under '4 Annotation file'
# This loads 11 lists, from 5-19Mb, and the full annotation matrix (0.5Gb)
load(paste0(data.dir, "fullannotInd.rda"))
rm(EXON1Ind, GENEBODYInd, ISLANDInd, NSHELFInd, NSHOREInd, SSHELFInd,
   SSHOREInd, TSS1500Ind, TSS200Ind, UTR3Ind, UTR5Ind)

# Turn off "stringsAsFactors"
annot <- as.data.frame(fullannot, stringsAsFactors = FALSE)


###  CPG Locations  ###
vars <- c("ILMNID", "CHR", "MAPINFO")
cpg.location <- annot[vars]
names(cpg.location) <- c("ILMNID", "chr", "MAPINFO")
cpg.location$MAPINFO <- as.numeric(cpg.location$MAPINFO)
cpg.location$chr <- paste0("chr", cpg.location$chr)
cpg.location$ILMNID <- as.factor(cpg.location$ILMNID)

cpg.location <- cpg.location[substr(cpg.location$ILMNID, 1, 2) != "rs", ]
anyNA(cpg.location)

# cpgLocation_df <- cpg.location
# devtools::use_data(cpgLocation_df)
saveRDS(cpg.location, paste0(data.dir, "my.cpg.location.RDS"))


###  Cluster by Chromosome  ###
vars <- c("ILMNID", "INFINIUM_DESIGN_TYPE", "CHR", "MAPINFO",
          "UCSC_REFGENE_NAME", "UCSC_REFGENE_GROUP", "UCSC_CPG_ISLANDS_NAME",
          "RELATION_TO_UCSC_CPG_ISLAND")
annots450k <- annot[vars]

# change colname 'MAPINFO' by name 'Coordinate_37'
colnames(annots450k) <- c("IlmnID", "Infinium_Design_Type", "CHR",
                          "Coordinate_37", "UCSC_RefGene_Name",
                          "UCSC_RefGene_Group", "UCSC_CpG_Islands_Name",
                          "Relation_to_UCSC_CpG_Island")

annots450k$Coordinate_37 <- as.numeric(annots450k$Coordinate_37)

# For Aclust, the locus must be sorted according to 'CHR' and
#   'Coordinate_37/MAPINFO' together in ascending order
annot_450k <- annots450k[order(annots450k$CHR, annots450k$Coordinate_37), ]
nrow(annot_450k) # 485577 cpgs
# Remove the cpgs (loci) that have no CHR (chromosome) name
annot_450k <- annot_450k[!is.na(annot_450k$CHR), ]
nrow(annot_450k) # 485512 cpgs

# make a list of annotation files for different chromosomes
chrome_annot_files <- lapply(1:22, function(chrom){
  
  out <- annot_450k[(annot_450k$CHR == chrom), ]
  subset(out, !duplicated(out[, 4]))
  
})

nrow(chrome_annot_files[[1]]) # 46857 cpgs

# Save our results and clean up the workspace.
saveRDS(chrome_annot_files,
         paste0(data.dir, "annotation-450K-Aclust-by-chromosome.rds"))
rm(fullannot, vars, annots450k, annot_450k)



######  2. Beta value files  ##################################################

beta_value1 <- read.csv(paste0(data.dir, "beta_simu.csv"),
                        row.names = 1, header = TRUE)
nrow(beta_value1) # 410752 cpgs



# Keep only the probes for which (0.05 < betaval < 0.95) in all probes. 
betaLogi_mat <- (as.matrix(beta_value1) > 0.05 & as.matrix(beta_value1) < 0.95)
betaKeep_logi <- rowSums(betaLogi_mat) > 0
beta.value.all <- beta_value1[betaKeep_logi, ]
nrow(beta.value.all) # 335039 cpgs


###  Save Beta Matrix  ###
# # Necessary for the data simulation step (2_simulatedata.R) and further steps
betaVals_mat <- beta.value.all

# devtools::use_data(betaVals_mat)

saveRDS(beta.value.all, paste0(data.dir, "beta.value.all.rds"))
write.csv(betaVals_mat,file = "/path/to/Beta_for_simu.csv")
load("beta.value.all.rds")
######  3. Assign cpgs to clusters  ###########################################

# This takes 12-13 minutes:
clusters_ls <- lapply(1:22, function(chrome){
  
  assign.to.clusters(betas = beta.value.all,
                     annot = as.data.table(chrome_annot_files[[chrome]]),
                     dist.thresh = 0.5, bp.merge = 200,
                     dist.type = "spearman", method = "complete")
  
})

results_ls <- unlist(clusters_ls, recursive = FALSE)
length(results_ls) # 287633

# select clusters with 5 cpgs or more
resultsTrim_ls <- results_ls[lengths(results_ls) >= 5]
length(resultsTrim_ls) # 2479

source("ConvertCPGList.R")

a <- Sys.time()
cluster_tab <- ConvertCPGList(cpgs_ls = resultsTrim_ls,
                              methylval_df = beta.value.all)
Sys.time() - a # 2.298284 min

rm(clusters_ls, results_ls, resultsTrim_ls)



######  4. Add ranges for clusters  ###########################################

clusterOrd_tab <- cluster_tab[order(cluster_tab$cluster, cluster_tab$probeID)]

###  merge with annotation to get position info  ###
annot.location <- annot[, c("CHR", "MAPINFO")]
annot.location$probeID <- row.names(annot.location)
annot.location <- annot.location[annot.location$probeID, ]
location_tab <- merge(clusterOrd_tab, annot.location, by = "probeID")
rm(cluster_tab, annot.location, clusterOrd_tab)


###  Start and End of Each Cluster  ###
location_tab$MAPINFO <- as.numeric(as.character(location_tab$MAPINFO))
# Order the location_tab row according to cluster id ('cluster') and 'MAPINFO'
#   together in ascending order
location_tab <-
  location_tab[order(location_tab$cluster, location_tab$MAPINFO), ]

# Start
start_df <- aggregate(MAPINFO ~ as.factor(cluster),
                      data = location_tab, function(x) min(x))
colnames(start_df) <- c("cluster", "start_position")
start_df$cluster <- as.numeric(as.character(start_df$cluster))

# End
end_df <- aggregate(MAPINFO ~ as.factor(cluster),
                    data = location_tab, function(x) max(x))
colnames(end_df) <- c("cluster", "end_position")
end_df$cluster <- as.numeric(as.character(end_df$cluster))


###  Merge with Methylation Values  ###
start_tab <- merge(location_tab, start_df, by = "cluster")
startEnd_tab <- merge(start_tab, end_df, by = "cluster")
rm(location_tab, start_df, end_df, start_tab)


###  Rename and Order Columns  ###
colnames(startEnd_tab)[1] <- "Clusternumber"
colnames(startEnd_tab)[2] <- "cpg"
startEnd_tab$coordinate_37 <- startEnd_tab$MAPINFO
startEnd_tab$chromosome <- paste0("chr", startEnd_tab$CHR)
startEnd_df <- as.data.frame(startEnd_tab)

impVars_df <- startEnd_df[c("Clusternumber", "cpg", "CHR", "MAPINFO",
                            "start_position", "end_position", "coordinate_37",
                            "chromosome")]
otherVars_df <- startEnd_df[setdiff(names(startEnd_df), colnames(impVars_df))]
startEndCPG_df <- cbind(impVars_df, otherVars_df)

# For the RunDMRcate() function to properly detect DMRs, the rownames of this
#   data frame must be equal to the CPGs for this data frame
rownames(startEndCPG_df) <- startEndCPG_df$cpg

# devtools::use_data(startEndCPG_df)
write.csv(startEndCPG_df,file = paste0(data.dir, "A-clust-results.csv"), row.names = FALSE)
rm(startEnd_tab, startEnd_df, impVars_df, otherVars_df)



