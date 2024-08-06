ConvertCPGList <- function(cpgs_ls, methylval_df){
  
  cluster_ls <- lapply(cpgs_ls, function(item){
    as.data.frame(methylval_df[item, ])
  })
  clusterRowLabel <- lapply(cpgs_ls, function(item){
    as.data.frame(rownames(methylval_df[item, ]))
  })
  
  cluster_tab <- rbindlist(cluster_ls,
                           idcol = "cluster",
                           use.names = TRUE)
  clusterLabel_tab <- rbindlist(clusterRowLabel,
                                idcol = "cluster",
                                use.names = TRUE)
  colnames(clusterLabel_tab)[2] <- "probeID"
  
  cbind(clusterLabel_tab[, 2], cluster_tab)
}