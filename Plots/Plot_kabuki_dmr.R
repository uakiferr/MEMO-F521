
library(ggplot2)
library(dplyr) 

DMRcate <- DMRcateResults_samplesize11_seed50_lambda250_C5[[1]]
idDMR <- idDMRResults_samplesize11_seed50_G500_FDRcpg0.001[[1]]
ProbeLasso <- ProbeLassoResults_samplesize11_seed50_adjPvalProbe0.01_meanLassoRd700_minDmrSep200[[1]]

DMRcate$method <- "DMRcate"
idDMR$method <- "idDMR"
ProbeLasso$method <- "ProbeLasso"
combined_df <- bind_rows(DMRcate, idDMR, ProbeLasso)

combined_df$method <- factor(combined_df$method ,levels = c("DMRcate", "idDMR", "ProbeLasso"))

get_box_stats <- function(y, upper_limit = NULL) {
  if (is.null(upper_limit)) {
    upper_limit <- max(y, na.rm = TRUE) * 1.15
  }
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

ggplot(combined_df, aes(x = method, y = dmr.n.cpgs, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#e74c3c", "#2ecc71", "#0099f0")) +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9) + 
  labs(x = "DMR detection methods",
       y = "Number of CpG site in the DMRs") +
  theme_classic()


ggplot(combined_df, aes(x = dmr.chr, y =dmr.pval)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6) +  # Jitter to avoid overplotting
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Optional: Add boxplot for summary
  labs(x = "Chromosome",
       y = "P-value") +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +  # Add horizontal line at 0.05
  theme_classic() +
  facet_wrap(~method)


idDMR_sorted <- idDMR[idDMR$dmr.n.cpgs > 2 & idDMR$dmr.pval < 0.05 & abs(idDMR$meandiff) > 0.1, ]
DMRcate_sorted <- DMRcate[DMRcate$dmr.n.cpgs > 2 & DMRcate$dmr.pval < 0.05 & abs(DMRcate$meandiff)> 0.1,  ]
Probelasso_sorted <-ProbeLasso[ProbeLasso$dmr.n.cpgs > 2 & ProbeLasso$dmr.pval < 0.05 & abs((ProbeLasso$betaAv_Control)-(ProbeLasso$betaAv_Kabuki)) > 0.1,  ]
combined_df <- bind_rows(idDMR_sorted, DMRcate_sorted, Probelasso_sorted)
combined_df$method <- factor(combined_df$method ,levels = c("DMRcate", "idDMR", "ProbeLasso"))

get_box_stats <- function(y, upper_limit = NULL) {
  if (is.null(upper_limit)) {
    upper_limit <- max(y, na.rm = TRUE) * 1.15
  }
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

ggplot(combined_df, aes(x = method, y = dmr.n.cpgs, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#e74c3c", "#2ecc71", "#0099f0")) +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9) + 
  labs(x = "DMR detection methods",
       y = "Number of CpG site in the DMRs") +
  theme_classic()


ggplot(combined_df, aes(x = dmr.chr, y =dmr.pval)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6) +  # Jitter to avoid overplotting
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Optional: Add boxplot for summary
  labs(x = "Chromosome",
       y = "P-value") +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +  # Add horizontal line at 0.05
  theme_classic() +
  facet_wrap(~method)



