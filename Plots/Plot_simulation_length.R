Bumphunter_01 <- BumphunterResults_delta0.1_seed340_pickQ0.95_maxGap750[[1]]
idDMR_01 <- idDMRResults_delta0.1_seed340_G1000_FDRcpg0.05[[1]]
ProbeLasso_01 <- ProbeLassoResults_delta0.1_seed340_adjPvalProbe0.01_meanLassoRd1000_minDmrSep500[[1]]
DMRcate_01 <- DMRcateResults_delta0.1_seed340_lambda1000_C2[[1]]

Bumphunter_01$method <- "Bumphunter"
idDMR_01$method <- "idDMR"
ProbeLasso_01$method <- "ProbeLasso"
DMRcate_01$method <- "DMRcate"

combined_01 <- bind_rows(Bumphunter_01,idDMR_01, ProbeLasso_01, DMRcate_01)
combined_01$method <- factor(combined_01$method ,levels = c("Bumphunter", "idDMR", "ProbeLasso","DMRcate"))

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

plot01 <- ggplot(combined_01, aes(x = method, y = dmr.n.cpgs, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#e74c3c", "#189AB4", "#CF9CD9",  "#2ecc71")) +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9) + 
  labs(y = "Number of CpG site in the DMRs") +
  theme_classic() + theme(axis.title.x = element_blank(), 
                          plot.title = element_blank(), 
                          axis.text.x =element_text(size=14), 
                          legend.position = "none")

Bumphunter_025 <- BumphunterResults_delta0.25_seed340_pickQ0.95_maxGap750[[1]]
idDMR_025 <- idDMRResults_delta0.25_seed340_G1000_FDRcpg0.05[[1]]
ProbeLasso_025 <- ProbeLassoResults_delta0.25_seed340_adjPvalProbe0.01_meanLassoRd1000_minDmrSep500[[1]]
DMRcate_025 <- DMRcateResults_delta0.25_seed340_lambda1000_C2[[1]]

Bumphunter_025$method <- "Bumphunter"
idDMR_025$method <- "idDMR"
ProbeLasso_025$method <- "ProbeLasso"
DMRcate_025$method <- "DMRcate"

combined_025 <- bind_rows(Bumphunter_025,idDMR_025, ProbeLasso_025, DMRcate_025)
combined_025$method <- factor(combined_025$method ,levels = c("Bumphunter", "idDMR", "ProbeLasso","DMRcate"))

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

plot025 <- ggplot(combined_025, aes(x = method, y = dmr.n.cpgs, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#e74c3c", "#189AB4", "#CF9CD9",  "#2ecc71")) +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9) + 
  labs(y = "Number of CpG site in the DMRs") +
  theme_classic()+ theme(axis.title.x = element_blank(), 
                         plot.title = element_blank(), 
                         axis.text.x =element_text(size=14), 
                         legend.position = "none")


Bumphunter_04 <- BumphunterResults_delta0.4_seed340_pickQ0.95_maxGap750[[1]]
idDMR_04 <- idDMRResults_delta0.4_seed340_G1000_FDRcpg0.05[[1]]
ProbeLasso_04 <- ProbeLassoResults_delta0.4_seed340_adjPvalProbe0.01_meanLassoRd1000_minDmrSep500[[1]]
DMRcate_04 <- DMRcateResults_delta0.4_seed340_lambda1000_C2[[1]]

Bumphunter_04$method <- "Bumphunter"
idDMR_04$method <- "idDMR"
ProbeLasso_04$method <- "ProbeLasso"
DMRcate_04$method <- "DMRcate"

combined_04 <- bind_rows(Bumphunter_04,idDMR_04, ProbeLasso_04, DMRcate_04)
combined_04$method <- factor(combined_04$method ,levels = c("Bumphunter", "idDMR", "ProbeLasso","DMRcate"))

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

plot04 <- ggplot(combined_04, aes(x = method, y = dmr.n.cpgs, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#e74c3c", "#189AB4", "#CF9CD9",  "#2ecc71")) +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9) + 
  labs(y = "Number of CpG site in the DMRs") +
  theme_classic() + theme(axis.title.x = element_blank(), 
                          plot.title = element_blank(), 
                          axis.text.x =element_text(size=14), 
                          legend.position = "none")


ggplot(combined_01, aes(x = dmr.chr, y =dmr.pval)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6) +  # Jitter to avoid overplotting
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Optional: Add boxplot for summary
  labs(x = "Chromosome",
       y = "P-value") +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +  # Add horizontal line at 0.05
  theme_classic() +
  facet_wrap(~method)

ggplot(combined_025, aes(x = dmr.chr, y =dmr.pval)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6) +  # Jitter to avoid overplotting
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Optional: Add boxplot for summary
  labs(x = "Chromosome",
       y = "P-value") +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +  # Add horizontal line at 0.05
  theme_classic() +
  facet_wrap(~method)

ggplot(combined_04, aes(x = dmr.chr, y =dmr.pval)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6) +  # Jitter to avoid overplotting
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Optional: Add boxplot for summary
  labs(x = "Chromosome",
       y = "P-value") +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +  # Add horizontal line at 0.05
  theme_classic() +
  facet_wrap(~method)




