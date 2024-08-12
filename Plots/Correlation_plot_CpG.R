library(corrplot)


colnames(Kabuki_DMR_to_find) <- Kabuki_DMR_to_find[1,]
for_Correlation <- Kabuki_DMR_to_find[6:24, ]
colnames(for_Correlation) <- colnames(Kabuki_DMR_to_find)
rownames(for_Correlation) <- for_Correlation$cpg

beta_value <- (for_Correlation[,9:30])
beta_value<- matrix(as.numeric(c(beta_value$GSM2562771, 
                                  beta_value$GSM2562772,
                                  beta_value$GSM2562773,
                                  beta_value$GSM2562774)),nrow = 19,ncol = 4 )
rownames(beta_value) <- rownames(for_Correlation) 
corr <- cor(t(beta_value)) 
corrplot(corr, method = "square", type = "full")
