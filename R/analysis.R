# load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/FEDRA_data_cleaned.RData")
# load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_bart_posterior.RData")
# load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_iowa_posteriors.RData")

# merge all
data_features_all=merge(data_bart_posterior,data_iowa_posteriors,by="IDK")
data_features_all=subset(data_features_all,select=-subjID)
data_features_all=merge(data_features_all,data_cleaned,by="IDK")
data_features_all$GRP <- as.factor(data_features_all$GRP)

unique(data_features_all$GRP)
head(data_features_all)

# Beta
ggplot(data_features_all, aes(x = GRP, y = beta, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of Beta by Group",
       x = "Group",
       y = "Beta") +
  theme(legend.position = "none")

# Rho
ggplot(data_features_all, aes(x = GRP, y = rho, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of Rho by Group",
       x = "Group",
       y = "Rho") +
  theme(legend.position = "none")

# Arew
ggplot(data_features_all, aes(x = GRP, y = Arew, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of Arew by Group",
       x = "Group",
       y = "Arew") +
  theme(legend.position = "none")

# Apun
ggplot(data_features_all, aes(x = GRP, y = Apun, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of Apun by Group",
       x = "Group",
       y = "Apun") +
  theme(legend.position = "none")

# K
ggplot(data_features_all, aes(x = GRP, y = K, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of K by Group",
       x = "Group",
       y = "K") +
  theme(legend.position = "none")

# betaF
ggplot(data_features_all, aes(x = GRP, y = betaF, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of betaF by Group",
       x = "Group",
       y = "betaF") +
  theme(legend.position = "none")

# betaP
ggplot(data_features_all, aes(x = GRP, y = betaP, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of betaP by Group",
       x = "Group",
       y = "betaP") +
  theme(legend.position = "none")
