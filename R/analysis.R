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

head(output_bart$BUGSoutput$sims.matrix[,1])
plot(output_bart$BUGSoutput$sims.matrix[,1],type='l')
hist(output_bart$BUGSoutput$sims.matrix[,1])

# Beta
library(ggplot2)
x11()
ggplot(data_features_all, aes(x = GRP, y = beta, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of Beta by Group",
       x = "Group",
       y = "Beta") +
  theme(legend.position = "none")

# Rho
x11()
ggplot(data_features_all, aes(x = GRP, y = rho, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of Rho by Group",
       x = "Group",
       y = "Rho") +
  theme(legend.position = "none")

# Arew
x11()
ggplot(data_features_all, aes(x = GRP, y = Arew, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of Arew by Group",
       x = "Group",
       y = "Arew") +
  theme(legend.position = "none")

# Apun
x11()
ggplot(data_features_all, aes(x = GRP, y = Apun, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of Apun by Group",
       x = "Group",
       y = "Apun") +
  theme(legend.position = "none")

# K
x11()
ggplot(data_features_all, aes(x = GRP, y = K, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of K by Group",
       x = "Group",
       y = "K") +
  theme(legend.position = "none")

# betaF
x11()
ggplot(data_features_all, aes(x = GRP, y = betaF, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of betaF by Group",
       x = "Group",
       y = "betaF") +
  theme(legend.position = "none")

# betaP
x11()
ggplot(data_features_all, aes(x = GRP, y = betaP, fill = GRP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # boxplot without showing outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # optional: individual points
  theme_minimal() +
  labs(title = "Boxplot of betaP by Group",
       x = "Group",
       y = "betaP") +
  theme(legend.position = "none")


# (1) Considera tutte le variabili, se al clustering uniamo una procedura di variable selection 
# eccetto le variabili con 'bart' (1) e 'iowa' (4) e ovviamente aggiungere beta, rho, Arew etc...
# hard clustering con un algoritmo COSA

# (2) Sparse k-means (hard) 

# (3) Fuzzy k-means (soft): una prob. di appartenere a ciascun gruppo a ciascun paziente 
# Da provare cosi' com'e' e magari fare PCA prima oppure combinare (pero non c'e' in letteratura, tipo fuzzy-COSA)

# (4) Misture gaussiane (soft) con (tutte le variabili) o senza selezione delle variabili (usa solo le k7 in questo caso) 
# oppure con trimming (opportuna gestione degli outlier), risultato ulteriore: rilevazione delle anomalie

# N.B. Trimming importante (trimmed-COSA?)




