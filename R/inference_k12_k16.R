
# Sparse k-means k12 ----------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/FEDRA_varsxclust_k12.Rdata")

# Togliere variabili con "iowa" e "bart" nel nome
# data_features_all_2=data_k12[,-grep("iowa|bart",colnames(data_features_all))]

# Escludi soggetti EXC
excl=which(data_k12$GRP=="EXC")
data_k12=data_k12[-excl,]

# Escludi IDK e GRP
features_k12=subset(data_k12, select = -c(IDK, GRP))

# Select K and sparsity parameter

# install.packages("sparcl")  # uncomment if you haven't installed sparcl yet
library(sparcl)

# 1) 
# Standardizza   
X <- features_k12
X <- scale(X)

# 3) Define range of cluster numbers and sparsity bounds
Ks      <- 2:5                          # try K = 2,3,4,5
p       <- ncol(X)
maxBound <- sqrt(p)                    # maximum wbound = sqrt(#features)
wbounds <- seq(from = 1.1,               # lower > 1 to induce sparsity
               to   = maxBound, 
               length.out = 10)         

# 4) Permutation test to select optimal (K, wbound) via GAP statistic
set.seed(123)
all_results <- lapply(Ks, function(k) {
  perm.out <- KMeansSparseCluster.permute(
    x       = X,
    K       = k,
    wbounds = wbounds,
    nperms  = 100,
    silent  = FALSE
  )
  # perm.out$gaps is a matrix length(wbounds) × 1 (since K fixed), or
  # sometimes 1 × length(wbounds). Coerce properly:
  gaps_vec <- as.vector(perm.out$gaps)
  data.frame(
    K      = k,
    wbound = wbounds,
    Gap    = gaps_vec
  )
})

# combine into one data.frame
gap_df <- do.call(rbind, all_results)

#-Inspect the head of the table ---
head(gap_df)

# Plot Gap vs wbound, one curve per K ---
ggplot(gap_df, aes(x = wbound, y = Gap, color = factor(K), shape = factor(K))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x     = "Sparsity parameter",
    y     = "Gap statistic",
    color = "K",
    shape = "K"
  ) +
  theme_minimal()

K= 2
wbound= 1.75

# Silhouette index
source("Utils_FEDRA.R")

res_sil <- evaluate_sparse_kmeans_silhouette(
  X       = X,
  K_vals  = Ks,
  wb_vals = wbounds,
  nstart  = 20,
  silent  = TRUE
)


res_sil
library(ggplot2)
ggplot(res_sil, aes(x = wbound,
                    y = avg_sil,
                    color = factor(K),
                    group = factor(K))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Dark2", name = "Number of\nClusters (K)") +
  labs(x = "Sparsity Bound (wbound)",
       y = "Average Silhouette",
       title = "Silhouette vs. wbound for Different K") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title     = element_text(face = "bold", hjust = 0.5)
  )

# Bad results with silhouette index

# 6) Run K-means with selected (K, wbound)

skm_=KMeansSparseCluster(X, K=K, wbounds = wbound, silent =
                                FALSE)
# 7) Inspect results

weights <- skm_[[1]]$ws
weights=weights/sum(weights) # normalize weights to sum to 1
features <- colnames(X)

library(ggplot2)

df_w <- data.frame(
  feature = features,
  weight  = weights
)

weights_plot=ggplot(df_w, aes(x = reorder(feature, weight), y = weight)) +
  geom_col() +
  coord_flip() +
  labs(
    #title = "Sparse k-Means Feature Weights",
    x     = "Feature",
    y     = "Weight"
  ) +
  theme_minimal()


weights_plot
# 8) Inspect cluster assignments

write.table(weights,file="weights_sparse_kmeans_k12.txt")

res_skm_=data.frame(
  data_k12,
  clust=skm_[[1]]$Cs
)

table(res_skm_$GRP, res_skm_$clust)

table(res_skm_$clust)/101*100

library(dplyr)
selected_feats <- df_w %>%
  filter(weight > 0) %>%
  pull(feature)

# Compute conditional means by cluster
state_means <- res_skm_ %>%     # replace with the name of your 2nd data.frame
  group_by(clust) %>%
  summarise(across(all_of(selected_feats),
                   ~ mean(.x, na.rm = TRUE),
                   .names = "mean_{.col}"))

print(state_means)

write.table(res_skm_,file="results_sparse_kmeans_k12.txt")

# Sparse k-means k12 WITHOUT IDK 13----------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/FEDRA_varsxclust_k12.Rdata")

# Togliere variabili con "iowa" e "bart" nel nome
# data_features_all_2=data_k12[,-grep("iowa|bart",colnames(data_features_all))]

# Escludi soggetti EXC
excl=which(data_k12$GRP=="EXC")
data_k12=data_k12[-excl,]

# Escludi IDK 13

idk13_row=which(data_k12$IDK=="#013")

data_k12=data_k12[-idk13_row,]

# Escludi IDK e GRP
features_k12=subset(data_k12, select = -c(IDK, GRP))

# Select K and sparsity parameter

# install.packages("sparcl")  # uncomment if you haven't installed sparcl yet
library(sparcl)

# 1) 
# Standardizza   
X <- features_k12
X <- scale(X)

# 3) Define range of cluster numbers and sparsity bounds
Ks      <- 2:5                          # try K = 2,3,4,5
p       <- ncol(X)
maxBound <- sqrt(p)                    # maximum wbound = sqrt(#features)
wbounds <- seq(from = 1.1,               # lower > 1 to induce sparsity
               to   = maxBound, 
               length.out = 10)         

# 4) Permutation test to select optimal (K, wbound) via GAP statistic
set.seed(123)
all_results <- lapply(Ks, function(k) {
  perm.out <- KMeansSparseCluster.permute(
    x       = X,
    K       = k,
    wbounds = wbounds,
    nperms  = 100,
    silent  = FALSE
  )
  # perm.out$gaps is a matrix length(wbounds) × 1 (since K fixed), or
  # sometimes 1 × length(wbounds). Coerce properly:
  gaps_vec <- as.vector(perm.out$gaps)
  data.frame(
    K      = k,
    wbound = wbounds,
    Gap    = gaps_vec
  )
})

# combine into one data.frame
gap_df <- do.call(rbind, all_results)

#-Inspect the head of the table ---
head(gap_df)

library(ggplot2)
# Plot Gap vs wbound, one curve per K ---
ggplot(gap_df, aes(x = wbound, y = Gap, color = factor(K), shape = factor(K))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x     = "Sparsity parameter",
    y     = "Gap statistic",
    color = "K",
    shape = "K"
  ) +
  theme_minimal()

K= 2
wbound= 1.75

# Silhouette index
source("Utils_FEDRA.R")

# res_sil <- evaluate_sparse_kmeans_silhouette(
#   X       = X,
#   K_vals  = Ks,
#   wb_vals = wbounds,
#   nstart  = 20,
#   silent  = TRUE
# )
# 
# 
# res_sil
# library(ggplot2)
# ggplot(res_sil, aes(x = wbound,
#                     y = avg_sil,
#                     color = factor(K),
#                     group = factor(K))) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   scale_color_brewer(palette = "Dark2", name = "Number of\nClusters (K)") +
#   labs(x = "Sparsity Bound (wbound)",
#        y = "Average Silhouette",
#        title = "Silhouette vs. wbound for Different K") +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "right",
#     plot.title     = element_text(face = "bold", hjust = 0.5)
#   )

# 6) Run K-means with selected (K, wbound)

skm_=KMeansSparseCluster(X, K=K, wbounds = wbound, silent =
                           FALSE)
# 7) Inspect results

weights <- skm_[[1]]$ws
weights=weights/sum(weights) # normalize weights to sum to 1
features <- colnames(X)

library(ggplot2)

df_w <- data.frame(
  feature = features,
  weight  = weights
)

weights_plot=ggplot(df_w, aes(x = reorder(feature, weight), y = weight)) +
  geom_col() +
  coord_flip() +
  labs(
    #title = "Sparse k-Means Feature Weights",
    x     = "Feature",
    y     = "Weight"
  ) +
  theme_minimal()


weights_plot
# 8) Inspect cluster assignments

write.table(weights,file="weights_sparse_kmeans_k12_without_IDK13.txt")

res_skm_=data.frame(
  data_k12,
  clust=skm_[[1]]$Cs
)

table(res_skm_$GRP, res_skm_$clust)

table(res_skm_$clust)/101*100

library(dplyr)
selected_feats <- df_w %>%
  filter(weight > 0) %>%
  pull(feature)

# Compute conditional means by cluster
state_means <- res_skm_ %>%     # replace with the name of your 2nd data.frame
  group_by(clust) %>%
  summarise(across(all_of(selected_feats),
                   ~ mean(.x, na.rm = TRUE),
                   .names = "mean_{.col}"))

print(state_means)

write.table(state_means,file="state_means_sparse_kmeans_k12_without_IDK13.txt")

write.table(res_skm_,file="results_sparse_kmeans_k12_without_IDK13.txt")



# Comparison with and without IDK 13 --------------------------------------

results_sparse_kmeans_k12 <- read.csv("D:/CNR/OneDrive - CNR/Brancati-Cortese/results_sparse_kmeans_k12.txt", sep="")

results_sparse_kmeans_k12_without_IDK13 <- read.csv("D:/CNR/OneDrive - CNR/Brancati-Cortese/data and results/results_excluding_IDK13/results_sparse_kmeans_k12_without_IDK13.txt", sep="")
# trim IDK 13 to compare

idk13=which(results_sparse_kmeans_k12$IDK=="#013")

results_sparse_kmeans_k12_compare=results_sparse_kmeans_k12[-idk13,]

table(results_sparse_kmeans_k12_compare$clust,results_sparse_kmeans_k12_without_IDK13$clust)

idx_change=which(results_sparse_kmeans_k12_compare$clust!=results_sparse_kmeans_k12_without_IDK13$clust)

data_idx_change=cbind(results_sparse_kmeans_k12_compare[idx_change,c("IDK","GRP","clust")],results_sparse_kmeans_k12_without_IDK13[idx_change,c("clust")])

colnames(data_idx_change)=c("IDK","GRP","clust_k12","clust_k12_without_IDKD13")
data_idx_change


# Sparse k-means k16 ------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/FEDRA_varsxclust_k12.Rdata")
# Togliere variabili con "iowa" e "bart" nel nome
data_k12=data_k12[,-grep("iowa|bart",colnames(data_k12))]

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_features_all.Rdata")
data_iowa_bart=data_features_all[,c("IDK","beta","rho","Arew","Apun",
                                    "K","betaF","betaP")]

data_k16=merge(data_k12, data_iowa_bart, by="IDK")

features_k16=subset(data_k16, select = -c(IDK, GRP))

# Select K and sparsity parameter

# install.packages("sparcl")  # uncomment if you haven't installed sparcl yet
library(sparcl)

# 1) 
# Standardizza   
X <- features_k16
X <- scale(X)

# 3) Define range of cluster numbers and sparsity bounds
Ks      <- 2:5                          # try K = 2,3,4,5
p       <- ncol(X)
maxBound <- sqrt(p)                    # maximum wbound = sqrt(#features)
wbounds <- seq(from = 1.1,               # lower > 1 to induce sparsity
               to   = maxBound, 
               length.out = 10)         

# 4) Permutation test to select optimal (K, wbound) via GAP statistic
set.seed(123)
all_results <- lapply(Ks, function(k) {
  perm.out <- KMeansSparseCluster.permute(
    x       = X,
    K       = k,
    wbounds = wbounds,
    nperms  = 100,
    silent  = FALSE
  )
  # perm.out$gaps is a matrix length(wbounds) × 1 (since K fixed), or
  # sometimes 1 × length(wbounds). Coerce properly:
  gaps_vec <- as.vector(perm.out$gaps)
  data.frame(
    K      = k,
    wbound = wbounds,
    Gap    = gaps_vec
  )
})

# combine into one data.frame
gap_df <- do.call(rbind, all_results)

#-Inspect the head of the table ---
head(gap_df)

# Plot Gap vs wbound, one curve per K ---
ggplot(gap_df, aes(x = wbound, y = Gap, color = factor(K), shape = factor(K))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x     = "Sparsity parameter",
    y     = "Gap statistic",
    color = "K",
    shape = "K"
  ) +
  theme_minimal()

K= 2
wbound= 1.75

# Silhouette index
source("Utils_FEDRA.R")
res_sil_k16 <- evaluate_sparse_kmeans_silhouette(
  X       = X,
  K_vals  = Ks,
  wb_vals = wbounds,
  nstart  = 20,
  silent  = TRUE
)

res_sil_k16
library(ggplot2)
ggplot(res_sil_k16, aes(x = wbound,
                    y = avg_sil,
                    color = factor(K),
                    group = factor(K))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Dark2", name = "Number of\nClusters (K)") +
  labs(x = "Sparsity Bound (wbound)",
       y = "Average Silhouette",
       title = "Silhouette vs. wbound for Different K") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title     = element_text(face = "bold", hjust = 0.5)
  )

# 6) Run K-means with selected (K, wbound)

skm_=KMeansSparseCluster(X, K=K, wbounds = wbound, silent =
                           FALSE)
# 7) Inspect results

weights <- skm_[[1]]$ws
weights=weights/sum(weights) # normalize weights to sum to 1
features <- colnames(X)


write.table(weights,file="weights_sparse_kmeans_k16.txt")

library(ggplot2)

df_w <- data.frame(
  feature = features,
  weight  = weights
)

weights_plot=ggplot(df_w, aes(x = reorder(feature, weight), y = weight)) +
  geom_col() +
  coord_flip() +
  labs(
    #title = "Sparse k-Means Feature Weights",
    x     = "Feature",
    y     = "Weight"
  ) +
  theme_minimal()


weights_plot
# 8) Inspect cluster assignments

res_skm_=data.frame(
  data_k16,
  clust=skm_[[1]]$Cs
)

table(res_skm_$GRP, res_skm_$clust)

table(res_skm_$clust)/101*100

library(dplyr)
selected_feats <- df_w %>%
  filter(weight > 0) %>%
  pull(feature)

# Compute conditional means by cluster
state_means <- res_skm_ %>%     # replace with the name of your 2nd data.frame
  group_by(clust) %>%
  summarise(across(all_of(selected_feats),
                   ~ mean(.x, na.rm = TRUE),
                   .names = "mean_{.col}"))

print(state_means)

write.table(res_skm_,file="results_sparse_kmeans_k16.txt")

# t clust -----------------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/FEDRA_varsxclust_k7.RData")
exc=data_k7[which(data_k7$GRP=="EXC"),1]
data_k7_cleaned=data_k7[-which(data_k7$IDK%in%exc),]
features_k7=data_k7_cleaned[,-(1:2)]

library(tclust)
# Dataset ristretto
n <- dim(features_k7)[1]
ctl_small=ctlcurves(features_k7, k = 2:5, alpha = (0:10) / n,
                    restr.fact = 12,parallel=T)

x11()
plot(
  ctl_small, 
  main = "CTL Curves",restr.fact =12)

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_features_all.Rdata")
# Togliere variabili con "iowa" e "bart" nel nome
data_features_all_2=data_features_all[,-grep("iowa|bart",colnames(data_features_all))]
features_k29=subset(data_features_all_2, select = -c(IDK, GRP))

# Full dataset
n <- dim(features_k29)[1]
ctl_large=ctlcurves(features_k29, k = 2:5, alpha = (0:10) / n,
                    restr.fact = 12,parallel=T)

x11()
plot(
  ctl_large, 
  main = "CTL Curves",restr.fact =12)

# Risultato interessante per k29, teniamo questo

tc_large=tclust(
  features_k29, k = 4, alpha = 0.02, nstart = 50,restr.fact = 12,opt='MIXT'
)

# Inspect cluster assignments

outl=which(tc_large$cluster==0)

outl

# Max probability
hist(apply(tc_large$posterior[-outl,],1,max))
# Fuzziness not informative

res_tcl_K4_alpha2=data.frame(
  data_features_all_2,
  clust=tc_large$cluster
)

# State cond means
state_means_all <- res_tcl_K4_alpha2 %>%
  group_by(clust) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

print(state_means_all)

table(res_tcl_K4_alpha2$GRP, res_tcl_K4_alpha2$clust)



# COSA k 12 --------------------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/FEDRA_varsxclust_k12.Rdata")

# Escludi soggetti EXC
excl=which(data_k12$GRP=="EXC")
data_k12=data_k12[-excl,]

# Escludi IDK e GRP
features_k12=subset(data_k12, select = -c(IDK, GRP))

X <- features_k12
X <- scale(X)

source("Utils_FEDRA.R")
library(cluster)
library(factoextra)

gap_COSA=COSA_gap(X,
                  K_grid    = 2:5,
                  zeta_grid = seq(0.01, 1.5, .1),
                  tol       = 1e-8,
                  n_outer   = 15,
                  alpha     = .1,
                  verbose   = FALSE,
                  n_cores   = 9,
                  B         = 50,
                  n_init=5)


gap_res=gap_COSA$gap_stats

library(ggplot2)
ggplot(gap_res, aes(x = zeta0, y = GAP, color = factor(K), group = K)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    #title = "GAP statistic vs zeta0",
    x = "Sparsity parameter",
    y = "Gap statistic",
    color = "K"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 13))


K=2
zeta0=0.2


cosak12=COSA2(Y=X,K=K,zeta0=zeta0,
              n_init=10,
              tol=1e-10,n_outer=20,alpha=.1,verbose=T)


round(cosak12$weights,2)


weights <- data.frame(cosak12$weights)
colnames(weights) <- colnames(X)
weights

write.table(weights,file="weights_cosak12.txt")

library(ggplot2)

weights_df <- as.data.frame(weights)
weights_df$cluster <- paste0("Cluster ", seq_len(nrow(weights_df)))

# reshape to long format
weights_df <- as.data.frame(weights)
weights_df$cluster <- paste0("Cluster ", seq_len(nrow(weights_df)))

# reshape to long format
library(dplyr)
library(tidyr)

weights_long <- pivot_longer(weights_df, 
                             cols = -cluster, 
                             names_to = "feature", 
                             values_to = "weight")


# plot with same y-axis scale across all features
ggplot(weights_long, aes(x = cluster, y = weight, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ feature, scales = "fixed") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold")) +
  labs(
    x = " ", 
       y = "Weight"
       #title = "Cluster-wise Feature Weights"
       )

res_cosak12=data.frame(
  data_k12,
  clust=cosak12$cluster
)

table(res_cosak12$GRP, res_cosak12$clust)

table(res_cosak12$clust)/101*100

library(dplyr)
th=.05
feature_totals <- colSums(weights)
selected_features <- names(feature_totals[feature_totals > th])

cluster_means <- res_cosak12 %>%
  group_by(clust) %>%
  summarise(across(all_of(selected_features), mean, na.rm = TRUE))

cluster_means

write.table(res_cosak12,file="results_res_cosak12.txt")

# COSA k 16 ---------------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/FEDRA_varsxclust_k12.Rdata")
# Togliere variabili con "iowa" e "bart" nel nome
data_k12=data_k12[,-grep("iowa|bart",colnames(data_k12))]

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_features_all.Rdata")
data_iowa_bart=data_features_all[,c("IDK","beta","rho","Arew","Apun",
                                    "K","betaF","betaP")]

data_k16=merge(data_k12, data_iowa_bart, by="IDK")

features_k16=subset(data_k16, select = -c(IDK, GRP))

# 1) 
# Standardizza   
X <- features_k16
X <- scale(X)

source("Utils_FEDRA.R")
library(cluster)
library(factoextra)

gap_COSA=COSA_gap(X,
                  K_grid    = 2:5,
                  zeta_grid = seq(0.01, 1.5, .1),
                  tol       = 1e-8,
                  n_outer   = 15,
                  alpha     = .1,
                  verbose   = FALSE,
                  n_cores   = 9,
                  B         = 50,
                  n_init=5)


gap_res=gap_COSA$gap_stats

library(ggplot2)
ggplot(gap_res, aes(x = zeta0, y = GAP, color = factor(K), group = K)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    #title = "GAP statistic vs zeta0",
    x = "Sparsity parameter",
    y = "Gap statistic",
    color = "K"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 13))


K=2
zeta0=0.6


cosak16=COSA2(Y=X,K=K,zeta0=zeta0,
              n_init=10,
              tol=1e-10,n_outer=20,alpha=.1,verbose=T)


round(cosak16$weights,2)


weights <- data.frame(cosak16$weights)
colnames(weights) <- colnames(X)
weights

write.table(weights,file="weights_cosak16.txt")

library(ggplot2)

weights_df <- as.data.frame(weights)
weights_df$cluster <- paste0("Cluster ", seq_len(nrow(weights_df)))

# reshape to long format
weights_df <- as.data.frame(weights)
weights_df$cluster <- paste0("Cluster ", seq_len(nrow(weights_df)))

# reshape to long format
weights_long <- pivot_longer(weights_df, 
                             cols = -cluster, 
                             names_to = "feature", 
                             values_to = "weight")

# plot with same y-axis scale across all features
ggplot(weights_long, aes(x = cluster, y = weight, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ feature, scales = "fixed") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold")) +
  labs(x = " ", y = "Weight"
       # ,
       # title = "Cluster-wise Feature Weights"
       )

res_cosak16=data.frame(
  data_k16,
  clust=cosak16$cluster
)

table(res_cosak16$GRP, res_cosak16$clust)

table(res_cosak16$clust)/101*100

library(dplyr)
th=.05
feature_totals <- colSums(weights)
selected_features <- names(feature_totals[feature_totals > th])

cluster_means <- res_cosak16 %>%
  group_by(clust) %>%
  summarise(across(all_of(selected_features), mean, na.rm = TRUE))

cluster_means

write.table(res_cosak16,file="results_res_cosak16.txt")
