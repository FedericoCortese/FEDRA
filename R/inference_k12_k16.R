
# Sparse k-means k12 ----------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/FEDRA_varsxclust_k12.Rdata")

# Togliere variabili con "iowa" e "bart" nel nome
# data_features_all_2=data_k12[,-grep("iowa|bart",colnames(data_features_all))]

# Escludi soggetti EXC
excl=which(data_k12$GRP=="EXC")
data_k12=data_k12[-excl,]

# Escludi IDK e GRP
features_k12=subset(data_k12, select = -c(IDK, GRP))

library(sparcl)

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

gap_COSA=COSA_gap(Y=X,
                     zeta_grid=seq(0.01,1.5,.1),
                     K_grid=2:5,
                     tol=NULL,n_outer=15,alpha=.1,verbose=F,n_cores=9)
