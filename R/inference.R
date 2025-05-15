
# Sparse k-means ----------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_features_all.Rdata")

# Togliere variabili con "iowa" e "bart" nel nome
data_features_all_2=data_features_all[,-grep("iowa|bart",colnames(data_features_all))]
features_k29=subset(data_features_all_2, select = -c(IDK, GRP))

library(sparcl)

# Select K and sparsity parameter

# install.packages("sparcl")  # uncomment if you haven't installed sparcl yet
library(sparcl)

# 1) 
#    
X <- features_k29

# 2) Standardize features (important for K-means)
X <- scale(X)

# 3) Define range of cluster numbers and sparsity bounds
Ks      <- 2:5                          # try K = 2,3,4,5,6
p       <- ncol(X)
maxBound <- sqrt(p)                    # maximum wbound = sqrt(#features)
wbounds <- seq(from = 1.1,               # lower > 1 to induce sparsity
               to   = maxBound, 
               length.out = 10)         # 20 candidate values

# 4) Permutation test to select optimal (K, wbound) via GAP statistic
set.seed(123)
all_results <- lapply(Ks, function(k) {
  perm.out <- KMeansSparseCluster.permute(
    x       = X,
    K       = k,
    wbounds = wbounds,
    nperms  = 50,
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
ggplot(gap_df, aes(x = wbound, y = Gap, color = factor(K))) +
  geom_line(size = 1) +
  labs(
    x     = "Sparsity bound (wbound)",
    y     = "Gap statistic",
    color = "K"
  ) +
  theme_minimal()

K= 4
wbound= 3

# 6) Run K-means with selected (K, wbound)

skm_K4_w3=KMeansSparseCluster(X, K=4, wbounds = 3, silent =
                                FALSE)
# 7) Inspect results

weights <- skm_K4_w3[[1]]$ws
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
    title = "Sparse K-Means Feature Weights (K=4, w=3)",
    x     = "Feature",
    y     = "Weight"
  ) +
  theme_minimal()

# 8) Inspect cluster assignments

res_skm_K4_w3=data.frame(
  data_features_all_2,
  clust=skm_K4_w3[[1]]$Cs
)

table(res_skm_K4_w3$GRP, res_skm_K4_w3$clust)


selected_feats <- df_w %>%
  filter(weight > 0) %>%
  pull(feature)

# Compute conditional means by cluster
state_means <- res_skm_K4_w3 %>%     # replace with the name of your 2nd data.frame
  group_by(clust) %>%
  summarise(across(all_of(selected_feats),
                   ~ mean(.x, na.rm = TRUE),
                   .names = "mean_{.col}"))

print(state_means)


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


# robust COSA -------------------------------------------------------------

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_features_all.Rdata")

# Togliere variabili con "iowa" e "bart" nel nome
data_features_all_2=data_features_all[,-grep("iowa|bart",colnames(data_features_all))]
features_k29=subset(data_features_all_2, select = -c(IDK, GRP))

source("Utils_FEDRA.R")

X=scale(features_k29)

gap_robCOSA=COSA_gap(Y=X,
                     zeta_grid=seq(0.01,1.5,.1),
                     K_grid=2:5,
                     tol=NULL,n_outer=15,alpha=.1,verbose=F,n_cores=8,
                     B=100,knn=10,c=2,M=NULL)

gap_res=gap_robCOSA$gap_stats

library(ggplot2)
ggplot(gap_res, aes(x = zeta0, y = GAP, color = factor(K), group = K)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "GAP statistic vs zeta0",
    x = expression(zeta[0]),
    y = "GAP",
    color = "K (number of clusters)"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 13))


best_cosa=gap_res[which.min(gap_res$GAP),]

cosa_K5_zeta21=robust_COSA(Y=X,zeta0=best_cosa$zeta0,K=best_cosa$K,tol=NULL,
                     n_outer=50,alpha=.1,verbose=T,knn=10,c=2,M=NULL)

cosa_K4_zeta21=robust_COSA(Y=X,zeta0=best_cosa$zeta0,K=4,tol=NULL,
                           n_outer=50,alpha=.1,verbose=T,knn=10,c=2,M=NULL)

# Features weights
df_w=data.frame(cosa_K4_zeta21$W)
colnames(df_w)=colnames(X)
df_w

colSums(df_w)

th <- 0.1  # Cambia questo valore come vuoi

# Aggiungi colonna Stato
df_w <- df_w %>% mutate(State = paste0("State_", row_number()))

# Reshape in formato long
df_long <- df_w %>%
  pivot_longer(-State, names_to = "Feature", values_to = "Weight")

# Calcola somma totale per feature
total_weights <- df_long %>%
  group_by(Feature) %>%
  summarise(total_weight = sum(Weight))

# Filtra le feature con somma > soglia
features_to_keep <- total_weights %>%
  filter(total_weight >= th) %>%
  pull(Feature)

# Filtro del dataset
df_filtered <- df_long %>%
  filter(Feature %in% features_to_keep)

# Plot
ggplot(df_filtered, aes(x = Feature, y = Weight, fill = State)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = paste("Feature Weights by State (threshold =", th, ")"),
       x = "Feature", y = "Weight") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")



#
sel_feat=unique(df_filtered$Feature)

res_cosa=data.frame(
  data_features_all_2[,c(sel_feat,"GRP")],
  clust=cosa_K4_zeta21$s,
  outl=I(cosa_K4_zeta21$v==0)
)

state_cond=res_cosa %>%
  select(where(is.numeric), clust) %>%
  group_by(clust) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  print(n = Inf)
state_cond

attach(res_cosa)

prop.table(table(clust[outl==F]))*100
table(GRP[outl==F],clust[outl==F])

# Outliers

res_outl=res_cosa[res_cosa$outl==T,c(features_to_keep,"GRP")]

state_mean_outl=res_outl %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

state_mean_outl
