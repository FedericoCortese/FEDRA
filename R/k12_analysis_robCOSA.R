load("~/Library/CloudStorage/OneDrive-CNR/Brancati-Cortese/data and results/FEDRA_varsxclust_k12.RData")

# TROVARE ALTERNATIVA A GAP E SILHOUETTE

excl=which(data_k12$GRP=="EXC")
data_k12=data_k12[-excl,]

# Escludi IDK e GRP
features_k12=subset(data_k12, select = -c(IDK, GRP))

# Select K and sparsity parameter

library(sparcl)

# 1) 
# Standardizza   
X <- features_k12
X <- scale(X)

source("Utils_sparse_robust_2.R")

zeta0=seq(0.05,.2,by=.025)
K=2:4
qt=c(0.99)

sil_results <- compute_silhouette_grid(
  X = X,
  zeta0 = zeta0,
  K = K,
  qt = qt,
  ncores=8
)

sil_results

gap_res <- compute_gap_grid(
  X = X,
  zeta0 = zeta0,
  qt = qt,
  K.max = 4,
  B = 20,
  ncores = 8
)

library(ggplot2)

ggplot(gap_res, aes(x = zeta0, y = gap, color = factor(K), group = K)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "zeta0",
    y = "GAP statistic",
    color = "K"
  ) +
  theme_minimal()

temp=robust_COSA(Y=X,
                 zeta0=zeta0[1],
                 K=K[1],
                 tol     = NULL,
                 n_init  = 5,
                 n_outer = 20,
                 n_inner = 10,
                 alpha   = 0.1,
                 verbose = T,
                 knn     = 10,
                 M       = NULL,
                 qt      = qt[1],
                 c       = NULL,
                 mif     = NULL,
                 outlier = TRUE,
                 truth   = NULL,
                 ncores=NULL)

#cluster
head(temp$s)

head(temp$W)

round(temp$W,2)
table(temp$s)
table(temp$s,data_k12$GRP)

