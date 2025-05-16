
load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_features_all.Rdata")

# Togliere variabili con "iowa" e "bart" nel nome
data_features_all_2=data_features_all[,-grep("iowa|bart",colnames(data_features_all))]
features_k29=subset(data_features_all_2, select = -c(IDK, GRP))



source("Utils_FEDRA.R")

X=scale(features_k29)

gap_robCOSA=COSA_gap(Y=X,
                     zeta_grid=seq(0.01,1.5,.1),
                     K_grid=2:5,
                     tol=NULL,n_outer=15,alpha=.1,verbose=F,n_cores=3,
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

library(tidyr)

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