# (i) Considera due dataset: ristretto (k7) e completo

# (ii) tclust, small dataset ------------------------------------------------------------------

# Dataset ristretto
#
load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/FEDRA_varsxclust_k7.RData")
exc=data_k7[which(data_k7$GRP=="EXC"),1]
data_k7_cleaned=data_k7[-which(data_k7$IDK%in%exc),]
features_k7=data_k7_cleaned[,-(1:2)]

## Selezione di alpha e K

# Dataset ristretto
n <- dim(features_k7)[1]
ctl_small=ctlcurves(features_k7, k = 2:6, alpha = (0:10) / n,
                    restr.fact = 12,parallel=T)

x11()
plot(
  ctl_small, 
  main = "CTL Curves",restr.fact =12)

# alpha=4%, K=4
tc_small=tclust(
  features_k7, k = 4, alpha = 0.04, nstart = 100,restr.fact = 12,opt='MIXT'
)

# Cluster assignements
tc_small_postprobs=data.frame(tc_small$posterior)
colnames(tc_small_postprobs)=c("p1","p2","p3","p4")
res_tc_small=data.frame(data_k7_cleaned,tc_small_postprobs)
res_tc_small$cluster=apply(res_tc_small[,c("p1","p2","p3","p4")],1,which.max)

# Outliers
tc_small_outliers=which(tc_small$cluster==0)
res_tc_small$cluster[tc_small_outliers]=0

# Cluster Assignements
table(res_tc_small$cluster)


# Cluster characterization
state_cond_medians=matrix(NA,nrow=ncol(features_k7),ncol=5)

for(i in 1:ncol(features_k7)){
  state_cond_medians[i,]=tapply(features_k7[,i],res_tc_small$cluster,median)
}
state_cond_medians=data.frame(state_cond_medians)
row.names(state_cond_medians)=colnames(features_k7)
colnames(state_cond_medians)=0:4

state_cond_medians

# Correspondance with original assignement
table(res_tc_small$cluster,res_tc_small$GRP)

# Representation on the PC space
tc_small.pca <- prcomp(features_k7, scale = TRUE)
res_tc_small$PC1=tc_small.pca$x[,1]
res_tc_small$PC2=tc_small.pca$x[,2]

# Plot
max_prob <- apply(res_tc_small[,c("p1","p2","p3","p4")], 1, max)
cex_values <- 1+(max_prob)^10
x11()
plot(res_tc_small[,c("PC1","PC2")], 
     col = res_tc_small$cluster + 1, 
     pch = 19,
     cex = cex_values,
     xlab = "PC1", 
     ylab = "PC2")

legend("topright", 
       legend = c("Outlier","Cluster 1","Cluster 2","Cluster 3","Cluster 4"), 
       col = 1:5, 
       pch = 19, 
       cex = 0.8)


# Plot interattivo
library(plotly)

# Calcolo della grandezza del punto
max_prob <- apply(res_tc_small[, c("p1", "p2", "p3", "p4")], 1, max)
# Parametri
size_min <- 5
size_max <- 20

# Crea un vettore vuoto
cex_values <- rep(NA, length(max_prob))

# Se max_prob è 0, assegno size 1 (praticamente invisibile)
cex_values[max_prob == 0] <- 1

# Se max_prob > 0, scala tra size_min e size_max
cex_values[max_prob > 0] <- size_min + 
  (max_prob[max_prob > 0] - min(max_prob[max_prob > 0])) / 
  (max(max_prob[max_prob > 0]) - min(max_prob[max_prob > 0])) * (size_max - size_min)

res_tc_small$cex_values=cex_values
res_tc_small$max_prob=max_prob

# Preparazione dei tooltip: un testo con tutte le info per ciascun soggetto
res_tc_small$tooltip <- apply(res_tc_small, 1, function(row) {
  paste0("ID: ", row["IDK"], "\n",
         "Group: ", row["GRP"], "\n",
         "cpt_com: ", row["cpt_com"], "\n",
         "wais_ds: ", row["wais_ds"], "\n",
         "wais_ss: ", row["wais_ss"], "\n",
         "wcst_pers: ", row["wcst_pers"], "\n",
         "bart_adj: ", row["bart_adj"], "\n",
         "iowa_risk_net: ", row["iowa_risk_net"], "\n",
         "esst_neri: ", row["esst_neri"], "\n",
         "Cluster: ", row["cluster"], "\n",
         "p1: ", round(as.numeric(row["p1"]), 3), "\n",
         "p2: ", round(as.numeric(row["p2"]), 3), "\n",
         "p3: ", round(as.numeric(row["p3"]), 3), "\n",
         "p4: ", round(as.numeric(row["p4"]), 3), "\n",
         "cex: ", row["cex_values"], "\n",
         "max_prob: ", row["max_prob"])
})


# Plot interattivo corretto
fig <- plot_ly(data = res_tc_small,
               x = ~PC1, 
               y = ~PC2,
               type = 'scatter',
               mode = 'markers',
               color = ~factor(cluster),
               colors = c("black", "red", "green", "blue", "orange"),
               text = ~tooltip,
               #marker = list(size = ~cex_values),
               hoverinfo = 'text')

fig <- fig %>% layout(title = "PCA plot interattivo",
                      xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"))

fig


# (iii) tclust, complete dataset ------------------------------------------------------------------
# Dataset completo
#
load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_features_all.Rdata")

# Togliere variabili con "iowa" e "bart" nel nome
data_features_all_2=data_features_all[,-grep("iowa|bart",colnames(data_features_all))]
features_k29=subset(data_features_all_2, select = -c(IDK, GRP))

library(tclust)

n <- dim(features_k29)[1]
ctl_large=ctlcurves(features_k29, k = 2:6, alpha = (0:10) / n,
                    restr.fact = 12,parallel=T)

x11()
plot(
  ctl_large, 
  main = "CTL Curves",restr.fact =12)

# alpha=4%, K=4
tc_large=tclust(
  features_k29, k = 4, alpha = 0.04, nstart = 100,restr.fact = 12,opt='MIXT'
)

# Cluster assignements
tc_large_postprobs=data.frame(tc_large$posterior)
colnames(tc_large_postprobs)=c("p1","p2","p3","p4")
res_tc_large=data.frame(data_features_all_2,tc_large_postprobs)
res_tc_large$cluster=apply(res_tc_large[,c("p1","p2","p3","p4")],1,which.max)

summary(res_tc_large[,c("p1","p2","p3","p4")])

# Outliers
tc_large_outliers=which(tc_large$cluster==0)
res_tc_large$cluster[tc_large_outliers]=0

# Cluster Assignements
table(res_tc_large$cluster)


# Cluster characterization
state_cond_medians=matrix(NA,nrow=ncol(features_k29),ncol=5)

for(i in 1:ncol(features_k29)){
  state_cond_medians[i,]=tapply(features_k29[,i],res_tc_large$cluster,median)
}
state_cond_medians=data.frame(state_cond_medians)
row.names(state_cond_medians)=colnames(features_k29)
colnames(state_cond_medians)=0:4

state_cond_medians

# Correspondance with original assignement
table(res_tc_large$cluster,res_tc_large$GRP)

# Representation on the PC space
tc_large.pca <- prcomp(features_k29, scale = TRUE)
res_tc_large$PC1=tc_large.pca$x[,1]
res_tc_large$PC2=tc_large.pca$x[,2]

# Plot
max_prob <- apply(res_tc_large[,c("p1","p2","p3","p4")], 1, max)
cex_values <- 1+(max_prob)^10
x11()
plot(res_tc_large[,c("PC1","PC2")], 
     col = res_tc_large$cluster + 1, 
     pch = 19,
     cex = cex_values,
     xlab = "PC1", 
     ylab = "PC2")

legend("topright", 
       legend = c("Outlier","Cluster 1","Cluster 2","Cluster 3","Cluster 4"), 
       col = 1:5, 
       pch = 19, 
       cex = 0.8)


# Plot interattivo
library(plotly)

# Calcolo della grandezza del punto
max_prob <- apply(res_tc_large[, c("p1", "p2", "p3", "p4")], 1, max)
# Parametri
size_min <- 5
size_max <- 20

# Crea un vettore vuoto
cex_values <- rep(NA, length(max_prob))

# Se max_prob è 0, assegno size 1 (praticamente invisibile)
cex_values[max_prob == 0] <- 1

# Se max_prob > 0, scala tra size_min e size_max
cex_values[max_prob > 0] <- size_min + 
  (max_prob[max_prob > 0] - min(max_prob[max_prob > 0])) / 
  (max(max_prob[max_prob > 0]) - min(max_prob[max_prob > 0])) * (size_max - size_min)

res_tc_large$cex_values=cex_values
res_tc_large$max_prob=max_prob

# Preparazione dei tooltip: un testo con tutte le info per ciascun soggetto
res_tc_large$tooltip <- apply(res_tc_large, 1, function(row) {
  paste0("ID: ", row["IDK"], "\n",
         "Group: ", row["GRP"], "\n",
         "beta: ", row["beta"], "\n",
         "Arew: ", row["Arew"], "\n",
         "Apun: ", row["Apun"], "\n",
         "K: ", row["K"], "\n",
         "betaF: ", row["betaF"], "\n",
         "betaP: ", row["betaP"], "\n",
         "cpt_com: ", row["cpt_com"], "\n",
         "wais_ds: ", row["wais_ds"], "\n",
         "wais_ss: ", row["wais_ss"], "\n",
         "wcst_pers: ", row["wcst_pers"], "\n",
         "esst_neri: ", row["esst_neri"], "\n",
         "Cluster: ", row["cluster"], "\n",
         "p1: ", round(as.numeric(row["p1"]), 3), "\n",
         "p2: ", round(as.numeric(row["p2"]), 3), "\n",
         "p3: ", round(as.numeric(row["p3"]), 3), "\n",
         "p4: ", round(as.numeric(row["p4"]), 3), "\n",
         "cex: ", row["cex_values"], "\n",
         "max_prob: ", row["max_prob"])
})


# Plot interattivo corretto
fig <- plot_ly(data = res_tc_large,
               x = ~PC1, 
               y = ~PC2,
               type = 'scatter',
               mode = 'markers',
               color = ~factor(cluster),
               colors = c("black", "red", "green", "blue", "orange"),
               text = ~tooltip,
               marker = list(size = ~cex_values),
               hoverinfo = 'text')

fig <- fig %>% layout(title = "PCA plot interattivo",
                      xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"))

fig
