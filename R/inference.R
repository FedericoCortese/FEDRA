library(cluster)
library(mclust)
library(tclust)
library(sparcl)

# load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_features_all.Rdata")

load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/FEDRA_varsxclust_k7.RData")
exc=data_k7[which(data_k7$GRP=="EXC"),1]

data_k7_cleaned=data_k7[-which(data_k7$IDK%in%exc),]
features_k7=data_k7_cleaned[,-(1:2)]

# How many groups?
length(unique(data_k7_cleaned$GRP))


# (0) Hyperparameters selection -------------------------------------------

# GAP statistics - tclust
factoextra::fviz_nbclust(features_k7,tclust,method = "gap_stat",
                         alpha = 0.03, nstart = 50,restr.fact = 12,opt='MIXT',k.max=6)

# Silhouette - tclust
factoextra::fviz_nbclust(features_k7,tclust,method = "silhouette",
                         alpha = 0.03, nstart = 50,restr.fact = 12,opt='MIXT',k.max=6)

# Curve CTL per diversi livelli di trimming e numero di cluster
n <- dim(features_k7)[1]
x11()
plot(ctlcurves(features_k7, k = 2:6, alpha = (0:10) / n), main = "CTL Curves",restr.fact =12)

# Sparse k means
# Selection of sparsity parameter

P=dim(features_k7)[2]

km.perm=sparcl::KMeansSparseCluster.permute(features_k7,K=4,wbounds=seq(1.1,sqrt(P),length.out=10))
sparcl::KMeansSparseCluster.permute(features_k7,K=5,wbounds=seq(1.1,sqrt(P),length.out=10))
sparcl::KMeansSparseCluster.permute(features_k7,K=6,wbounds=seq(1.1,sqrt(P),length.out=10))

# It always selects 1.27175

factoextra::fviz_nbclust(features_k7,KMeansSparseCluster,
                         method = "gap_stat",
                         wbounds=km.perm$bestw,
                         k.max=6)

cluster::clusGap(features_k7,KMeansSparseCluster,K.max=6,
                 wbounds=km.perm$bestw)

# (1) Sparse k-means ------------------------------------------------------

library(RSKC)

P=dim(features_k7)[2]


cluster::clusGap(features_k7,RSKC::RSKC,K.max=6,
                 alpha=.03, L1 = sqrt(P), nstart = 100,
                 silent=F, scaling = T, correlation = FALSE)


prova=RSKC(features_k7, ncl=4, alpha=.05, L1 = sqrt(P), nstart = 100,
           silent=TRUE, scaling = T, correlation = FALSE)

# (2) trimmed GMM with 7 essential features ---------------------------------------------

library(tclust)  # Per clustering robusto

x11()
par(mfrow=c(3,3))
for(i in 1:ncol(features_k7)){
  boxplot(features_k7[,i], main=colnames(features_k7)[i], ylab="value", xlab="feature",col=i+1)
}

x11()
par(mfrow=c(3,3))
for(i in 1:ncol(features_k7)){
  hist(features_k7[,i], main=colnames(features_k7)[i], ylab="value", xlab="feature",col=i+1)
}


## 6 groups, 2-4% outliers ##
tc6_3 <- tclust(features_k7, k = 6, alpha = 0.03, nstart = 100,restr.fact = 12,opt='MIXT')

plot(tc6_3)

# Outliers
tc6_3_outliers=which(tc6_3$cluster==0)

# Cluster weights

tc6_3_postprobs=data.frame(tc6_3$posterior)
colnames(tc6_3_postprobs)=c("p1","p2","p3","p4","p5","p6")
res_tc6_3=data.frame(data_k7_cleaned,tc6_3_postprobs)
res_tc6_3$cluster=apply(res_tc6_3[,c("p1","p2","p3","p4","p5","p6")],1,which.max)
res_tc6_3$cluster[tc6_3_outliers]=0

head(round(res_tc6_3[,c("p1","p2","p3","p4","p5","p6")],3))
head(res_tc6_3$cluster)

table(res_tc6_3$GRP,res_tc6_3$cluster)

# Representation on the PC space
tc6_3.pca <- prcomp(features_k7, scale = TRUE)

tc6_3.pca.eig.val<-factoextra::get_eigenvalue(tc6_3.pca)
tc6_3.pca.eig.val

# Add two first PC to the data
res_tc6_3$PC1=tc6_3.pca$x[,1]
res_tc6_3$PC2=tc6_3.pca$x[,2]

# Plot
max_prob <- apply(res_tc6_3[,c("p1","p2","p3","p4","p5","p6")], 1, max)
cex_values <- 1+(max_prob)^10
x11()
plot(res_tc6_3[,c("PC1","PC2")], 
     col = res_tc6_3$cluster + 1, 
     pch = 19,
     cex = cex_values,
     xlab = "PC1", 
     ylab = "PC2", 
     main = "PCA of the features")

legend("topright", 
       legend = c("Outlier","Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6"), 
       col = 1:5, 
       pch = 19, 
       cex = 0.8)


## 4 groups, 2-4% outliers ##
tc4_3 <- tclust(features_k7, k = 4, alpha = 0.03, nstart = 100,restr.fact = 12,opt='MIXT')

plot(tc4_3)

# Outliers
tc4_3_outliers=which(tc4_3$cluster==0)

# Cluster weights

tc4_3_postprobs=data.frame(tc4_3$posterior)
colnames(tc4_3_postprobs)=c("p1","p2","p3","p4")
res_tc4_3=data.frame(data_k7_cleaned,tc4_3_postprobs)
res_tc4_3$cluster=apply(res_tc4_3[,c("p1","p2","p3","p4")],1,which.max)
res_tc4_3$cluster[tc4_3_outliers]=0

head(round(res_tc4_3[,c("p1","p2","p3","p4")],3))
head(res_tc4_3$cluster)

table(res_tc4_3$GRP,res_tc4_3$cluster)

# Representation on the PC space
tc4_3.pca <- prcomp(features_k7, scale = TRUE)

tc4_3.pca.eig.val<-factoextra::get_eigenvalue(tc4_3.pca)
tc4_3.pca.eig.val

# Add two first PC to the data
res_tc4_3$PC1=tc6_3.pca$x[,1]
res_tc4_3$PC2=tc6_3.pca$x[,2]

# Plot
max_prob <- apply(res_tc4_3[,c("p1","p2","p3","p4")], 1, max)
cex_values <- 1+(max_prob)^10
x11()
plot(res_tc4_3[,c("PC1","PC2")], 
     col = res_tc4_3$cluster + 1, 
     pch = 19,
     cex = cex_values,
     xlab = "PC1", 
     ylab = "PC2", 
     main = "PCA of the features")

legend("topright", 
       legend = c("Outlier","Cluster 1","Cluster 2","Cluster 3","Cluster 4"), 
       col = 1:5, 
       pch = 19, 
       cex = 0.8)

