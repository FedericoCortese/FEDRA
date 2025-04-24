#install.packages("devtools")
library(devtools) # for loading the function install_github
#install_github("mkampert/rCOSA") # install rCOSA
library(rCOSA) # load rCOSA

data(ApoE3) # ?ApoE3
cosa_rslts <- cosa2(ApoE3)
hierclust(cosa_rslts$D) 

set.seed(123); N <- 100; P <- 50;
X <- matrix(rnorm(N*P), nrow = N, ncol = P)
#i <- sample(x = 1:N) # i conform Figure 2
i=1:30
#k <- sample(x = 1:P) # k conform Figure 2
k=1:45

X[i[1:15], k[1:15]] <- X[i[1:15], k[1:15]]*0.2 + 1.5
X[i[1:15], k[16:30]] <- X[i[1:15], k[16:30]]*0.2 - 1.5
X[i[16:30], k[16:30]] <- X[i[16:30], k[16:30]]*0.2 - 1.5
X[i[16:30], k[31:45]] <- X[i[16:30], k[31:45]]*0.2 + 1.5
X <- data.frame(scale(X))
dim(X)
Y=X

true_clust=matrix(0,nrow=nrow(Y),ncol=ncol(Y))
true_clust[1:15,1:15]=1
true_clust[1:15,16:30]=2
true_clust[16:30,16:30]=2
true_clust[16:30,31:45]=1

x11()
par(mfrow=c(5,5))
for(i in 1:25){
  plot(Y[,i],col=true_clust[,i]+1,pch=19,main=paste("Feature",i))
}

x11()
par(mfrow=c(5,5))
for(i in 26:50){
  plot(Y[,i],col=true_clust[,i]+1,pch=19,main=paste("Feature",i))
}

# s=sample(1:2,nrow(Y),replace=T)
# W=matrix(runif(ncol(Y)*2),nrow=2,ncol=ncol(Y))
# W=W/rowSums(W)

X <- sapply(X, function(x) as.numeric(x))
cosa.rslt <- cosa2(X)
