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
# s=sample(1:2,nrow(Y),replace=T)
# W=matrix(runif(ncol(Y)*2),nrow=2,ncol=ncol(Y))
# W=W/rowSums(W)

X <- sapply(X, function(x) as.numeric(x))
cosa.rslt <- cosa2(X)
