library(parallel)
library(Rcpp)
#library(gower)
library(poliscidata)
library(StatMatch)
library(cluster)
library(dplyr)
library(tidyr)

source("Utils_robCOSA.R")

TT=100
P=10

zeta0=seq(0.05,0.4,by=.05)
alpha=.1
K=3
lambda=seq(0,2,.25)

tol=1e-16
verbose=F
nseed=50

c=c(5,7.5,10)

mu=3
rho=0
nu=10
pers = 0
K_true=3
perc_out=.05
out_sigma=50

hp=expand.grid(
  seed=1:nseed,
  zeta0=zeta0,
  K=K,
  lambda=lambda,
  TT=TT,
  P=P,
  c=c,
  K_true=K_true
)

ncores=parallel::detectCores()-1

rel_=list()
rel_[[1]]=c(1,4,5)
rel_[[2]]=c(2,4,5)
rel_[[3]]=c(3,4,5)

thres_out=0
thres_feat_weight=.02

start <- Sys.time()

res_list_K3 <- mclapply(seq_len(nrow(hp)), function(i) {
  tryCatch({
    
    seed <- hp$seed[i]
    zeta0 <- hp$zeta0[i]
    K <- hp$K[i]
    lambda <- hp$lambda[i]
    TT <- hp$TT[i]
    P <- hp$P[i]
    c <- hp$c[i]
    K_true <- hp$K_true[i]
    
    simDat <- sim_data_stud_t(
      seed = seed,
      TT = TT,
      P = P,
      Pcat = NULL,
      Ktrue = K_true,
      mu = mu,
      rho = rho,
      nu = nu,
      pers = pers
    )
    
    simDat_sparse <- simulate_sparse_hmm(
      Y = simDat$SimData,
      rel_ = rel_,
      true_stat = simDat$mchain,
      perc_out = perc_out,
      out_sigma = out_sigma,
      seed = seed
    )
    
    fit <- robust_COSA(
      Y = as.matrix(simDat_sparse$Y),
      zeta0 = zeta0,
      K = K,
      tol = 1e-8,
      n_init = 10,
      n_outer = 10,
      alpha = 0.1,
      verbose = T,
      knn = 10,
      c = c,
      M = NULL,
      hd=F,
      outlier=T
    )
    
    est_s <- fit$s
    est_s[fit$v <= thres_out] <- 0
    W_ind <- fit$W > thres_feat_weight
    truth <- simDat_sparse$truth
    W_truth <- simDat_sparse$W_truth
    ARI_s <- mclust::adjustedRandIndex(est_s, truth)
    ARI_s
    table(est_s, truth)
    
    temp=pam(Y[,1],3)
    table(temp$clustering,truth)
    
    ARI_W <- tryCatch(
      mclust::adjustedRandIndex(W_ind, W_truth),
      error = function(e) 0
    )
    
    list(
      seed = seed,
      zeta0 = zeta0,
      K = K,
      lambda = lambda,
      TT = TT,
      P = P,
      c = c,
      W = fit$W,
      s = est_s,
      v = fit$v,
      truth = truth,
      ARI_s = ARI_s,
      ARI_W = ARI_W
    )
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    list(error = e$message)
  })
}, mc.cores = ncores)

end <- Sys.time()