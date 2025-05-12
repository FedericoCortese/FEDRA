initialize_states <- function(Y, K) {
  n <- nrow(Y)
  
  ### Repeat the following few times?
  centr_indx=sample(1:n, 1)
  centroids <- Y[centr_indx, , drop = FALSE]  # Seleziona il primo centroide a caso
  
  closest_dist <- as.matrix(cluster::daisy(Y, metric = "gower"))
  closest_dist <- closest_dist[centr_indx,]
  
  for (i in 2:K) {
    prob <- closest_dist / sum(closest_dist)
    next_centr_indx <- sample(1:n, 1, prob = prob)
    next_centroid <- Y[next_centr_indx, , drop = FALSE]
    centroids <- rbind(centroids, next_centroid)
  }
  
  # Faster solution 
  dist_matrix <- StatMatch::gower.dist(Y, centroids)
  init_stats <- apply(dist_matrix, 1, which.min)
  
  return(init_stats)
}

weight_inv_exp_dist <- function(Y,
                                #Ymedoids,
                                #index_medoids, 
                                s, 
                                W, zeta) {
  TT <- nrow(Y)
  P <- ncol(Y)
  
  # 1. Normalizzazione Gower
  # range_Y <- apply(Y, 2, function(col) {
  #   r <- max(col) - min(col)
  #   if (r == 0) 1 else r
  # })
  
  sk=apply(Y,2,function(x)IQR(x)/1.35)
  
  # 2. Genera indici delle coppie (i < j)
  pairs <- combn(TT, 2)
  i_idx <- pairs[1, ]
  j_idx <- pairs[2, ]
  n_pairs <- ncol(pairs)

  # 3. Estrai le righe corrispondenti
  Yi <- Y[i_idx, , drop = FALSE]
  Yj <- Y[j_idx, , drop = FALSE]
  diff <- abs(Yi - Yj)
  #sk=apply(diff,2,function(x)sum(x)/TT^2)
  
  diff=sweep(diff, 2, sk, FUN = "/")
  #diff=sweep(diff, 2, range_Y, FUN = "/")
  

  # 4. Estrai direttamente i pesi W[si, ] e W[sj, ] in blocco
  W_si <- W[s[i_idx], , drop = FALSE]
  W_sj <- W[s[j_idx], , drop = FALSE]
  # W_si <- W[i_idx, , drop = FALSE]
  # W_sj <- W[j_idx, , drop = FALSE]
  max_w <- pmax(W_si, W_sj)
  
  # 5. Calcola la distanza finale
  weighted_exp <- exp(-diff / zeta) * max_w
  dist_vals <- -zeta * log(rowSums(weighted_exp))
  
  # 6. Ricostruzione della matrice simmetrica
  mat <- matrix(0, TT, TT)
  mat[cbind(i_idx, j_idx)] <- dist_vals
  mat[cbind(j_idx, i_idx)] <- dist_vals
  
  return(mat)
}

# WCD=function(s,Y,K){
#   Tk=table(s)
#   wcd=matrix(0,nrow=K,ncol=ncol(Y))
#   for(k in unique(s)){
#     for(p in 1:P){
#       temp=Y[s==k,p]
#       dist_temp=abs(outer(temp, temp, "-"))
#       dist_temp[upper.tri(dist_temp)] <- 0
#       dist_temp=dist_temp/diff(range(temp))
#       wcd[k,p]=sum(dist_temp)/Tk[k]^2
#     }
#   }
#   return(wcd)
# }

weight_inv_exp_dist_medoids <- function(Y, Ymedoids, s, W, zeta) {
  indx_num <- sapply(Y, is.numeric)
  
  TT <- nrow(Y)
  K <- nrow(Ymedoids)
  P <- ncol(Y)
  
  Y_orig <- Y
  Ymedoids_orig <- Ymedoids
  
  # Continuous variables
  Y <- Y[, indx_num, drop = FALSE]
  Ymedoids_cont <- Ymedoids[, indx_num, drop = FALSE]
  sk <- apply(Y, 2, function(x) IQR(x) / 1.35)
  
  # Pre-allocate output matrix
  mat <- matrix(0, nrow = TT, ncol = K)
  
  # Ciclo su tutti i medoids
  for (k in 1:K) {
    diff_cont <- abs(sweep(Y, 2, as.numeric(Ymedoids_cont[k,]), FUN = "-"))
    
    diff_cont <- sweep(diff_cont, 2, sk, FUN = "/") # normalizzazione
    diff <- matrix(0, nrow = TT, ncol = P)
    diff[, indx_num] <- as.matrix(diff_cont[, names(indx_num)[indx_num]])
    
    # Categorical variables
    if (sum(indx_num) != P) {
      Y_cat <- Y_orig[, !indx_num, drop = FALSE]
      medoid_cat <- Ymedoids_orig[k, !indx_num, drop = FALSE]
      diff_cat <- sweep(as.matrix(Y_cat), 2, as.matrix(medoid_cat), FUN = "!=") * 1
      diff[, !indx_num] <- as.matrix(diff_cat[, names(!indx_num)[!indx_num]])
    }
    
    # Pesatura con pesi W
    W_si <- W[s, , drop = FALSE]
    W_sj <- matrix(W[k, ], nrow = TT, ncol = ncol(W), byrow = TRUE)
    max_w <- pmax(W_si, W_sj)
    
    weighted_exp <- exp(-diff / zeta) * max_w
    dist_vals <- -zeta * log(rowSums(weighted_exp))
    
    mat[, k] <- dist_vals
  }
  
  return(mat)  # matrice T x K
}

WCD=function(s,Y,K){
  #TT <- nrow(Y)
  P <- ncol(Y)
  
  wcd=matrix(0,nrow=K,ncol=P)
  
  # 1. Normalizzazione Gower
  # range_Y <- apply(Y, 2, function(col) {
  #   r <- max(col) - min(col)
  #   if (r == 0) 1 else r
  # })
  sk=apply(Y,2,function(x)IQR(x)/1.35)
  
  for(i in 1:K){
    Ys=Y[s==i,]
    TTk <- nrow(Ys)
    pairs <- combn(TTk, 2)
    i_idx <- pairs[1, ]
    j_idx <- pairs[2, ]
    n_pairs <- ncol(pairs)
    
    Yi <- Ys[i_idx, , drop = FALSE]
    Yj <- Ys[j_idx, , drop = FALSE]
    diff <- abs(Yi - Yj)
    #sk=apply(diff,2,function(x)sum(x)/TTk^2)
    
    diff=sweep(diff, 2, sk, FUN = "/")
   # diff=sweep(diff, 2, range_Y, FUN = "/")
    
    for(p in 1:P){
      mat <- matrix(0, TTk, TTk)
      mat[cbind(i_idx, j_idx)] <- diff[,p]
      mat[cbind(j_idx, i_idx)] <- diff[,p]
      wcd[i,p]=mean(apply(mat,1,median))
    }
    # 
    #wcd[i,]=colSums(diff)/TTk^2
    #wcd[i,]=apply(diff,2,median)/TTk
  }
  return(wcd)
  
}

sim_data_stud_t=function(seed=123,
                         TT,
                         P,
                         Pcat,
                         Ktrue=3,
                         mu=1.5,
                         rho=0,
                         nu=4,
                         phi=.8,
                         pers=.95){
  
  
  MU=seq(mu, -mu, length.out=Ktrue)
  
  # Markov chain simulation
  x <- numeric(TT)
  Q <- matrix(rep((1-pers)/(Ktrue-1),Ktrue*Ktrue), 
              ncol = Ktrue,
              byrow = TRUE)
  diag(Q)=rep(pers,Ktrue)
  init <- rep(1/Ktrue,Ktrue)
  set.seed(seed)
  x[1] <- sample(1:Ktrue, 1, prob = init)
  for(i in 2:TT){
    x[i] <- sample(1:Ktrue, 1, prob = Q[x[i - 1], ])
  }
  
  # Continuous variables simulation
  Sigma <- matrix(rho,ncol=P,nrow=P)
  diag(Sigma)=1
  
  Sim = matrix(0, TT, P * Ktrue)
  SimData = matrix(0, TT, P)
  
  set.seed(seed)
  for(k in 1:Ktrue){
    # u = MASS::mvrnorm(TT,rep(mu[k],P),Sigma)
    u = mvtnorm::rmvt(TT, sigma = (nu-2)*Sigma/nu, df = nu, delta = rep(MU[k],P))
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:TT) {
    k = x[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
    #SimDataCat[i, ] = SimCat[i, (Pcat * k - Pcat + 1):(Pcat * k)]
  }
  
  SimData=data.frame(SimData)
  
  if(!is.null(Pcat)){
    for (j in 1:Pcat) {
      SimData[, j] <- get_cat_t(SimData[, j], x, MU, phi=phi, df = nu)
      SimData[, j]=factor(SimData[, j],levels=1:Ktrue)
    }  
  }
  
  
  
  return(list(
    SimData=SimData,
    mchain=x,
    TT=TT,
    P=P,
    K=Ktrue,
    Ktrue=Ktrue,
    pers=pers, 
    seed=seed))
  
}

# apply_noise_by_cluster <- function(Y, s, feat_list) {
#   Y_noised <- Y
#   K <- length(feat_list)
#   
#   for (k in 1:K) {
#     cluster_rows <- which(s == k)
#     if (length(cluster_rows) <= 1) next  # nulla da mescolare se solo una riga
#     all_features <- seq_len(ncol(Y))
#     irrelevant_feats <- setdiff(all_features, feat_list[[k]])
#     
#     for (j in irrelevant_feats) {
#       # mescola le osservazioni della colonna j solo tra le righe del cluster k
#       Y_noised[cluster_rows, j] <- sample(Y[, j],size=length(Y_noised[cluster_rows, j]))
#     }
#   }
#   
#   return(Y_noised)
# }

zeta0=0.3
alpha=.1
K=2
tol=1e-9
n_outer=15
verbose=T

TT=500
P=10

simDat=sim_data_stud_t(seed=123,
                       TT=TT,
                       P=P,
                       Pcat=NULL,
                       Ktrue=2,
                       mu=5,
                       rho=0,
                       nu=100,
                       phi=.8,
                       pers=0.1)

Y=simDat$SimData
true_stat=simDat$mchain

#plot(Y[,2],col=true_stat,pch=19)

# feat_list=list()
# feat_list[[1]]=1:10
# feat_list[[2]]=5:15
# feat_list[[3]]=10:20
# Y_noised=apply_noise_by_cluster(Y,simDat$mchain,feat_list)
# Y=Y_noised


nu=4
# State 1, only features 1,2 and 3 are relevant
indx=which(true_stat!=1)
Sigma <- matrix(0,ncol=P-3,nrow=P-3)
diag(Sigma)=5
Y[indx,-(1:3)]=mvtnorm::rmvt(length(indx), 
                             sigma = (nu-2)*Sigma/nu, 
                             df = nu, delta = rep(0,P-3))

# State 2, only features 3,4 and 5 are relevant
indx=which(true_stat!=2)
Y[indx,-(3:5)]=mvtnorm::rmvt(length(indx), 
                             sigma = (nu-2)*Sigma/nu, 
                             df = nu, delta = rep(0,P-3))

# # State 3, only features 5,6 and 7 are relevant
# indx=which(true_stat!=3)
# Y[indx,-(5:7)]=mvtnorm::rmvt(length(indx), 
#                              sigma = (nu-2)*Sigma/nu, 
#                              df = nu, delta = rep(0,P-3))

# Other features are not relevant
# Sigma <- matrix(0,ncol=5,nrow=5)
# diag(Sigma)=5
# Y[,6:10]=mvtnorm::rmvt(TT, 
#                        sigma = (nu-2)*Sigma/nu, 
#                        df = nu, delta = rep(0,5))

Sigma <- matrix(0,ncol=P-5,nrow=P-5)
diag(Sigma)=5
Y[,6:P]=mvtnorm::rmvt(TT, 
                       sigma = (nu-2)*Sigma/nu, 
                       df = nu, delta = rep(0,P-5))


# Introduce outliers
set.seed(1)
N_out=TT*0.02
t_out=sample(1:TT,size=N_out)
Y[t_out,]=Y[t_out,]+rnorm(N_out*P,0,10)

truth=simDat$mchain
truth[t_out]=0
  
x11()
par(mfrow=c(4,3))
for (i in 1:P) {
  plot(Y[, i], col=truth+1, pch=19,ylab=i)
}


COSA=function(Y,zeta0,K,tol=NULL,n_outer=20,alpha=.1,verbose=F){
  P=ncol(Y)
  TT=nrow(Y)
  
  # best_loss <- NULL
  # best_s <- NULL
  # best_W = NULL
  
  W=matrix(1/P,nrow=K,ncol=P)
  W_old=W
  
  zeta=zeta0
  
  s=initialize_states(Y,K)
  # Ymedoids=cluster::pam(Y,k=K)
  # s=Ymedoids$clustering
  # Ymedoids=Ymedoids$medoids
  #s=sample(1:K,TT,replace = T)
  
  for (outer in 1:n_outer){
    
    ## Clustering
    #for(inner in 1:n_inner){
      
      #Compute distances
      DW=weight_inv_exp_dist(Y,
                             s,
                             W,zeta)
      medoids=cluster::pam(x=DW,k=K,diss=TRUE)
      #Ymedoids=Y[medoids$medoids,]
      s=medoids$clustering
      
      # Compute weights
      
      Spk=WCD(s,Y,K)
      wcd=exp(-Spk/zeta0)
      W=wcd/rowSums(wcd)
      
    #}
    
    w_loss=sum(W*Spk)  
    eps_W=mean((W-W_old)^2)
    if (!is.null(tol)) {
      if (eps_W < tol) {
        break
      }
    }
    
    W_old=W
    zeta=zeta+alpha*zeta0
    
    # print(W)
    # print(zeta)
    # print(Spk)
    # print(zeta0)
    # print(range(DW))
    
    if (verbose) {
      cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
    }
    
  }
  
  
  return(list(W=W,s=s,medoids=medoids,w_loss=w_loss))
}

COSA_hd=function(Y,zeta0,K,tol=NULL,n_outer=20,alpha=.1,verbose=F,Ts=NULL){
  P=ncol(Y)
  TT=nrow(Y)
  
  if(is.null(Ts)){
    Ts=round(TT/2)
  }
  
  # best_loss <- NULL
  # best_s <- NULL
  # best_W = NULL
  
  W=matrix(1/P,nrow=K,ncol=P)
  W_old=W
  
  zeta=zeta0
  
  s=initialize_states(Y,K)
  # Ymedoids=cluster::pam(Y,k=K)
  # s=Ymedoids$clustering
  # Ymedoids=Ymedoids$medoids
  #s=sample(1:K,TT,replace = T)
  
  for (outer in 1:n_outer){
    
    subsample=sample(1:TT,Ts,replace = F)
    Ys=Y[subsample,]
    ss=s[subsample]
    
    ## Clustering
    #for(inner in 1:n_inner){
    
    #Compute distances
    DW_1=weight_inv_exp_dist(Ys,
                             ss,
                             W,zeta)
    medoids=cluster::pam(x=DW_1,k=K,diss=TRUE)
    Ymedoids=Ys[medoids$medoids,]
    s=medoids$clustering
    
    loss_by_state=weight_inv_exp_dist_medoids(Y, Ymedoids, s, W, zeta)
    s=apply(loss_by_state,1,which.min)
    # Compute weights
    
    Spk=WCD(s,Y,K)
    wcd=exp(-Spk/zeta0)
    W=wcd/rowSums(wcd)
    
    #}
    
    w_loss=sum(W*Spk)  
    eps_W=mean((W-W_old)^2)
    if (!is.null(tol)) {
      if (eps_W < tol) {
        break
      }
    }
    
    W_old=W
    zeta=zeta+alpha*zeta0
    
    # print(W)
    # print(zeta)
    # print(Spk)
    # print(zeta0)
    # print(range(DW))
    
    if (verbose) {
      cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
    }
    
  }
  
  
  return(list(W=W,s=s,medoids=medoids,w_loss=w_loss))
}

COSA_gap=function(Y,
                  zeta_grid=seq(0,1,.1),
                  K_grid=2:6,
                  tol=NULL,n_outer=20,alpha=.1,verbose=F,n_cores=NULL,
                  B=10, Ts=NULL){
  
  # B is the number of permutations
  # Ts is the sample size for the subsampling (as in CLARA), if NULL it uses the full sample
  
  grid <- expand.grid(K = K_grid, zeta0 = zeta_grid, b = 0:B)
  
  library(foreach)
  library(doParallel)
  
  if(is.null(n_cores)){
    n_cores <- parallel::detectCores() - 1
  } 
  
  # Set up cluster
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  if(is.null(Ts)){
    results_list <- foreach(i = 1:nrow(grid), .combine = 'list',
                            .packages = c("cluster"),
                            .multicombine = TRUE,
                            .export = c("Y", "COSA", "WCD", "weight_inv_exp_dist",
                                        "initialize_states",
                                        "grid", "tol", "n_outer", "alpha")) %dopar% {
                                          K_val <- grid$K[i]
                                          zeta_val <- grid$zeta0[i]
                                          b <- grid$b[i]
                                          
                                          set.seed(b + 1000 * i)
                                          
                                          if (b == 0) {
                                            Y_input <- Y
                                            permuted <- FALSE
                                          } else {
                                            Y_input <- apply(Y, 2, sample)
                                            permuted <- TRUE
                                          }
                                          
                                          res <- COSA(Y_input, zeta0 = zeta_val, K = K_val, tol = tol,
                                                      n_outer = n_outer, alpha = alpha, verbose = FALSE)
                                          
                                          list(
                                            meta = data.frame(K = K_val, zeta0 = zeta_val, 
                                                              loss = res$w_loss, permuted = permuted),
                                            cosa = if (!permuted) list(K = K_val, zeta0 = zeta_val, 
                                                                       W = res$W, s = res$s, 
                                                                       medoids = res$medoids$medoids) else NULL
                                          )
                                        }
  }
  else{
    results_list <- foreach(i = 1:nrow(grid), .combine = 'list',
                            .packages = c("cluster"),
                            .multicombine = TRUE,
                            .export = c("Y","COSA_hd", "WCD", "weight_inv_exp_dist",
                                        "initialize_states",
                                        "grid", "tol", "n_outer", "alpha","Ts")) %dopar% {
                                          K_val <- grid$K[i]
                                          zeta_val <- grid$zeta0[i]
                                          b <- grid$b[i]
                                          
                                          set.seed(b + 1000 * i)
                                          
                                          if (b == 0) {
                                            Y_input <- Y
                                            permuted <- FALSE
                                          } else {
                                            Y_input <- apply(Y, 2, sample)
                                            permuted <- TRUE
                                          }
                                          
                                          res <- COSA_hd(Y_input, zeta0 = zeta_val, K = K_val, tol = tol,
                                                      n_outer = n_outer, alpha = alpha, verbose = FALSE,
                                                      Ts=Ts)
                                          
                                          list(
                                            meta = data.frame(K = K_val, zeta0 = zeta_val, 
                                                              loss = res$w_loss, permuted = permuted),
                                            cosa = if (!permuted) list(K = K_val, zeta0 = zeta_val, 
                                                                       W = res$W, s = res$s, 
                                                                       medoids = res$medoids$medoids) else NULL
                                          )
                                        }
  }
  
  
  stopCluster(cl)
  
  # Flatten results
  meta_df <- do.call(rbind, lapply(results_list, `[[`, "meta"))
  cosa_results <- Filter(Negate(is.null), lapply(results_list, `[[`, "cosa"))
  
  
  # Compute GAP
  library(dplyr)
  gap_stats <- meta_df %>%
    group_by(K, zeta0) %>%
    summarise(
      log_O = log(loss[!permuted]),
      log_O_star_mean = mean(log(loss[permuted])),
      se_log_O_star=sd(log(loss[permuted])),
      GAP = log_O_star_mean - log_O,
      .groups = 'drop'
    )
  
  return(list(
    gap_stats = gap_stats,
    cosa_results = cosa_results
  ))
  
}

temp=COSA_gap(Y,zeta_grid=seq(0.1,1,length.out=10),
                  K_grid=2:3,
                  tol=1e-6,n_outer=10,alpha=.1,verbose=F,
                  B=10,Ts=round(TT*.2),n_cores=20)

ggplot2::ggplot(temp$gap_stats, aes(x = zeta0, y = GAP, color = factor(K), group = K)) +
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


temp_hd=COSA_hd(Y,.1,2,tol=1e-7,n_outer=10,alpha=.1,verbose=F,
                Ts=round(TT*.2))

table(temp_hd$s,simDat$mchain)
round(temp_hd$W,3)

library(DescTools)

lof_star=function(x,knn){
  lof_x=DescTools::LOF(x, knn)
  mean_lof=mean(lof_x)
  sd_lof=sd(lof_x)
  lof_st=(lof_x-mean_lof)/sd_lof
  return(lof_st)
}

v_1=function(x,knn,c=2,M=NULL){
  
  lof_st=lof_star(x,knn)
  
  if(is.null(M)){
    M=median(lof_st)+mad(lof_st)
  }
  
  v=rep(1,dim(x)[1])
  v[lof_st>=c]=0
  indx=which(M<lof_st&lof_st<c)
  v[indx]=(1-((lof_st[indx]-M)/(c-M))^2)^2
  return(v)
}

robust_COSA=function(Y,zeta0,K,tol,n_outer=20,alpha=.1,verbose=F,trim_lev=.1){
  P=ncol(Y)
  TT=nrow(Y)
  Y_orig=Y
  
  # best_loss <- NULL
  # best_s <- NULL
  # best_W = NULL
  
  W=matrix(1/P,nrow=K,ncol=P)
  W_old=W
  
  #v1=rep(0,TT)
  
  zeta=zeta0
  
  #s=initialize_states(Y,K)
  # Ymedoids=cluster::pam(Y,k=K)
  # s=Ymedoids$clustering
  # Ymedoids=Ymedoids$medoids
  s=sample(1:K,TT,replace = T)
  
  for (outer in 1:n_outer){
    
    ## Clustering
    #for(inner in 1:n_inner){
    
    #Compute distances
    DW=weight_inv_exp_dist(Y,
                           s,
                           W,zeta)
    
    # k med
    medoids=cluster::pam(x=DW,k=K,diss=TRUE)
    indx_med=medoids$medoids
    Ymedoids=Y[indx_med,]
    s=medoids$clustering
    
    # fuzzy k med
    # fkm=fclust::FKM.med (Y, k=K, m=1.5)
    # indx_med=fkm$medoid
    # Ymedoids=Y[fkm$medoid,]
    
    # Compute weights
    
    Spk=WCD(s,Y,K)
    wcd=exp(-Spk/zeta0)
    W=wcd/rowSums(wcd)
    
    # Trimming
    # Numero totale di osservazioni da escludere
    n_trim <- ceiling(TT * trim_lev)
    
    # Calcola la distanza di ogni osservazione dal proprio medoide
    d_from_medoid <- sapply(1:TT, function(i) {
      cluster <- s[i]
      medoid_index <- indx_med[cluster]
      DW[i, medoid_index]
    })
    
    # Ordina le distanze decrescenti e prendi gli n_trim indici più lontani
    indx_trim <- order(d_from_medoid, decreasing = TRUE)[1:n_trim]
    
    
    # LOF
    
    # for(i in 1:K){
    #   indx=which(s==i)
    #   Ys=Y[indx,]
    #   wYs=sweep(Ys, 2, W[i,], `*`)
    #   v1[indx]=v_1(wYs,knn=5,c=2)
    #   v2[indx]=v_1(Ys,knn=5,c=2)
    #   v1[indx]=pmin(v1[indx],v2[indx])
    #   Ys=sweep(Ys, 1, v1[indx], `*`)
    #   Y[indx,]=Ys
    # }
    
    
    #}
    
    eps_W=mean((W-W_old)^2)
    if (!is.null(tol)) {
      if (eps_W < tol) {
        break
      }
    }
    
    W_old=W
    zeta=zeta+alpha*zeta0
    
    print(W)
    # print(zeta)
    # print(Spk)
    # print(zeta0)
    # print(range(DW))
    
    if (verbose) {
      cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
    }
    
  }
  return(list(W=W,s=s,medoids=medoids,v1=v1))
}

plot_W=function(W){
  library(reshape)
  df <- as.data.frame(W)
  df$Cluster <- factor(paste0("Cluster_", 1:nrow(df)))
  
  # Riorganizziamo in formato lungo
  df_long <- melt(df, id.vars = "Cluster", variable.name = "Feature", value.name = "Weight")
  
  # Converti Feature in fattore per ordinare le colonne
  #df_long$Feature <- factor(df_long$Feature, levels = paste0("V", 1:ncol(fit$W)))
  
  # Bar plot
  library(ggplot2)
  p=ggplot2::ggplot(df_long, aes(x = variable, y = value, fill = Cluster)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Feature Weights by Cluster",
         x = "Feature",
         y = "Weight") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}


