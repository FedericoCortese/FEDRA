evaluate_sparse_kmeans_silhouette <- function(X,
                                              K_vals,
                                              wb_vals,
                                              nstart = 20,
                                              silent = TRUE) {
  # 1) precompute full distance matrix once
  dist_X <- dist(X)
  
  # 2) build all (K, wbound) combinations
  combos <- expand.grid(K = K_vals, wbound = wb_vals)
  
  # 3) one lapply over all rows of combos
  res_list <- lapply(seq_len(nrow(combos)), function(i) {
    k <- combos$K[i]
    w <- combos$wbound[i]
    
    # run sparse k-means for this single w
    skm <- KMeansSparseCluster(X, K = k,
                               wbounds = w,
                               nstart   = nstart,
                               silent   = silent)
    # skm is a list of length 1 when `wbounds` is scalar
    cl <- skm[[1]]$Cs
    
    # compute silhouette and take its mean
    sil <- silhouette(cl, dist_X)
    avg_sil <- mean(sil[, "sil_width"])
    
    data.frame(K        = k,
               wbound   = w,
               avg_sil  = avg_sil)
  })
  
  # 4) bind into one data.frame and return
  results <- do.call(rbind, res_list)
  rownames(results) <- NULL
  return(results)
}


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

# weight_inv_exp_dist <- function(Y,
#                                 #Ymedoids,
#                                 #index_medoids, 
#                                 s, 
#                                 W, zeta) {
#   TT <- nrow(Y)
#   P <- ncol(Y)
#   
#   # 1. Normalizzazione Gower
#   # range_Y <- apply(Y, 2, function(col) {
#   #   r <- max(col) - min(col)
#   #   if (r == 0) 1 else r
#   # })
#   
#   sk=apply(Y,2,function(x)IQR(x)/1.35)
#   
#   # 2. Genera indici delle coppie (i < j)
#   pairs <- combn(TT, 2)
#   i_idx <- pairs[1, ]
#   j_idx <- pairs[2, ]
#   n_pairs <- ncol(pairs)
# 
#   # 3. Estrai le righe corrispondenti
#   Yi <- Y[i_idx, , drop = FALSE]
#   Yj <- Y[j_idx, , drop = FALSE]
#   diff <- abs(Yi - Yj)
#   #sk=apply(diff,2,function(x)sum(x)/TT^2)
#   
#   diff=sweep(diff, 2, sk, FUN = "/")
#   #diff=sweep(diff, 2, range_Y, FUN = "/")
#   
# 
#   # 4. Estrai direttamente i pesi W[si, ] e W[sj, ] in blocco
#   W_si <- W[s[i_idx], , drop = FALSE]
#   W_sj <- W[s[j_idx], , drop = FALSE]
#   # W_si <- W[i_idx, , drop = FALSE]
#   # W_sj <- W[j_idx, , drop = FALSE]
#   max_w <- pmax(W_si, W_sj)
#   
#   # 5. Calcola la distanza finale
#   weighted_exp <- exp(-diff / zeta) * max_w
#   dist_vals <- -zeta * log(rowSums(weighted_exp))
#   
#   # 6. Ricostruzione della matrice simmetrica
#   mat <- matrix(0, TT, TT)
#   mat[cbind(i_idx, j_idx)] <- dist_vals
#   mat[cbind(j_idx, i_idx)] <- dist_vals
#   
#   return(mat)
# }

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

# WCD=function(s,Y,K){
#   #TT <- nrow(Y)
#   P <- ncol(Y)
#   
#   wcd=matrix(0,nrow=K,ncol=P)
#   
#   # 1. Normalizzazione Gower
#   # range_Y <- apply(Y, 2, function(col) {
#   #   r <- max(col) - min(col)
#   #   if (r == 0) 1 else r
#   # })
#   sk=apply(Y,2,function(x)IQR(x)/1.35)
#   
#   for(i in 1:K){
#     Ys=Y[s==i,]
#     TTk <- nrow(Ys)
#     pairs <- combn(TTk, 2)
#     i_idx <- pairs[1, ]
#     j_idx <- pairs[2, ]
#     n_pairs <- ncol(pairs)
#     
#     Yi <- Ys[i_idx, , drop = FALSE]
#     Yj <- Ys[j_idx, , drop = FALSE]
#     diff <- abs(Yi - Yj)
#     #sk=apply(diff,2,function(x)sum(x)/TTk^2)
#     
#     diff=sweep(diff, 2, sk, FUN = "/")
#    # diff=sweep(diff, 2, range_Y, FUN = "/")
#     
#     for(p in 1:P){
#       mat <- matrix(0, TTk, TTk)
#       mat[cbind(i_idx, j_idx)] <- diff[,p]
#       mat[cbind(j_idx, i_idx)] <- diff[,p]
#       wcd[i,p]=mean(apply(mat,1,median))
#     }
#     # 
#     #wcd[i,]=colSums(diff)/TTk^2
#     #wcd[i,]=apply(diff,2,median)/TTk
#   }
#   return(wcd)
#   
# }

library(DescTools)

lof_star=function(x,knn){
  lof_x=DescTools::LOF(x, knn)
  mean_lof=mean(lof_x)
  sd_lof=sd(lof_x)
  lof_st=(lof_x-mean_lof)/sd_lof
  return(lof_st)
}

v_1=function(x,knn=10,c=2,M=NULL){
  
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

# Rcpp --------------------------------------------------------------------

# library(Rcpp)
# Rcpp::sourceCpp("weight_inv_exp_dist.cpp")
# Rcpp::sourceCpp("wcd.cpp")



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



COSA=function(Y,zeta0,K,tol=NULL,n_outer=20,alpha=.1,verbose=F){
  
  # Simple version of COSA (no outlier detection)
  
  library(Rcpp)
  Rcpp::sourceCpp("weight_inv_exp_dist.cpp")
  Rcpp::sourceCpp("wcd.cpp")
  
  P=ncol(Y)
  TT=nrow(Y)
  
  # best_loss <- NULL
  # best_s <- NULL
  # best_W = NULL
  
  W=matrix(1/P,nrow=K,ncol=P)
  W_old=W
  
  zeta=zeta0
  
  s=initialize_states(Y,K)
  
  for (outer in 1:n_outer){
    
    DW=weight_inv_exp_dist(as.matrix(Y),
                           s,
                           W,zeta)
    medoids=cluster::pam(x=DW,k=K,diss=TRUE)
    
    s=medoids$clustering
    
    # Compute weights
    
    Spk=WCD(s,as.matrix(Y),K)
    wcd=exp(-Spk/zeta0)
    W=wcd/rowSums(wcd)
    
    
    w_loss=sum(W*Spk)  
    eps_W=mean((W-W_old)^2)
    if (!is.null(tol)) {
      if (eps_W < tol) {
        break
      }
    }
    
    W_old=W
    zeta=zeta+alpha*zeta0
    
    if (verbose) {
      cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
    }
    
  }
  
  
  return(list(W=W,s=s,medoids=medoids,w_loss=w_loss))
}

COSA2=function(Y,K,zeta0,n_init=10,tol=NULL,n_outer=20,alpha=.1,verbose=F){
  
  # Same as COSA but with multiple initialization
  
  Y=as.matrix(Y)
  
  library(Rcpp)
  Rcpp::sourceCpp("weight_inv_exp_dist.cpp")
  Rcpp::sourceCpp("wcd.cpp")
  
  P=ncol(Y)
  TT=nrow(Y)
  
  
  # Function to do one initialization run
  run_one <- function(init_id) {
    # 1) init
    s      <- initialize_states(Y, K)
    W=matrix(1/P,nrow=K,ncol=P)
    W_old=W
    zeta   <- zeta0
    
    # 2) outer loop
    for (outer in seq_len(n_outer)) {
      DW   <- weight_inv_exp_dist(Y, s, W, zeta)
      med  <- pam(x = DW, k = K, diss = TRUE)
      s    <- med$clustering
      medoids_indx=med$medoids
      
      Spk  <- WCD(s, Y, K)
      wcd  <- exp(-Spk / zeta0)
      W    <- wcd / rowSums(wcd)
      
      w_loss <- sum(W * Spk)
      eps_W  <- mean((W - W_old)^2)
      if (verbose) {
        message(sprintf("[init %2d | iter %3d] eps_W=%.2e w_loss=%.2e",
                        init_id, outer, eps_W, w_loss))
      }
      if(!is.null(tol)){
        if (eps_W < tol) break 
      }
      
      W_old <- W
      zeta  <- zeta + alpha * zeta0
    }
    
    list(
      init       = init_id,
      cluster = s,
      weights    = W,
      zeta       = zeta,
      loss       = w_loss,
      iter       = outer,
      diss       = DW,
      medoids    = Y[medoids_indx,],
      medoids_indx = medoids_indx
    )
  }
  
  # Run n_init times with lapply
  results <- lapply(seq_len(n_init), run_one)
  
  # Extract losses and pick best
  losses <- vapply(results, `[[`, numeric(1), "loss")
  best   <- which.min(losses)
  best_res <- results[[best]]
  
  return(best_res)
}


robust_COSA=function(Y,zeta0,K,tol=NULL,
                     n_outer=20,alpha=.1,verbose=F,knn=10,c=2,M=NULL){
  
  # Da sistemare come in COSA2
  # Robust version of COSA
  library(Rcpp)
  Rcpp::sourceCpp("weight_inv_exp_dist.cpp")
  Rcpp::sourceCpp("wcd.cpp")
  
  P=ncol(Y)
  TT=nrow(Y)
  
  # best_loss <- NULL
  # best_s <- NULL
  # best_W = NULL
  
  W=matrix(1/P,nrow=K,ncol=P)
  W_old=W
  
  zeta=zeta0
  
  s=initialize_states(Y,K)
  
  for (outer in 1:n_outer){
    
    ## Clustering
    #for(inner in 1:n_inner){
    
    v1=v_1(W[s,]*Y,knn=knn,c=c,M=M)
    v2=v_1(Y,knn=knn,c=c,M=M)
    
    v=apply(cbind(v1,v2),1,min)
    
    #Compute distances
    DW=weight_inv_exp_dist(as.matrix(Y * v),
                           s,
                           W,zeta)
    medoids=cluster::pam(x=DW,k=K,diss=TRUE)
    #Ymedoids=Y[medoids$medoids,]
    s=medoids$clustering
    
    # Compute weights
    
    Spk=WCD(s,as.matrix(Y * v),K)
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
  
  
  return(list(W=W,s=s,medoids=medoids,w_loss=w_loss,v=v))
}

robust_COSA2=function(Y,K,zeta0,n_init=10,tol=NULL,n_outer=20,alpha=.1,verbose=F,
                      knn=10, c=2, M=NULL){
  
  # Same as COSA but with multiple initialization
  
  Y=as.matrix(Y)
  
  library(Rcpp)
  Rcpp::sourceCpp("weight_inv_exp_dist.cpp")
  Rcpp::sourceCpp("wcd.cpp")
  
  P=ncol(Y)
  TT=nrow(Y)
  
  
  # Function to do one initialization run
  run_one <- function(init_id) {
    # 1) init
    s      <- initialize_states(Y, K)
    W=matrix(1/P,nrow=K,ncol=P)
    W_old=W
    zeta   <- zeta0
    
    # 2) outer loop
    for (outer in seq_len(n_outer)) {
      
      v1=v_1(W[s,]*Y,knn=knn,c=c,M=M)
      v2=v_1(Y,knn=knn,c=c,M=M)
      
      v=apply(cbind(v1,v2),1,min)
      
      DW   <- weight_inv_exp_dist(as.matrix(Y * v), s, W, zeta)
      med  <- pam(x = DW, k = K, diss = TRUE)
      s    <- med$clustering
      
      Spk  <- WCD(s, as.matrix(Y * v), K)
      wcd  <- exp(-Spk / zeta0)
      W    <- wcd / rowSums(wcd)
      
      w_loss <- sum(W * Spk)
      eps_W  <- mean((W - W_old)^2)
      if (verbose) {
        message(sprintf("[init %2d | iter %3d] eps_W=%.2e w_loss=%.2e",
                        init_id, outer, eps_W, w_loss))
      }
      if(!is.null(tol)){
        if (eps_W < tol) break 
      }
      
      W_old <- W
      zeta  <- zeta + alpha * zeta0
    }
    
    list(
      init       = init_id,
      cluster = s,
      weights    = W,
      zeta       = zeta,
      loss       = w_loss,
      iter       = outer,
      diss       = DW,
      v=v
    )
  }
  
  # Run n_init times with lapply
  results <- lapply(seq_len(n_init), run_one)
  
  # Extract losses and pick best
  losses <- vapply(results, `[[`, numeric(1), "loss")
  best   <- which.min(losses)
  best_res <- results[[best]]
  
  return(best_res)
}

# COSA_hd=function(Y,zeta0,K,tol=NULL,n_outer=20,alpha=.1,verbose=F,Ts=NULL){
#   P=ncol(Y)
#   TT=nrow(Y)
#   
#   if(is.null(Ts)){
#     Ts=round(TT/2)
#   }
#   
#   # best_loss <- NULL
#   # best_s <- NULL
#   # best_W = NULL
#   
#   W=matrix(1/P,nrow=K,ncol=P)
#   W_old=W
#   
#   zeta=zeta0
#   
#   s=initialize_states(Y,K)
#   # Ymedoids=cluster::pam(Y,k=K)
#   # s=Ymedoids$clustering
#   # Ymedoids=Ymedoids$medoids
#   #s=sample(1:K,TT,replace = T)
#   
#   for (outer in 1:n_outer){
#     
#     subsample=sample(1:TT,Ts,replace = F)
#     Ys=Y[subsample,]
#     ss=s[subsample]
#     
#     ## Clustering
#     #for(inner in 1:n_inner){
#     
#     #Compute distances
#     DW_1=weight_inv_exp_dist(Ys,
#                              ss,
#                              W,zeta)
#     medoids=cluster::pam(x=DW_1,k=K,diss=TRUE)
#     Ymedoids=Ys[medoids$medoids,]
#     s=medoids$clustering
#     
#     loss_by_state=weight_inv_exp_dist_medoids(Y, Ymedoids, s, W, zeta)
#     s=apply(loss_by_state,1,which.min)
#     # Compute weights
#     
#     Spk=WCD(s,Y,K)
#     wcd=exp(-Spk/zeta0)
#     W=wcd/rowSums(wcd)
#     
#     #}
#     
#     w_loss=sum(W*Spk)  
#     eps_W=mean((W-W_old)^2)
#     if (!is.null(tol)) {
#       if (eps_W < tol) {
#         break
#       }
#     }
#     
#     W_old=W
#     zeta=zeta+alpha*zeta0
#     
#     # print(W)
#     # print(zeta)
#     # print(Spk)
#     # print(zeta0)
#     # print(range(DW))
#     
#     if (verbose) {
#       cat(sprintf('Outer iteration %d: %.6e\n', outer, eps_W))
#     }
#     
#   }
#   
#   
#   return(list(W=W,s=s,medoids=medoids,w_loss=w_loss))
# }


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

COSA_gap <- function(Y,
                     K_grid    = 2:6,
                     zeta_grid = seq(0.01, 1, .1),
                     tol       = NULL,
                     n_outer   = 20,
                     alpha     = .1,
                     verbose   = FALSE,
                     n_cores   = NULL,
                     B         = 10,
                     n_init=10) {
  
  require(foreach)
  require(doParallel)
  require(dplyr)
  
  # 1) build your parameter grid
  grid <- expand.grid(K = K_grid, zeta0 = zeta_grid, b = 0:B)
  
  # 2) set up cores
  if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # 3) parallel loop, with tryCatch around each call
  results_list <- foreach(i = seq_len(nrow(grid)),
                          .packages = c("cluster","Rcpp"),
                          .export   = c("Y","COSA2","initialize_states",
                                        "tol","n_outer",
                                        "alpha","n_init"),
                          .errorhandling = "pass") %dopar% {
                            
                            params <- grid[i, ]
                            K_val  <- params$K
                            zeta0  <- params$zeta0
                            b      <- params$b
                            
                            set.seed(1000*i + b)
                            
                            # permute or not
                            Y_input <- if (b == 0) Y else apply(Y, 2, sample)
                            permuted <- (b != 0)
                            
                            # tryCatch so errors don’t blow everything up
                            out <- tryCatch({
                              res <- COSA2(Y_input,
                                           zeta0   = zeta0,
                                           K       = K_val,
                                           tol     = tol,
                                           n_outer = n_outer,
                                           alpha   = alpha,
                                           verbose = FALSE,
                                           n_init  = n_init)
                              
                              list(
                                meta = data.frame(K        = K_val,
                                                  zeta0    = zeta0,
                                                  loss     = res$loss,
                                                  permuted = permuted),
                                cosa = if (!permuted)
                                  list(K       = K_val,
                                       zeta0   = zeta0,
                                       W       = res$W,
                                       s       = res$cluster
                                  )
                                else NULL,
                                error = NA_character_
                              )
                              
                            }, error = function(e) {
                              # on error, record the params and the message
                              list(
                                meta  = data.frame(K        = K_val,
                                                   zeta0    = zeta0,
                                                   loss     = NA_real_,
                                                   permuted = permuted),
                                cosa  = NULL,
                                error = e$message
                              )
                            })
                            
                            out
                          }
  
  stopCluster(cl)
  
  # 4) pull out all the meta–data and cosa results
  meta_df <- do.call(rbind, lapply(results_list, `[[`, "meta"))
  # cosa_results <- Filter(Negate(is.null),
  #                        lapply(results_list, `[[`, "cosa"))
  
  # 5) inspect errors, if any
  errors <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(x$meta, error = x$error, stringsAsFactors = FALSE)
  })) %>% filter(!is.na(error))
  
  if (nrow(errors) && verbose) {
    message("Some iterations failed. Inspect 'errors' data frame.")
    print(errors)
  }
  
  # 6) compute GAP statistics
  gap_stats <- meta_df %>%
    group_by(K, zeta0) %>%
    summarise(
      log_O             = log(loss[!permuted]),
      log_O_star_mean   = mean(log(loss[permuted]), na.rm = TRUE),
      se_log_O_star     = sd(log(loss[permuted]), na.rm = TRUE),
      GAP               = log_O_star_mean - log_O,
      .groups           = "drop"
    )
  
  list(
    gap_stats    = gap_stats,
    errors       = if (nrow(errors)) errors else NULL
  )
}
# 
# robust_COSA_gap <- function(Y,
#                      zeta_grid = seq(0.01, 1, .1),
#                      K_grid    = 2:6,
#                      tol       = NULL,
#                      n_outer   = 20,
#                      alpha     = .1,
#                      verbose   = FALSE,
#                      n_cores   = NULL,
#                      B         = 10,
#                      knn       = 10,
#                      c         = 2,
#                      M         = NULL) {
#   
#   require(foreach)
#   require(doParallel)
#   require(dplyr)
#   
#   # 1) build your parameter grid
#   grid <- expand.grid(K = K_grid, zeta0 = zeta_grid, b = 0:B)
#   
#   # 2) set up cores
#   if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
#   cl <- makeCluster(n_cores)
#   registerDoParallel(cl)
#   
#   # 3) parallel loop, with tryCatch around each call
#   results_list <- foreach(i = seq_len(nrow(grid)),
#                           .packages = c("cluster","Rcpp","DescTools"),
#                           .export   = c("Y","robust_COSA","initialize_states",
#                                         "v_1","lof_star","tol","n_outer",
#                                         "alpha","knn","c","M"),
#                           .errorhandling = "pass") %dopar% {
#                             
#                             params <- grid[i, ]
#                             K_val  <- params$K
#                             zeta0  <- params$zeta0
#                             b      <- params$b
#                             
#                             set.seed(1000*i + b)
#                             
#                             # permute or not
#                             Y_input <- if (b == 0) Y else apply(Y, 2, sample)
#                             permuted <- (b != 0)
#                             
#                             # tryCatch so errors don’t blow everything up
#                             out <- tryCatch({
#                               res <- robust_COSA(Y_input,
#                                                  zeta0   = zeta0,
#                                                  K       = K_val,
#                                                  tol     = tol,
#                                                  n_outer = n_outer,
#                                                  alpha   = alpha,
#                                                  verbose = FALSE,
#                                                  knn     = knn,
#                                                  c       = c,
#                                                  M       = M)
#                               
#                               list(
#                                 meta = data.frame(K        = K_val,
#                                                   zeta0    = zeta0,
#                                                   loss     = res$w_loss,
#                                                   permuted = permuted),
#                                 cosa = if (!permuted)
#                                   list(K       = K_val,
#                                        zeta0   = zeta0,
#                                        W       = res$W,
#                                        s       = res$s,
#                                        medoids = res$medoids$medoids,
#                                        v       = res$v)
#                                 else NULL,
#                                 error = NA_character_
#                               )
#                               
#                             }, error = function(e) {
#                               # on error, record the params and the message
#                               list(
#                                 meta  = data.frame(K        = K_val,
#                                                    zeta0    = zeta0,
#                                                    loss     = NA_real_,
#                                                    permuted = permuted),
#                                 cosa  = NULL,
#                                 error = e$message
#                               )
#                             })
#                             
#                             out
#                           }
#   
#   stopCluster(cl)
#   
#   # 4) pull out all the meta–data and cosa results
#   meta_df <- do.call(rbind, lapply(results_list, `[[`, "meta"))
#   cosa_results <- Filter(Negate(is.null),
#                          lapply(results_list, `[[`, "cosa"))
#   
#   # 5) inspect errors, if any
#   errors <- do.call(rbind, lapply(results_list, function(x) {
#     data.frame(x$meta, error = x$error, stringsAsFactors = FALSE)
#   })) %>% filter(!is.na(error))
#   
#   if (nrow(errors) && verbose) {
#     message("Some iterations failed. Inspect 'errors' data frame.")
#     print(errors)
#   }
#   
#   # 6) compute GAP statistics
#   gap_stats <- meta_df %>%
#     group_by(K, zeta0) %>%
#     summarise(
#       log_O             = log(loss[!permuted]),
#       log_O_star_mean   = mean(log(loss[permuted]), na.rm = TRUE),
#       se_log_O_star     = sd(log(loss[permuted]), na.rm = TRUE),
#       GAP               = log_O_star_mean - log_O,
#       .groups           = "drop"
#     )
#   
#   list(
#     gap_stats    = gap_stats,
#     cosa_results = cosa_results,
#     errors       = if (nrow(errors)) errors else NULL
#   )
# }
