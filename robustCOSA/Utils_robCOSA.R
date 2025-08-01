simulate_sparse_hmm <- function(Y,
                                rel_,
                                true_stat,
                                perc_out   = 0.02,
                                out_sigma  = 100,
                                seed       = NULL) {
  # Y         : T x P data matrix or data.frame
  # rel_      : list of length K; rel_[[k]] is a vector of features in state k
  # true_stat : integer vector length T with values in 1:K
  # perc_out  : fraction of rows to turn into outliers
  # out_sigma : sd of Gaussian noise added for outliers
  # seed      : optional RNG seed for reproducibility
  
  if(!is.null(seed)) set.seed(seed)
  Y <- as.matrix(Y)
  TT  <- nrow(Y)
  P  <- ncol(Y)
  
  # 1) invert the rel_ list: for each feature p, which states mention p?
  inv_rel <- invert_rel(rel_, P)
  
  # 2) Irrelevant features = those never mentioned in rel_
  irrelevant <- which(vapply(inv_rel, length, integer(1)) == 0)
  if(length(irrelevant) > 0) {
    # permute their rows globally
    Y[, irrelevant] <- Y[sample(TT), irrelevant]
  }
  
  # 3) Relevant features = those that appear in at least one state
  relevant <- which(vapply(inv_rel, length, integer(1)) > 0)
  for(p in relevant) {
    relevant_states <- inv_rel[[p]]
    # indices of rows belonging to any of those states
    idx_in_state <- which(true_stat %in% relevant_states)
    # rows not in those states:
    idx_out_state <- setdiff(seq_len(TT), idx_in_state)
    if(length(idx_out_state) > 1) {
      # permute only those rows of column p
      Y[idx_out_state, p] <- sample(Y[idx_out_state, p])
    }
  }
  
  # 4) Introduce outliers
  N_out <- ceiling(TT * perc_out)
  t_out <- sort(sample(seq_len(TT), N_out))
  # add Gaussian noise to all features of those rows
  Y[t_out, ] <- Y[t_out, ] + matrix(rnorm(N_out * P, 0, out_sigma),
                                    nrow = N_out, ncol = P)
  
  # 5) update truth: set outlier rows to 0
  new_truth <- true_stat
  new_truth[t_out] <- 0L
  
  W_truth <- matrix(FALSE, nrow = K, ncol = P)
  
  for (k in seq_len(K)) {
    W_truth[k, rel_[[k]]] <- TRUE
  }
  
  list(
    Y          = Y,
    truth      = new_truth,
    out_indices = t_out,
    W_truth = W_truth
  )
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

invert_rel <- function(rel_, P) {
  # rel_: list of length K, each rel_[[k]] is a vector of indices in 1:P
  # P    : total number of elements
  #
  # returns a list inv of length P, where
  #   inv[[i]] = all k such that i %in% rel_[[k]]
  
  # initialize with empty integer vectors
  inv <- vector("list", P)
  for(i in seq_len(P)) inv[[i]] <- integer(0)
  
  # for each group k, append k to every member i in rel_[[k]]
  for(k in seq_along(rel_)) {
    members <- rel_[[k]]
    for(i in members) {
      inv[[i]] <- c(inv[[i]], k)
    }
  }
  
  inv
}

order_states_condMed=function(y,s){
  
  # This function organizes states by assigning 1 to the state with the smallest conditional median for vector y
  # and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  
  #Slong=c(t(S))
  # condMeans=sort(tapply(y,Slong,mean,na.rm=T))
  condMed=sort(tapply(y,s,median,na.rm=T))
  
  states_temp=match(s,names(condMed))
  
  #states_temp=matrix(states_temp,nrow=nrow(S),byrow = T)
  
  return(states_temp)
}

Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}

robust_COSA <- function(Y,
                               zeta0,
                               K,
                               tol     = 1e-16,
                               n_init  = 5,
                               n_outer = 20,
                               alpha   = 0.1,
                               verbose = FALSE,
                               knn     = 10,
                               c       = 10,
                               M       = NULL,
                               mif=NULL,
                               hd=F,
                               n_hd=NULL,
                               outlier=T) {
  
  P  <- ncol(Y)
  TT <- nrow(Y)
  
  library(Rcpp)
  
  Rcpp::sourceCpp("robCOSA.cpp")
  
  if(!is.matrix(Y)){
    Y=as.matrix(
      data.frame(
        lapply(Y, function(col) {
          if (is.factor(col)||is.character(col)) as.integer(as.character(col))
          else              as.numeric(col)
        })
      )
    )
  }
  
  # Check if which features are constant
  constant_features <- apply(Y, 2, function(col) length(unique(col)) == 1)
  
  # Transform into noise (?)
  Y[,constant_features]=rnorm(TT,sd=10)
  
  # Standardize
  Y=scale(Y)
  
  Gamma <- lambda * (1 - diag(K))
  
  run_one <- function(init_id) {
    # 1) initial W, zeta, s
    W        <- matrix(1/P, nrow=K, ncol=P)
    W_old    <- W
    zeta     <- zeta0
    loss_old <- Inf
    #s = sample(1:K,TT,replace=T)
    s        <- initialize_states(Y, K)
    
    for (outer in seq_len(n_outer)) {
      
      if(!outlier){
        v=rep(1,TT)
      }
      else{
        v1 <- v_1(W[s, , drop=FALSE] * Y, knn=knn, c=c, M=M)
        v2 <- v_1(Y,                     knn=knn, c=c, M=M)
        v  <- pmin(v1, v2)
      }
      
      
      if(hd){
        if(is.null(n_hd)){
          n_hd=500
        }
        
        sel_idx=sort(sample(1:TT,n_hd,replace=F))
        Y_search=Y[sel_idx,]
        
        Y_search=as.matrix(Y_search*v[sel_idx])
        
      }
      
      else{
        Y_search=as.matrix(Y * v)
        sel_idx=1:TT
      }
      
      DW      <- weight_inv_exp_dist(Y_search, s[sel_idx], W, zeta)
      pam_out <- cluster::pam(DW, k=K, diss=TRUE)
      #medoids <- pam_out$id.med
      medoids=sel_idx[pam_out$id.med]
      s_old <- s
      # 4) build loss-by-state
      if(!hd){
        loss_by_state <- DW[, medoids, drop=FALSE]  # TT x K
        s=pam_out$clustering
        loss=as.numeric(pam_out$objective[2])
      }
      
      else{
        # Not working
        loss_by_state <- weight_inv_exp_dist(Y=as.matrix(Y * v),
                                             s=s,W=W,zeta=zeta,
                                             medoids=medoids)
        s=apply(loss_by_state,1,which.min)
        loss  <- loss_by_state[1, s[1]]
      }
      
      
      # 7) must have all K states or revert
      if (length(unique(s)) < K) {
        s <- s_old
        break
      }
      
      # 8) loss‐convergence
      if (!is.null(tol) && (loss_old - loss) < tol) break
      loss_old <- loss
      
      # 9) update W via WCD + exp
      if(hd){
        Spk <- WCD(s[sel_idx], as.matrix(Y_search * v[sel_idx]), K)
      }
      else{
        Spk <- WCD(s, as.matrix(Y * v), K)
      }
      
      wcd <- exp(-Spk / zeta0)
      W   <- wcd / rowSums(wcd)
      
      # 10) W‐convergence
      epsW <- mean((W - W_old)^2)
      if (!is.null(tol) && epsW < tol) break
      W_old <- W
      
      # 11) bump zeta
      zeta <- zeta + alpha * zeta0
      
      if (verbose) {
        cat(sprintf("init %2d, outer %2d → loss=%.4e, epsW=%.4e, zeta=%.3f\n",
                    init_id, outer, loss, epsW, zeta))
      }
    }
    
    list(W      = W,
         s      = s,
         medoids= medoids,
         v      = v,
         loss   = loss)
  }
  
  # run n_init times, pick the one with smallest loss
  res_list <- lapply(seq_len(n_init), run_one)
  losses   <- vapply(res_list, `[[`, numeric(1), "loss")
  best_run <- res_list[[ which.min(losses) ]]
  
  best_s   <- best_run$s
  best_loss<- best_run$loss
  best_W = best_run$W
  best_medoids  <- Y[best_run$medoids,]
  best_v <- best_run$v
  
  # Most important features (mif)
  if(is.null(mif)){
    mif=which.max(apply(best_W,2,sum))
  }
  
  # Re‐order states based on most important feature state-conditional median
  new_best_s <- order_states_condMed(Y[, mif], best_s)
  
  tab <- table(best_s, new_best_s)
  new_order <- apply(tab, 1, which.max)
  
  best_W <- best_W[new_order,]
  
  ret_list=list(W = best_W,
                s = new_best_s,
                medoids = best_medoids,
                v = best_v,
                loss = best_loss,
                zeta0 = zeta0,
                c = c,
                knn=knn,
                M = M)
  
  return(ret_list)
  
}