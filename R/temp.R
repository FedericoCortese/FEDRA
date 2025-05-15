source("Utils_FEDRA.R")

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

prv=robust_COSA(Y=Y,zeta0=.2,K=2,tol=NULL,
            n_outer=20,alpha=.1,
            verbose=F,knn=10,c=2,M=NULL)


temp=COSA_gap(Y,zeta_grid=seq(0.1,1,length.out=3),
              K_grid=2:3,
              tol=1e-4,n_outer=10,alpha=.1,verbose=F,
              B=3,n_cores=3,knn=10,c=2,M=NULL)


library(ggplot2)
ggplot(temp$gap_stats, aes(x = zeta0, y = GAP, color = factor(K), group = K)) +
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
