
# Big dataset -------------------------------------------------------------

#load("D:/CNR/OneDrive - CNR/Brancati-Cortese/data and results/data_features_all.Rdata")
load("C:/Users/federico/OneDrive - CNR/Brancati-Cortese/data and results/data_features_all.Rdata")

# Togliere variabili con "iowa" e "bart" nel nome
data_features_all_2=data_features_all[,-grep("iowa|bart",colnames(data_features_all))]
features_k29=subset(data_features_all_2, select = -c(IDK, GRP))

source("Utils_FEDRA.R" )

features_k29_scaled=apply(features_k29,2,scale)

# Set K=4
COSA_4=COSA(Y=features_k29_scaled, zeta0=0.2, 
            K=4, alpha=0.1,verbose=T,n_outer=20)

W=COSA_4$W
colnames(W)=colnames(features_k29)

round(W,2)

res_COSA4=data.frame(data_features_all_2,clust=COSA_4$s)

tapply(res_COSA4$wcst_fail,res_COSA4$clust,mean)

table(res_COSA4$clust,res_COSA4$GRP)
