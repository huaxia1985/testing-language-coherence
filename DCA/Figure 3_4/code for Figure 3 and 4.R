#the dca package is available from Williams et al. 2010, but the package is not published on R, so we need to manually load the package
source("/MuDICA/dca_package/loadStuff.R")

setwd("~/DCA/Figure 3_4")

#read in all data metrices
data_mat= read.csv("data.csv",header=T)
group_mat= read.csv("group.csv",header=T)
variable_mat=read.csv("variable.csv",header=T,row.name=1)
factor_mat = read.csv("factor.csv",header=T,row.name=1)

#choose PC axes to plot
axis1 = 1
axis2 = 2

#perform DCA
X = as.matrix(data_mat)
Y = as.matrix(group_mat)
Z = as.matrix(variable_mat)

cpsu <- colorPoints(Y)

ret_DCA <- dca(X,Y,Z)

F = ret_DCA$f
G = ret_DCA$g
l = ret_DCA$l
S = ret_DCA$S
spp = ret_DCA$spp
c = ret_DCA$c
tau = ret_DCA$tau

ret_Sup <- supplementaryObservationPoints(l,X,F)
fobs = ret_Sup$sup
profile = ret_Sup$profile
deltainv = ret_Sup$i

ret_SCos <- squaredCosines(l,fobs)
fo2 <- ret_SCos$fo2
d_obs <- ret_SCos$d_obs
c_obs <- ret_SCos$c_obs

assignments <- obsToCenter(Y,fobs,G)
Dsup <- assignments$sup
assigned <- assignments$assigned
confusion  <- assignments$confusion

#plotting variable contribution
var_contr1 <- Z%*%(ret_DCA$cj[,axis1])  #contribution to PC1
var_contr2 <- Z%*%(ret_DCA$cj[,axis2])  #contribution to PC2
idx1 <- which(var_contr1>(1/dim(var_contr1)[1]))
idx2 <- which(var_contr2>(1/dim(var_contr2)[1]))
par(mfrow=c(2,1))
barplot(var_contr1[idx1][order(var_contr1[idx1])],names.arg=row.names(var_contr1)[idx1[order(var_contr1[idx1])]],cex.names=0.1)
barplot(var_contr2[idx2][order(var_contr2[idx2])],names.arg=row.names(var_contr2)[idx2[order(var_contr2[idx2])]],cex.names=0.1)

#plot factor level contribution
fact_contri1 <- as.matrix(factor_mat) %*% ret_DCA$ci[,axis1]   #contribution to PC1
fact_contri2 <- as.matrix(factor_mat) %*% ret_DCA$ci[,axis2]   #contribution to PC2
idx1 <- which(fact_contri1>(1/dim(fact_contri1)[1]))
idx2 <- which(fact_contri2>(1/dim(fact_contri2)[1]))
par(mfrow=c(2,1))
barplot(fact_contri1[idx1][order(fact_contri1[idx1])],names.arg=colnames(fact_contri1)[idx1[order(fact_contri1[idx1])]],cex.names=0.3)
barplot(fact_contri2[idx2][order(fact_contri2[idx2])],names.arg=colnames(fact_contri2)[idx2[order(fact_contri2[idx2])]],cex.names=0.3)