#the dca package is available from Williams et al. 2010, but the package is not published on R, so we need to manually load the package
source("/MuDICA/dca_package/loadStuff.R")

setwd("~/DCA/Figure 2")

#read in all data metrices
data_mat= read.csv("data.csv",header=T)
group_mat= read.csv("group.csv",header=T)
variable_mat=read.csv("variable.csv",header=T,row.name=1)
#factor_mat = read.csv("factor.csv",header=T,row.name=1)

#choose confidence interval for tolerance ellipses around groups
interval=0.95
intervalPercentage = paste(100*interval," %")

#choose PC axes to plot
axis1 = 1
axis2 = 2
pc = "Principal Component: "

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

#plot PC space and groups
pc_text1 = paste(pc,axis1)
pc_text1 = paste(pc_text1, " tau: %")
pc_text1 = paste(pc_text1, round(tau[axis1],2))
pc_text2 = paste(pc,axis2)
pc_text2 = paste(pc_text2, " tau: %")
pc_text2 = paste(pc_text2, round(tau[axis2],2))

mmHelp <- minmaxHelper(G,fobs)
minMaxMatrix <- mmHelp$tmm
mmHelp <- minmaxHelper(minMaxMatrix,F)
minMaxMatrix <- mmHelp$tmm
minMaxList <- mmHelp$mml
minMaxList$minx = minMaxList$minx * 1.1
minMaxList$maxx = minMaxList$maxx * 1.1
minMaxList$miny = minMaxList$miny * 1.1
minMaxList$maxy = minMaxList$maxy * 1.1		

#plot speakers (dots) and groups (ellipses)
main="Generations";sub="Fixed Effect";xlab=pc_text1;ylab=pc_text2;plot_ellipses=0.95;
ellipsePlotter(fobs,G,cps=cpsu$cps,axis1,axis2,classes=Y,uniqueClasses=as.matrix(cpsu$uc),minMaxList,interval=plot_ellipses,main=main,sub=sub,xlab=xlab,ylab=ylab)
plotHelper(fobs,G,cps=cpsu$cps,axis1,axis2,classes=Y,uniqueClasses=as.matrix(cpsu$uc),TRUE,plot_points=1,plot_centers=1,minMaxList,main,sub,xlab,ylab)

#plot variants with higher than average contribution
idx1 <- which(ret_DCA$cj[,axis1]>(1/dim(ret_DCA$cj)[1]))  #important variables for PC1
idx2 <- which(ret_DCA$cj[,axis2]>(1/dim(ret_DCA$cj)[1]))  #important variables for PC2
idx <- intersect(idx1,idx2)
idx11 <- setdiff(idx1,idx)
idx22 <- setdiff(idx2,idx)
color1 <- numeric(length(idx1))
pch1 <- numeric(length(idx1))
type <- read.csv("type.csv",row.name=1)  #variable types
for (i in 1:length(idx1)){
	if (type[names(which(Z[,idx1[i]]==1)),]=="NL") {
		color1[i] <- "blueviolet"
		pch1[i] <- 17
	}
	if (type[names(which(Z[,idx1[i]]==1)),]=="VL") {
		color1[i] <- "blueviolet"
		pch1[i] <- 16
	}
	if (type[names(which(Z[,idx1[i]]==1)),]=="NG") {
		color1[i] <- "blueviolet"
		pch1[i] <- 2
	}
	if (type[names(which(Z[,idx1[i]]==1)),]=="VG") {
		color1[i] <- "blueviolet"
		pch1[i] <- 1
	}
	if (type[names(which(Z[,idx1[i]]==1)),]=="C") {
		color1[i] <- "forestgreen"
		pch1[i] <- 17
	}
}
color2 <- numeric(length(idx2))
pch2 <- numeric(length(idx2))
for (i in 1:length(idx2)){
    if (type[names(which(Z[,idx2[i]]==1)),]=="NL") {
        color2[i] <- "blueviolet"
        pch2[i] <- 17
    }
    if (type[names(which(Z[,idx2[i]]==1)),]=="VL") {
        color2[i] <- "blueviolet"
        pch2[i] <- 16
    }
    if (type[names(which(Z[,idx2[i]]==1)),]=="NG") {
        color2[i] <- "blueviolet"
        pch2[i] <- 2
    }
    if (type[names(which(Z[,idx2[i]]==1)),]=="VG") {
        color2[i] <- "blueviolet"
        pch2[i] <- 1
    }
    if (type[names(which(Z[,idx2[i]]==1)),]=="C") {
        color2[i] <- "forestgreen"
        pch2[i] <- 17
    }
}
points(ret_DCA$f[idx1,c(axis1,axis2)],cex=1,col=color1,pch=pch1)
points(ret_DCA$f[idx11,c(axis1,axis2)],cex=1,pch="-")
text(ret_DCA$f[idx1,c(axis1,axis2)]-0.05,rownames(ret_DCA$f)[idx1],cex=0.7)
points(ret_DCA$f[idx2,c(axis1,axis2)],cex=1,col=color2,pch=pch2)
points(ret_DCA$f[idx22,c(axis1,axis2)],cex=1,pch="|")
text(ret_DCA$f[idx2,c(axis1,axis2)]-0.05,rownames(ret_DCA$f)[idx2],cex=0.7)
