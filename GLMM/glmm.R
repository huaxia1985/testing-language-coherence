library(MCMCglmm)
setwd("~/GLMM")

#prepare data
data <- read.csv("data.csv")
var.name <- colnames(data)[-c(1:4)]
data$FAMILY <- as.factor(data$FAMILY)
data$AGE <- 2020-data$DOB
idx <- which(data$GENERATION=="Gen1")
#set reference level as Generation 1
data$GENERATION <- as.factor(data$GENERATION)

#fit GLMM models to each variable
m <- vector("list",length(var.name))
DIC <- numeric(length(var.name))
for (i in 1:length(var.name)) {
	var <- paste("data$",var.name[i],sep="")
	data$tmp <- eval(parse(text=var))
	#set reference pattern as the most often used patterns in Generation 1
	ref <- table(data$tmp[idx])
	ref <- names(ref)[which(ref==max(ref))[1]]
	data$tmp <- as.factor(data$tmp)
	data$tmp <- relevel(data$tmp,ref=ref)
	k <- length(levels(data$tmp))
	if (k>2) {
	IJ <- (1/k)*(diag(k-1)+matrix(1,k-1,k-1))
	#Model1: fixed variable: generation; random variable: family
	prior <- list(R=list(V=IJ,fix=1,n=k-1),G=list(G1=list(V=diag(k-1),n=k-1)),B=list(mu=rep(0,3*(k-1)),V=kronecker(IJ,diag(3))*(1+pi^2/3)))
	m1 <- MCMCglmm(tmp ~ trait:GENERATION-1,
	              random = ~us(trait):FAMILY,
	              rcov = ~us(trait):units,
	              prior = prior,
	              burnin=15000,
	              nitt=40000,
	              family="categorical",
	              data=data,verbose=F)
	#Model2: fixed variable: age + intercept; random variable: family
	prior <- list(R=list(V=IJ,fix=1,n=k-1),G=list(G1=list(V=diag(k-1),n=k-1)),B=list(mu=rep(0,2*(k-1)),V=kronecker(IJ,diag(2))*(1+pi^2/3)))
	m2 <- MCMCglmm(tmp ~ AGE:trait+trait,
	              random = ~us(trait):FAMILY,
	              rcov = ~us(trait):units,
	              prior = prior,
	              burnin=15000,
	              nitt=40000,
	              family="categorical",
	              data=data,verbose=F)
	#Model3: fixed variable: age + intercept; random variable: family
	prior <- list(R=list(V=IJ,fix=1,n=k-1),G=list(G1=list(V=diag(k-1),n=k-1)),B=list(mu=rep(0,4*(k-1)),V=kronecker(IJ,diag(4))*(1+pi^2/3)))
	m[[i]] <- MCMCglmm(tmp ~ trait:GENERATION+AGE:trait-1,
	              random = ~us(trait):FAMILY,
	              rcov = ~us(trait):units,
	              prior = prior,
	              burnin=15000,
	              nitt=40000,
	              family="categorical",
	              data=data,verbose=F)
	} else {
		prior <- list(R=list(V=1,fix=1),G=list(G1=list(V=1,nu=0.002)),B=list(mu=rep(0,3),V=diag(3)*(1+pi^2/3)))
		m1 <- MCMCglmm(tmp ~ GENERATION-1,
	              random = ~FAMILY,
	              prior = prior,
	              burnin=15000,
	              nitt=40000,
	              family="categorical",
	              data=data,verbose=F)
	    prior <- list(R=list(V=1,fix=1),G=list(G1=list(V=1,nu=0.002)),B=list(mu=rep(0,2),V=diag(2)*(1+pi^2/3)))
		m2 <- MCMCglmm(tmp ~ AGE+1,
	              random = ~FAMILY,
	              prior = prior,
	              burnin=15000,
	              nitt=40000,
	              family="categorical",
	              data=data,verbose=F)
	    prior <- list(R=list(V=1,fix=1),G=list(G1=list(V=1,nu=0.002)),B=list(mu=rep(0,4),V=diag(4)*(1+pi^2/3)))
		m[[i]] <- MCMCglmm(tmp ~ GENERATION+AGE-1,
	              random = ~FAMILY,
	              prior = prior,
	              burnin=15000,
	              nitt=40000,
	              family="categorical",
	              data=data,verbose=F)
	}
	#calculate DIC difference between model I and II
	DIC[i] <- m1$DIC - m2$DIC
}

#get pattern information
#there are 5 key tables, with key1.csv gives the key for variables with one variants, key2.csv gives the key for variables with two variants, etc.
#Using key2.csv as an example, each variable is a row, column V1 and V2 give the first and the second variant types of each variable, columns V1.1 and V1.2 give pattern 1 for each variable, which is 0,1, meaning that the pattern only uses the second variant of the variable. Similarly, columns V2.1 and V2.2 give pattern 2 for each variable, which is 1,1, meaning that the pattern uses both variants.
variable <- vector("list",1)
for (i in 1:5) {
	filename <- paste0("~/key",i,".csv")
	key <- read.csv(filename)
	if (i==1) {
		variable[[1]]$type <- key[1,2:(i+1)]
		variable[[1]]$n <- (dim(key)[2]-(i+1))/i
		variable[[1]]$pattern <- matrix(key[1,(i+2):dim(key)[2]],i,variable[[1]]$n)
		names(variable) <- key[1,1]
	}
	for (j in ifelse(i==1,2,1):dim(key)[1]) {
		tmp <- vector("list",1)
		tmp[[1]]$type <- key[j,2:(i+1)]
		tmp[[1]]$n <- (dim(key)[2]-(i+1))/i
		tmp[[1]]$pattern <- matrix(key[j,(i+2):dim(key)[2]],i,tmp[[1]]$n)
		names(tmp) <- key[j,1]
		variable <- c(variable,tmp)
	}
}

#extract regression coefficients
family.contri <- NULL
gen2.contri <- NULL
gen3.contri <- NULL
age.contri <- NULL
rowname <- NULL
type <- NULL
for (i in c(1:length(var.name))) {
	var <- paste("data$",var.name[i],sep="")
	data$tmp <- eval(parse(text=var))
	ref <- table(data$tmp[idx])
	ref <- names(ref)[which(ref==max(ref))[1]]
	data$tmp <- as.factor(data$tmp)
	data$tmp <- relevel(data$tmp,ref=ref)
	k <- length(levels(data$tmp))
	tmp <- variable[[which(names(variable)==var.name[i])]]
	ref <- as.numeric(ref)
	pat <- c(1:dim(tmp$pattern)[2])[-ref]
	for (j in 1:(k-1)) {
		a <- m[[i]]$VCV[,1+(j-1)*k]/(m[[i]]$VCV[,1+(j-1)*k]+m[[i]]$VCV[,1+(k-1)^2+(j-1)*k]+pi^2/3)
		a <- summary(a)
		ic <- a$statistics[1]
		family.contri <- c(family.contri,ic)
		c2 <- ((16*sqrt(3))/(15*pi))^2
		a <- m[[i]]$Sol/sqrt(1+c2*m[[i]]$VCV[,1+(k-1)^2+(j-1)*k])
		a <- summary(a)
		gen2.contri <- c(gen2.contri,(a$statistics[k-1+j,1]-a$statistics[j,1])/sqrt(a$statistics[k-1+j,2]^2+a$statistics[j,2]^2))
		gen3.contri <- c(gen3.contri,(a$statistics[2*(k-1)+j,1]-a$statistics[k-1+j,1])/sqrt(a$statistics[2*(k-1)+j,2]^2+a$statistics[k-1+j,2]^2))
		age.contri <- c(age.contri,a$statistics[3*(k-1)+j,1]/a$statistics[3*(k-1)+j,2])
		tmp2 <- unlist(tmp$pattern[,pat[j]])-unlist(tmp$pattern[,ref])
		type <- c(type,paste(sapply(c(1:length(tmp2))[which(tmp2!=0)], function (j) paste0(ifelse(tmp2[j]>0,"+","-"),as.character(tmp$type[[j]]))),collapse=""))
		rowname <- c(rowname,var.name[i])
	}
}

#divide patterns into smaller number of types
typeG <- grep(pattern="-G",type,fixed=T)
typeGp <- grep(pattern="+G",type,fixed=T)
typeK <- unique(c(grep(pattern="-K",type,fixed=T),grep(pattern="-E",type,fixed=T),grep(pattern="-I",type,fixed=T)))
typeKp <- unique(c(grep(pattern="+K",type,fixed=T),grep(pattern="+E",type,fixed=T),grep(pattern="+I",type,fixed=T)))
type3 <- type
type3[typeK] <- "-K"
type3[typeKp] <- "+K"
type3[typeGp] <- "+G"
type3[typeG] <- "-G"

#plot Figure 5
library(ggplot2)
rowname2 <- gsub("_.*","",rowname)
rowname2 <- gsub("F","",rowname2)
dat <- data.frame(varname=rowname2,family.contri,gen2.contri,gen3.contri,age.contri,type=type3)
gg <- ggplot(dat,aes(x=gen2.contri,y=gen3.contri,xmin=-3.8,xmax=3.8,ymin=-3.8,ymax=3.8))
gg+geom_point(aes(color=type,size=family.contri))+geom_text(data=subset(dat,gen2.contri>=1.96 |gen3.contri>=1.96),aes(gen2.contri,gen3.contri,label=varname,cex=0.5))+geom_text(data=subset(dat,gen2.contri<=-1.96 |gen3.contri<=-1.96),aes(gen2.contri,gen3.contri,label=varname,cex=0.5))
