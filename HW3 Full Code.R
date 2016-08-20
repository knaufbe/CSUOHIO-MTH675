#Create normal random variables X1~N(0,1) and X2={-X1 if -1 leq X1 leq 1, X1 otherwise
#Show that X2~N(0,1)
#Show that X1 and X2 do not have a bivariate normal distribution.


x1<-rnorm(1000,0,1)
x2<-x1
index<-which(abs(x1)>1)
x2[index]=-x2[index]
x<-cbind(x1,x2)
mvqq<-function(x){
	p <- ncol(x); n<-nrow(x)
	S <- var(x)
	xbar <- apply(x, 2, mean)
	D2 <- mahalanobis(x, xbar, S)
	qqplot(qchisq(ppoints(n), df=p), D2, pch=21, bg="blue")
	abline(0,1,col=2)
}
par(mfrow = c(1, 3))
qqnorm(x[,1],main="X1"); qqline(x[,1],col=2)
qqnorm(x[,1],main="X2"); qqline(x[,1],col=2)
mvqq(x);title("X1 and X2 Chi-square Q-Q Plot")

##Sketching solid ellipsoids for three given matrices

A = matrix(c(5, 4, 4, 5),2,2)
B = matrix(c(5, -4, -4, 5),2,2)
C = matrix(c(3,0,0,3),2,2)
install.packages("car",dependencies=TRUE)
SA <- function(X,add=FALSE,data.plot=TRUE)
{
	# sample mean and covariance
	#==============================
	xbar<-c(2,1)
	par(bg='white')
	library(car)
	ellipse(xbar,X,1,add=add,xlab="X1",ylab="X2",grid=FALSE,fill=T) #c=1
	#if (data.plot)
	#points(X[,1],X[,2],pch=20,col=4)
	# eigendecomposition
	#===============================
	e<-eigen(X)
	arrows(xbar[1],xbar[2],xbar[1]+e$vectors[1,1]*sqrt(e$values[1]),xbar[2]+e$vectors[2,1]*sqrt(e$values[1]),length=.1,col='red',lwd=2)
	arrows(xbar[1],xbar[2],xbar[1]+e$vectors[1,2]*sqrt(e$values[2]),xbar[2]+e$vectors[2,2]*sqrt(e$values[2]),length=.1,col='green',lwd=2)
	e
}

SA(A);title("Matrix 1");curve(x-1,-1,4.5,add=TRUE,type="l", lty=2, lwd=1,col='blue');curve(-x+3,-1,4.5,add=TRUE,type="l", lty=2, lwd=1,col='orange');legend(-.4,3,c(expression(y==x-1), expression(y==-x+3), expression(sqrt(lambda[1])*e[1]==bgroup("[",list(frac(3*sqrt(2),2), frac(3*sqrt(2),2)),"]")),expression(sqrt(lambda[2])*e[2]==bgroup("[",list(frac(sqrt(2),2), frac(-sqrt(2),2)),"]"))), cex=.7, col=c('blue','orange','red','green'),lty=c(2,2,1,1), lwd=c(2,2,1,1),bty = "n");
SA(B);title("Matrix 2");curve(x-1,-1,4.5,add=TRUE,type="l", lty=2, lwd=1,col='blue');curve(-x+3,-1,4.5,add=TRUE,type="l", lty=2, lwd=1,col='orange');legend("left",c(expression(y==x-1), expression(y==-x+3), expression(sqrt(lambda[1])*e[1]==bgroup("[",list(frac(-3*sqrt(2),2), frac(3*sqrt(2),2)),"]")),expression(sqrt(lambda[2])*e[2]==bgroup("[",list(frac(-sqrt(2),2), frac(-sqrt(2),2)),"]"))), cex=.7, col=c('blue','orange','red','green'),lty=c(2,2,1,1), lwd=c(2,2,1,1),bty = "n");
SA(C);title("Matrix 3");abline(v=2,type="l", lty=2, lwd=1,col='blue');abline(h=1,type="l", lty=2, lwd=1,col='orange');legend("topleft",c(expression(x==2), expression(y==1), expression(sqrt(lambda[1])*e[1]==bgroup("[",list(0, sqrt(3)),"]")),expression(sqrt(lambda[2])*e[2]==bgroup("[",list(-sqrt(3), 0),"]"))), cex=.85, col=c('blue','orange','red','green'),lty=c(2,2,1,1), lwd=c(2,2,1,1),bty = "n");


## Full analysis of data
x<-read.table("air_pollution.txt")
xnam<-cbind("Wind","Solar Radiation","CO","NO","NO2","O3","HC")
colnames(x)<-xnam

#summary statistics
summary(x)


# box plots
par(mfrow = c(1, 7))
boxplot(x[,1],main="Wind")
boxplot(x[,2],main="Solar Radiation")
boxplot(x[,3],main="CO")
boxplot(x[,4],main="NO")
boxplot(x[,5],main="NO2")
boxplot(x[,6],main="O3")
boxplot(x[,7],main="HC")


# histogram plots
par(mfrow = c(1, 7))
hist(x[,1],main="Wind",xlab="",col="blue", border="pink")
hist(x[,2],main="Solar Radiation",xlab="",col="blue", border="pink")
hist(x[,3],main="CO",xlab="",col="blue", border="pink")
hist(x[,4],main="NO",xlab="",col="blue", border="pink")
hist(x[,5],main="NO2",xlab="",col="blue", border="pink")
hist(x[,6],main="O3",xlab="",col="blue", border="pink")
hist(x[,7],main="HC",xlab="",col="blue", border="pink")


# Q-Q plots 
par(mfrow = c(2, 4))
qqnorm(x[,1],main="Wind"); qqline(x[,1],col=2)
qqnorm(x[,2],main="Solar Radiation"); qqline(x[,2],col=2)
qqnorm(x[,3],main="CO"); qqline(x[,3],col=2)
qqnorm(x[,4],main="NO"); qqline(x[,4],col=2)
qqnorm(x[,5],main="NO2"); qqline(x[,5],col=2)
qqnorm(x[,6],main="O3"); qqline(x[,6],col=2)
qqnorm(x[,7],main="HC"); qqline(x[,7],col=2)

#shapiro tests
lshap <- lapply(x, shapiro.test)
lres <- sapply(lshap, `[`, c("statistic","p.value"))
lres


### Multivariate chi-square q-q plot ####
mvqq<-function(x){
	p <- ncol(x); n<-nrow(x)
	S <- var(x)
	xbar <- apply(x, 2, mean)

	D2 <- mahalanobis(x, xbar, S)
	qqplot(qchisq(ppoints(n), df=p), D2, pch=21, bg="blue")
	abline(0,1,col=2)
}

mvqq(x)
title("A Chi-square Q-Q Plot")

#BoxCox Transformation
summary(p3<-powerTransform(x))
coef(p3, round=TRUE)
transformedY <- bcPower(x,coef(p3,round=TRUE))
transformedY

##Check for multivariate normality
mvqq(transformedY)
title("A Chi-square Q-Q Plot")


#redo summary and shape
summary(transformedY)


# histogram plots
par(mfrow = c(1, 7))
hist(transformedY[,1],main="Wind",xlab="",col="blue", border="pink")
hist(transformedY[,2],main="Solar Radiation",xlab="",col="blue", border="pink")
hist(transformedY[,3],main="CO",xlab="",col="blue", border="pink")
hist(transformedY[,4],main="log(NO)",xlab="",col="blue", border="pink")
hist(transformedY[,5],main="log(NO2)",xlab="",col="blue", border="pink")
hist(transformedY[,6],main="log(O3)",xlab="",col="blue", border="pink")
hist(transformedY[,7],main="HC",xlab="",col="blue", border="pink")

#correlations
cor(transformedY)

#covariance matrix
cov(transformedY)

## 2-d scatter plots
pairs(transformedY)

# Pricipal Components Analysis
# entering raw data and extracting PCs
# from the correlation matrix
fit <- princomp(transformedY, cor=TRUE)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
fit$scores # the principal components
par(mfrow=c(1,2))
plot(fit) # scree plot
biplot(fit)
	
