source("AWNCut_fun.R")
#This sets up the initial parameters for the simulation.
lambda <- seq(2,6,1) #Tuning parameter lambda
Tau <- seq(0.2,0.8,0.2) #Tuning parameter tau 
n=90; n1=30; n2=30; n3=n-n1-n2 #Sample size
p1=300; p2=500; r1=280; r2=480; #Number of variables and noises in each dataset 
K=3; #Number of clusters 
mu=1; #Mean of the marginal distribution
u1=0.5; #Range of enties in the coefficient matrix

library(mvtnorm)
epsilon <- matrix(rnorm(n*(p1-r1),0,1), n, (p1-r1)) #Generation of random error in the regression model

Sigma1 <- matrix(rep(0.8,(p1-r1)^2),(p1-r1),(p1-r1)) #Generation of the covariance matrix
diag(Sigma1) <- 1

T1 <- matrix(rmvnorm(n1,mean=rep(-mu,(p1-r1)),sigma=Sigma1),n1,(p1-r1)) #Generation of the original distribution of the three clusters
T2 <- matrix(rmvnorm(n2,mean=rep(0,(p1-r1)),sigma=Sigma1),n2,(p1-r1))
T3 <- matrix(rmvnorm(n3,mean=rep(mu,(p1-r1)),sigma=Sigma1),n3,(p1-r1))

X1 <- sign(T1)*(exp(abs(T1))) #Generation of signals in X
X2 <- sign(T2)*(exp(abs(T2)))
X3 <- sign(T3)*(exp(abs(T3)))
ep1 <- (matrix(rnorm(n*r1,0,1),n,r1)) #Generation of noises in X
X <- rbind(X1,X2,X3)

beta1 <- matrix(runif((p1-r1)*(p2-r2),-u1,u1),(p1-r1),(p2-r2)) #Generation of the coefficient matrix
Z <- X%*%beta1+epsilon #Generation of signals in Z
ep2 <- (matrix(rnorm(n*r2,0.5,1),n,r2)) #Generation of noises in Z

X <- cbind(X,ep1)
Z <- cbind(Z,ep2)
#AWNCut method
Tune1 <- AWNcut.TuningSelection(X, Z, K, lambda, Tau, B=100, L=1000)
AWNcut.result <- AWNcut(X,Z,K=3, Tune1$lam, Tune1$tau,B=300,L=1000)
Cs <- AWNcut.result[[1]]$Cs
Ws <- AWNcut.result[[1]]$ws

#NCutXZ method
Tune2 <- NcutXZ.TuningSelection(X, Z, K, Tau, B=100, L=1000)
NCutXZ.result <- NcutXZ(X,Z,K=3, Tune1$tau,B=300,L=1000)
Cs <- NCutXZ.result[[1]]$Cs

#NCut method
NCut.result <- Ncut(X,K=3,B=300,L=1000)
Cs <- NCut.result$Cs

#WNCut method
WNcut.result <- WNcut(X,K=3,B=300,L=1000)
Cs <- WNcut.result$Cs
Ws <- WNcut.result$ws