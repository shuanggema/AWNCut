#' Cluster the rows of X into K clusters using the AWNCut method.
#' 
#' This function will output the clustering and feature selection result of the AWNCut method 
#' as well as the value of tuning parameters and the objective function.
#' @value lambda is the value of tuning parameter lambda for the result.
#' @value tau is the value of tuning parameter tau for the result.
#' @value Cs is the clustering result.
#' @value ws is the feature selection reuslt.
#' @value OP.value is the final value of the objective function.
#' 
#' @param X is an n x p1 matrix of n observations and p1 variables. The rows of 
#'         X into K clusters using the AWNCut method.
#' @param Z is an n x p2 matrix of n observations and p2 variables. Z is the assistant dataset.
#' @param K is the number of clusters.
#' @param lambda is a vector of tuning parameter lambda in the objective function.
#' @param Tau is a vector of tuning parameter tau in the objective function.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
#' 
#' @examples
#' set.seed(123456)
#' This sets up the initial parameters for the simulation.
#' lambda <- seq(2,6,1) #Tuning parameter lambda
#' Tau <- seq(0.2,0.8,0.2) #Tuning parameter tau
#' 
#' n=90; n1=30; n2=30; n3=n-n1-n2 #Sample size
#' p1=300; p2=500; r1=280; r2=480; #Number of variables and noises in each dataset 
#' 
#' K=3; #Number of clusters
#' 
#' mu=1; #Mean of the marginal distribution
#' u1=0.5; #Range of enties in the coefficient matrix
#'
#' library(mvtnorm)
#'  epsilon <- matrix(rnorm(n*(p1-r1),0,1), n, (p1-r1)) #Generation of random error in the regression model
#'
#'  Sigma1 <- matrix(rep(0.8,(p1-r1)^2),(p1-r1),(p1-r1)) #Generation of the covariance matrix
#'  diag(Sigma1) <- 1
#'
#'  T1 <- matrix(rmvnorm(n1,mean=rep(-mu,(p1-r1)),sigma=Sigma1),n1,(p1-r1)) #Generation of the original distribution of the three clusters
#'  T2 <- matrix(rmvnorm(n2,mean=rep(0,(p1-r1)),sigma=Sigma1),n2,(p1-r1))
#'  T3 <- matrix(rmvnorm(n3,mean=rep(mu,(p1-r1)),sigma=Sigma1),n3,(p1-r1))
#' 
#'  X1 <- sign(T1)*(exp(abs(T1))) #Generation of signals in X
#'  X2 <- sign(T2)*(exp(abs(T2)))
#'  X3 <- sign(T3)*(exp(abs(T3)))
#'  ep1 <- (matrix(rnorm(n*r1,0,1),n,r1)) #Generation of noises in X
#'  X <- rbind(X1,X2,X3)
#'
#'  beta1 <- matrix(runif((p1-r1)*(p2-r2),-u1,u1),(p1-r1),(p2-r2)) #Generation of the coefficient matrix
#'  Z <- X%*%beta1+epsilon #Generation of signals in Z
#'  ep2 <- (matrix(rnorm(n*r2,0.5,1),n,r2)) #Generation of noises in Z
#'
#'  X <- cbind(X,ep1)
#'  Z <- cbind(Z,ep2)
AWNcut <- function(X, Z, K, lambda, Tau, B=500, L=1000)  
{
  X <- scale(X)
  Z <- scale(Z)
  #Generate a Cartesian product of the two tuning parameters and try all possible conbinations
  Para <- as.data.frame(cbind(rep(lambda,each=length(Tau)),rep(Tau,length(lambda))))
  out <- list()
  for(para in 1:nrow(Para)){
  
  #Initialization
  lam <- Para[para,1] 
  tau <- Para[para,2]
  p1 <- ncol(X)
  p2 <- ncol(Z)
  w1 <- rep(1/sqrt(p1), p1)
  w2 <- rep(1/sqrt(p2), p2) 
  b <- 0
  ws.old <- c(w1,w2)
  ws <- rep(0, p1+p2)
  Cs.old <- matrix(rep(0,nrow(Z)*K),nrow(Z),K)
  for(i in 1:nrow(Z)){
    Cs.old[i,sample(K,1)] <- 1
  }
  
  while((b<=B)||(sum(ws-ws.old)/sum(ws.old)>=10e-4)){
    b <- b+1
	#Calculate the weight datasets
    wm1 <- AWNcut.W(X, Z, ws.old)
    WX1 <- wm1[[1]]
    WZ1 <- wm1[[2]]
	
	#Compute the value of the objective function using the old clustering and feature selection results
    a1 <- AWNcut.OP(X, Z, WX1, WZ1, Cs.old, tau)
    OP.value.old <- a1$TOP+lam*sum(ws.old*a1$Cor.perfeature)/(p1+p2)
	
	#Update the clustering and feature selection results
    Cs <- AWNcut.UpdateCs(WX1, WZ1, K, Cs.old)
    ws <- AWNcut.UpdateWs(X, Z, K, WX1, WZ1, b, Cs, ws.old, tau)

	#Calculate the weight datasets using updated weights
    wm2 <- AWNcut.W(X, Z, ws)
    WX2 <- wm2[[1]]
    WZ2 <- wm2[[2]]
	
	#Compute the value of the objective function using the updated clustering and feature selection results
    a2 <- AWNcut.OP(X, Z, WX2, WZ2, Cs, tau)
    OP.value <- a2$TOP+lam*sum(ws*a2$Cor.perfeature)/(p1+p2)
    if(OP.value<=OP.value.old){
       des <- rbinom(1,1,Prob(OP.value, OP.value.old, L, b))
       if(des==1){
          Cs.old <- Cs
          ws.old <- ws
       }else{
          Cs <- Cs.old
          ws <- ws.old
       }
    }else{
      Cs.old <- Cs
      ws.old <- ws
    }
  }
  out[[para]] <- list(lambda=lam, tau=tau, Cs=Cs.old, ws=ws.old, OP.value=OP.value)
  }
  return(out)
}

#' This function calculates the probability to accept the updated result 
#' when updated value of the objectove function is smaller than the old one.
#' Note that argu2 must be greater than argu1.
#' @value The output is the proobability to accept the updated result.
#'
#' @param argu1 is the updates value of the objective function
#' @param argu2 is the old value of the objective function
#' @param b is the iteration times.
Prob <- function(argu1, argu2, L, b) 
{                                    
  T <- L*log(b+1)
  return(exp(-(argu2-argu1)/T))
}

#' This function calculates the weighed datasets for the AWNCut method.
#' @value WX is the weighed dataset X
#' @value WZ is the weighed dataset Z
#' 
#' @param X is a n x p1 matrix
#' @param Z is a n x p2 matrix
#' @param ws is a vector of weights for both X and Z datasets.
AWNcut.W <- function(X, Z, ws){
  p1 <- ncol(X)
  p2 <- ncol(Z)
  #Calculate matrix W
  DistX <- as.matrix(dist(t(t(X)*ws[1:p1]), upper=T, diag=T))
  diag(DistX) <- 1
  WX <- DistX^(-1)
  DistZ <- as.matrix(dist(t(t(Z)*ws[(1+p1):(p1+p2)]), upper=T, diag=T))
  diag(DistZ) <- 1
  WZ <- DistZ^(-1)
  return(list(WX, WZ))
}

#' This function updates the clustering result of the AWNCut method.
#' @value Cs is the updated clustering result
#'
#' @param WX is the weighed dataset X.
#' @param WZ is the weighed dataset Z.
#' @param K is the number of clusters.
#' @param Cs is the clustering result from the las iteration.
AWNcut.UpdateCs <- function(WX, WZ, K, Cs)
{
  P <- NULL
  Cs.rec <- NULL
  for(i in 1:ncol(Cs)){
    # Calculate cutvol/cut for each dataset and record the pair 
	# that has the largest Euclidean distance in each dataset.
    T1 <- WX[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WX[which(Cs[,i]==1),which(Cs[,i]==0)]
    T3 <- WZ[which(Cs[,i]==1),which(Cs[,i]==1)]
    T4 <- WZ[which(Cs[,i]==1),which(Cs[,i]==0)]
    P <- c(P, sum(T1[upper.tri(T1)])/sum(T2)+sum(T3[upper.tri(T3)])/sum(T4))
    Cs.rec <- rbind(Cs.rec,cbind(which(T1==min(T1),arr.ind=T)[1,],which(T3==min(T3),arr.ind=T)[1,]))
  }
  P <- P/sum(P)
   #Select K(+) and K(-) according to cutvol/cut, K(+) is the first element in S and K(-) is the last element in S
  S <- sample(K,K,replace=F,prob=P)
   #Select the observation that needs to be rearranged and record its location in the clustering result of the last iteration
  Ob.rec <- sample(Cs.rec[S[1],],1)
  ob.rec <- which(Cs[,S[1]]==1)[Ob.rec]
  
  Cs.new <- Cs
  Cs.new[ob.rec,S[1]] <- 0
  Cs.new[ob.rec,S[K]] <- 1
  #To make sure each cluster has at least two observations
  t <- sum(apply(Cs.new,2,max)) 
  l <- min(apply(Cs.new,2,function(x){return(length(which(x==1)))}))
  if((t==K)&&(l>=2)){
    return(Cs.new)
  }else{
    return(Cs)
  }
}

#' This function updates the feature selection result of the AWNCut method.
#' @value The output is a vector of the standardized updated weights. 
#' 
#' @param X is an n x p1 matrix.
#' @param Z is an n x p2 matrix.
#' @param K is the number of clusters.
#' @param WX is the weighed dataset X.
#' @param WZ is the weighed dataset Z.
#' @param b is the iteration times.
#' @param Cs is the old clustering result.
#' @param ws is the old feature selection result.
#' @param tau is a vector of tuning parameter tau in the objective function.
AWNcut.UpdateWs <- function(X, Z, K, WX, WZ, b, Cs, ws, tau)
{
  n <- nrow(X)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  #Compute the pace of each iteration according to b
  pace.x <- 1/(10*ceiling(b/10)*sqrt(p1))
  pace.z <- 1/(10*ceiling(b/10)*sqrt(p2))
  pace1 <- rep(-1, p1)
  pace2 <- rep(-1, p2)
  #Calculate the average correlation of each feature
  temp <- AWNcut.OP(X, Z, WX, WZ, Cs, tau)$Cor.perfeature
  #Use Kmeans method to cluster the average weight into 2 clusters in each dataset
  temp1 <- kmeans(temp[1:p1],2)
  temp2 <- kmeans(temp[(1+p1):(p1+p2)],2)
  #For features in the cluster that has a higher level of average correlation, update their weight by adding the pace
  pace1[which(temp1$cluster==which.max(temp1$center))] <- 1
  pace2[which(temp2$cluster==which.max(temp2$center))] <- 1
  pace <- c(pace1,pace2)
  ws.new <- ws+c(pace*c(rep(pace.x,p1),rep(pace.z,p2)))
  #To make sure that each weight is nonngative
  for(i in 1:length(ws.new)){
    ws.new[i] <- max(ws.new[i],0)
  } 
  return(c(ws.new[1:p1]/l2n(ws.new[1:p1]),ws.new[(p1+1):(p1+p2)]/l2n(ws.new[(p1+1):(p1+p2)])))
}

#' This function calculate the value of the objective function.
#' @value OP1 is the value of the NCut measure in X.
#' @value OP2 is the value of the NCut measure in Z.
#' @value TOP is the sum of the NCut measure in both X and Z.
#' @value Cor.X is a vector of the average correlation for X.
#' @value Cor.Z is a vector of the average correlation for Z.
#' @value Cor.perfeature is a combination of the average correlation for X and Z.
#'
#' @param X is a n x p1 matrix.
#' @param Z is a n x p2 matrix.
#' @param WX is the weighed dataset X.
#' @param WZ is the weighed dataset Z.
#' @param Cs is clustering result.
#' @param tau is tuning parameter in the objective function.
AWNcut.OP <- function(X, Z, WX, WZ, Cs, tau)
{
  n <- nrow(X)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  #Calculate cuts and cutvols for each dataset
  cutvolX <- NULL
  cutX <- NULL
  for(i in 1:ncol(Cs)){
    T1 <- WX[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WX[which(Cs[,i]==1),which(Cs[,1]==0)]
    cutvolX <- c(cutvolX, sum(T1[upper.tri(T1)]))
    cutX <- c(cutX,sum(T2))
  }
  OP1 <- cutvolX/cutX
  
  cutvolZ <- NULL
  cutZ <- NULL
  for(i in 1:ncol(Cs)){
    T1 <- WZ[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WZ[which(Cs[,i]==1),which(Cs[,1]==0)]
    cutvolZ <- c(cutvolZ, sum(T1[upper.tri(T1)]))
    cutZ <- c(cutZ,sum(T2))
  }
  OP2 <- cutvolZ/cutZ
  #Calculate average correlations per feature
  Cor.X <- 0
  Cor.Z <- 0
  for(i in 1:ncol(Cs)){
    T <- cor(X[which(Cs[,i]==1),],Z[which(Cs[,i]==1),])
    T[is.na(T)] <- 0
    Cor.X <- Cor.X+apply(abs(T),1,mean)
    Cor.Z <- Cor.Z+apply(abs(T),2,mean)
  }
  return(list(OP1=sum(OP1), OP2=sum(OP2), TOP=sum(OP1)+tau*sum(OP2), Cor.X=Cor.X, Cor.Z=Cor.Z, Cor.perfeature=c(Cor.X,Cor.Z)))
}

#' This function standardize a vector
#' @value The output is a standardized vector.
#'
#' @param vec is a vector.
l2n <- function(vec){
  return(sqrt(sum(vec^2)))
}

#' This fuction output the clustering and feature selection result of the WNCut method.
#' @value Cs is the clustering result.
#' @value ws is the feature selection reuslt.
#' @value OP.value is the final value of the objective function.
#' 
#' @param X is a n x p matrix with n observations and p variables.
#' @param K is the number of clusters.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
WNcut <- function(X, K, B=500, L=1000)
{
  X <- scale(X)
  out <- list()
  
  #Initialization
  p1 <- ncol(X)
  b <- 0
  ws.old <- rep(1/sqrt(p1), p1)
  ws <- rep(0, p1)
  Cs.old <- matrix(rep(0,nrow(X)*K),nrow(X),K)
  for(i in 1:nrow(X)){
    Cs.old[i,sample(K,1)] <- 1
  }
  while((b<=B)||(sum(ws-ws.old)/sum(ws.old)>=10e-4)){
    b <- b+1
	
	#Calculate the weighed dataset
    WX1 <- WNcut.W(X, ws.old)
	
	#Compute the value of the objective function using the old clustering and feature selection results
    OP.value.old <- WNcut.OP(X, WX1, Cs.old)
	
	#Update the clustering result
    Cs <- WNcut.UpdateCs(WX1,K,Cs.old)
	
	#Update the feature selection result
    ws<- WNcut.UpdateWs(X, b, ws.old)
	
	#Calculate the weighed dataset using the updated weights
    WX2 <- WNcut.W(X, ws)
	
	#Compute the value of the objective function using the updated clustering and feature selection results
    OP.value <- WNcut.OP(X, WX2, Cs)
    if(OP.value<=OP.value.old){
       des <- rbinom(1,1,Prob(OP.value, OP.value.old, L, b))
       if(des==1){
          Cs.old <- Cs
          ws.old <- ws
       }else{
          Cs <- Cs.old
          ws <- ws.old
       }
    }else{
      Cs.old <- Cs
      ws.old <- ws
    }
  }
  out <- list(Cs=Cs.old, ws=ws.old, OP.value=OP.value)
  return(out)
}

#' This function calculate the weight dataset for the WNCut method.
#' @value WX is the weighed dataset.
#'
#' @param X is a n x p1 matrix.
#' @param ws is a vector of weights.
WNcut.W <- function(X, ws){
  p1 <- ncol(X)
  DistX <- as.matrix(dist(t(t(X)*ws[1:p1]), upper=T, diag=T))
  diag(DistX) <- 1
  WX <- DistX^(-1)
  return(WX)
}

#' This function updates the clustering result for the WNCut method.
#' @value Cs is the updated clustering result.
#'
#' @param WX is a weighed dataset X.
#' @param K is the number of clusters.
#' @param Cs is the old clustering result.
WNcut.UpdateCs <- function(WX, K, Cs)
{
  P <- NULL
  Cs.rec <- NULL
  for(i in 1:ncol(Cs)){
    TX <- WX[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WX[which(Cs[,i]==1),which(Cs[,i]==0)]
    P <- c(P, sum(TX[upper.tri(TX)])/sum(T2))
    Cs.rec <- rbind(Cs.rec,which(TX==min(TX),arr.ind=T)[1,])
  }
  P <- P/l2n(P)
  S <- sample(K,K,replace=F,prob=P) 
  Ob.rec <- sample(Cs.rec[S[1],],1)
  ob.rec <- which(Cs[,S[1]]==1)[Ob.rec]
  Cs.new <- Cs
  Cs.new[ob.rec,S[1]] <- 0
  Cs.new[ob.rec,S[K]] <- 1
  t <- sum(apply(Cs.new,2,max)) #To make sure each cluster has at least one observation.
  l <- min(apply(Cs.new,2,function(x){return(length(which(x==1)))}))
  if((t==K)&&(l>=2)){
    return(Cs.new)
  }else{
    return(Cs)
  }
}

#' This function updates the weights for the WNCut method.
#' @value The output is a vector of the standardized updated weights.
#'
#' @param X is a n x p matrix.
#' @param b is the iteration times.
#' @param ws is the old feature selection reuslt.
WNcut.UpdateWs <- function(X, b, ws)
{
  n <- nrow(X)
  p1 <- ncol(X)
  pace.x <- 1/(10*ceiling(b/10)*sqrt(p1))
  pace <- NULL
  for(i in 1:p1){
   #Without the assistant dataset, there is no way to calculate average correlation.
   #So whether the weight of each feature should add or minus the pace is decided randomly and 
   #the feature selection result is judged by the comparison of the value of objective function
    pace <- c(pace,rbinom(1,1,0.5))
  }
  pace[pace==0] <- -1
  ws.new <- ws+pace*pace.x
  for(i in 1:length(ws.new)){
    ws.new[i] <- max(ws.new[i],0)
  } 
  return(ws.new/l2n(ws.new))
}

#' This function calculates the value of the objective function of the WNCut method.
#' @value OPX is the value of the NCut measure in X. 
#'            It's also the value of the objective function.
#'
#' @param X is a n x p matrix.
#' @param WX is the weighed dataset X.
#' @param Cs is the clustering result.
WNcut.OP <- function(X, WX, Cs)
{
  n <- nrow(X)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  #Calculate cut and cutvol
  cutvolX <- NULL
  cutX <- NULL
  for(i in 1:ncol(Cs)){
    T1 <- WX[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WX[which(Cs[,i]==1),which(Cs[,1]==0)]
    cutvolX <- c(cutvolX, sum(T1[upper.tri(T1)]))
    cutX <- c(cutX,sum(T2))
  }
  OPX <- cutvolX/cutX
  return(sum(OPX))
}

#' This function output the clustering result of the NCutXZ method, 
#' which use both datasets for clustering without feature selection.
#' @value tau is the value of tuning parameter tau for the result.
#' @value Cs is the clustering result.
#' @value OP.value is the final value of the objective function.
#' 
#' @param X is an n x p1 matrix of n observations and p1 variables. The rows of 
#'          X into K clusters using the NCutXZ method.
#' @param Z is an n x p2 matrix of n observations and p2 variables. Z is the assistant dataset.
#' @param K is the number of clusters.
#' @param Tau is a vector of tuning parameter tau which decides the 
#'            weight of each dataset in the compution of the objective function.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
NcutXZ <- function(X, Z, K, Tau, B=500, L=1000)  
{
  X <- scale(X)
  Z <- scale(Z)

  out <- list()
  for(tau in Tau){
  
  #Initialization
  p1 <- ncol(X)
  p2 <- ncol(Z)
  w1 <- rep(1/sqrt(p1), p1)
  w2 <- rep(1/sqrt(p2), p2) 
  b <- 0
  ws <- c(w1,w2)
  Cs.old <- matrix(rep(0,nrow(Z)*K),nrow(Z),K)
  for(i in 1:nrow(Z)){
    Cs.old[i,sample(K,1)] <- 1
  }
  
  #Calculate the weighed dataset using the uniform weights
  wm1 <- AWNcut.W(X, Z, ws)
  WX1 <- wm1[[1]]
  WZ1 <- wm1[[2]]
  while((b<=B)){
    b <- b+1
	
	#Calculate the value of the objective function using the old clustering result
    a1 <- NcutXZ.OP(X, Z, WX1, WZ1, Cs.old, tau)
    OP.value.old <- a1$TOP
	
	#Update the clustering result
    Cs <- AWNcut.UpdateCs(WX1, WZ1, K, Cs.old)
	
	#Calculate the value of the objective function using updated clustering result	
    a2 <- NcutXZ.OP(X, Z, WX1, WZ1, Cs, tau)
    OP.value <- a2$TOP
    if(OP.value<=OP.value.old){
       des <- rbinom(1,1,Prob(OP.value, OP.value.old, L, b))
       if(des==1){
          Cs.old <- Cs
       }else{
          Cs <- Cs.old
       }
    }else{
      Cs.old <- Cs
    }
  }
  out[[which(Tau==tau)]] <- list(tau=tau, Cs=Cs.old, OP.value=OP.value)
  }
  return(out)
}

#' This function calculates the value of the objective function for the NCutXZ method.
#' @value OP1 is the value of the NCut measure in X.
#' @value OP2 is the value of the NCut measure in Z.
#' @value TOP is the sum of the NCut measure in both X and Z.
#'
#' @param X is a n x p1 matrix.
#' @param Z is a n x p2 matrix.
#' @param WX is the weighed dataset X.
#' @param WZ is the weighed dataset Z.
#' @param Cs is clustering result.
#' @param tau is tuning parameter in the objective function.
NcutXZ.OP <- function(X, Z, WX, WZ, Cs, tau)
{
  n <- nrow(X)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  #Calculate cut and cutvol for both datasets
  cutvolX <- NULL
  cutX <- NULL
  for(i in 1:ncol(Cs)){
    T1 <- WX[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WX[which(Cs[,i]==1),which(Cs[,1]==0)]
    cutvolX <- c(cutvolX, sum(T1[upper.tri(T1)]))
    cutX <- c(cutX,sum(T2))
  }
  OP1 <- cutvolX/cutX
  
  cutvolZ <- NULL
  cutZ <- NULL
  for(i in 1:ncol(Cs)){
    T1 <- WZ[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WZ[which(Cs[,i]==1),which(Cs[,1]==0)]
    cutvolZ <- c(cutvolZ, sum(T1[upper.tri(T1)]))
    cutZ <- c(cutZ,sum(T2))
  }
  OP2 <- cutvolZ/cutZ
  return(list(OP1=sum(OP1), OP2=sum(OP2), TOP=sum(OP1)+tau*sum(OP2)))
}

#' This function outputs the clustering result of the NCut method.
#' @value Cs is the clustering result.
#' @value OP.value is the final value of the objective function.
#' 
#' @param X is an n x p1 matrix of n observations and p1 variables. The rows of 
#'          X into K clusters using the AWNCut method.
#' @param K is the number of clusters.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
Ncut <- function(X, K, B=500, L=1000)
{

  X <- scale(X)
  out <- list()

  #Initialization
  p1 <- ncol(X)
  b <- 0
  Cs.old <- matrix(rep(0,nrow(X)*K),nrow(X),K)
  for(i in 1:nrow(X)){
    Cs.old[i,sample(K,1)] <- 1
  }
  
  #Calculate the weighed dataset using the uniform weights
  WX1 <- WNcut.W(X, rep(1/sqrt(p1),p1))
  while(b<=B){
    b <- b+1
	
	#Calculate the value of the objective function using the old clustering result
    OP.value.old <- Ncut.OP(X, WX1, Cs.old)
	
	#Update the clustering result
    Cs <- Ncut.UpdateCs(WX1, K, Cs.old)
	
	#Calculate the value of the objective function using the updated clustering result
    OP.value <- Ncut.OP(X, WX1, Cs)
    if(OP.value<=OP.value.old){
       des <- rbinom(1,1,Prob(OP.value, OP.value.old, L, b))
       if(des==1){
          Cs.old <- Cs
       }else{
          Cs <- Cs.old
       }
    }else{
      Cs.old <- Cs
    }
  }
  out <- list(Cs=Cs.old, OP.value=OP.value)
  return(out)
}

#' This function updates the clustering result of the NCut method.
#' @value Cs is the updated clustering result.
#'
#' @param WX is a weighed dataset X.
#' @param K is the number of clusters.
#' @param Cs is the old clustering result.
Ncut.UpdateCs <- function(WX, K, Cs)
{
  P <- NULL
  Cs.rec <- NULL
  for(i in 1:ncol(Cs)){
    TX <- WX[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WX[which(Cs[,i]==1),which(Cs[,i]==0)]
    P <- c(P, sum(TX[upper.tri(TX)])/sum(T2))
    Cs.rec <- rbind(Cs.rec,which(TX==min(TX),arr.ind=T)[1,])
  }
  P <- P/l2n(P)
  S <- sample(K,K,replace=F,prob=P)
  Ob.rec <- sample(Cs.rec[S[1],],1)
  ob.rec <- which(Cs[,S[1]]==1)[Ob.rec]
  Cs.new <- Cs
  Cs.new[ob.rec,S[1]] <- 0
  Cs.new[ob.rec,S[K]] <- 1
  t <- sum(apply(Cs.new,2,max)) #To make sure each cluster has at least one observation.
  l <- min(apply(Cs.new,2,function(x){return(length(which(x==1)))}))
  if((t==K)&&(l>=2)){
    return(Cs.new)
  }else{
    return(Cs)
  }
}

#' This function calculates the value of the objective function for the NCut method.
#' @value OPX is the value of the NCut measure in X. 
#'            It's also the value of the objective function.
#'
#' @param X is a n x p matrix.
#' @param WX is the weighed dataset X.
#' @param Cs is the clustering result.
Ncut.OP <- function(X, WX, Cs)
{
  n <- nrow(X)
  p1 <- ncol(X)
  p2 <- ncol(Z)
  #Calculate cut and cutvol
  cutvolX <- NULL
  cutX <- NULL
  for(i in 1:ncol(Cs)){
    T1 <- WX[which(Cs[,i]==1),which(Cs[,i]==1)]
    T2 <- WX[which(Cs[,i]==1),which(Cs[,1]==0)]
    cutvolX <- c(cutvolX, sum(T1[upper.tri(T1)]))
    cutX <- c(cutX,sum(T2))
  }
  OPX <- cutvolX/cutX
  return(sum(OPX))
}

#' This function outputs the value of DBI for a clustering result.
#' @value The output is the value of DBI for a clustering result.
#' 
#' @param X is a n x p matrix with n observations and p variables.
#' @param K is the number of clusters.
#' @param Cs is a n x K matrix containing the clustering result of X. The entries
#'        in Cs is either 1 or 0. If the ith observation is in the kth cluster, then Cs(i,k)=1,
#'        otherwise, it equals to 0.
#' @param ws is a vector of the weights for the p variables. the length of ws equals p.
DBI <- function(X, K, Cs, ws){
  X <- scale(X)
  p1 <- ncol(X)
  w1 <- ws[1:p1]/l2n(ws[1:p1]) 
  WX <- as.matrix(dist(t(t(X)*w1), upper=T, diag=T))
  
  #Calculate both cut and cuvol, cutvol is the diagonal of the matrix Cut, cut is the reat of the entries
  Cut <- matrix(0,K,K)
  for(i in 1:K){
    for(j in 1:K){
      T <- WX[which(Cs[,i]==1),which(Cs[,j]==1)]
      Cut[i,j] <- sum(T)
    }
  }
  cutvol <- diag(Cut)/2
  DBI <- NULL
  for(i in 1:K){
    t <- NULL
    for(j in (1:K)[-i]){
      t <- c(t,Cut[i,j]/(cutvol[i]+cutvol[j]))
    }
    DBI <- c(DBI, max(t))
  }
  return(mean(DBI[DBI<Inf])) #When a clustering result contains two clusters that has only one observation, the value of DBI will be Inf.
}

#' This function outputs the selection of tuning parameters for the AWNCut method.
#' @value num is the position of the max DBI.
#' @value Table is the Table of the DBI for all possible combination of the parameters.
#' @value lam is the best choice of tuning parameter lambda.
#' @value tau is the best choice of tuning parameter tau.
#'
#' @param X is an n x p1 matrix of n observations and p1 variables. The rows of 
#'             X into K clusters using the AWNCut method.
#' @param Z is an n x p2 matrix of n observations and p2 variables. Z is the assistant dataset.
#' @param K is the number of clusters.
#' @param lambda is a vector of tuning parameter lambda in the objective function.
#' @param Tau is a vector of tuning parameter tau in the objective function.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
AWNcut.TuningSelection <- function(X, Z, K, lambda, Tau, B=500, L=1000){
  out <- AWNcut(X, Z, K, lambda, Tau, B, L=1000)
  Para <- as.data.frame(cbind(rep(lambda,each=length(Tau)),rep(Tau,length(lambda))))
  dbi <- NULL
  for(i in 1:nrow(Para)){
    Cs <- out[[i]]$Cs
    ws <- out[[i]]$ws
    dbi <- c(dbi, DBI(cbind(X,Z),K,Cs,ws))
  }
  return(list(num=which.max(dbi),Table=t(cbind(Para,dbi)), lam=Para[which.max(dbi),1], tau=Para[which.max(dbi),2], DBI=max(dbi)))
}

#' This function outputs the selection of tuning parameters for the NCutXZ method.
#' @value num is the position of the max DBI.
#' @value Table is the Table of the DBI for all possible parameters.
#' @value tau is the best choice of tuning parameter tau.
#'
#' @param X is an n x p1 matrix of n observations and p1 variables. The rows of 
#'             X into K clusters using the NCutXZ method.
#' @param Z is an n x p2 matrix of n observations and p2 variables. Z is the assistant dataset.
#' @param K is the number of clusters.
#' @param Tau is a vector of tuning parameter tau in the objective function.
#' @param B is the number of iterations in the simulated annealing algorithm.
#' @param L is the temperature coefficient in the simulated annealing algorithm.
NcutXZ.TuningSelection <- function(X, Z, K, Tau, B=500, L=1000){
  out <- NcutXZ(X, Z, K, Tau, B, L=1000)
  dbi <- NULL
  for(i in 1:length(Tau)){
    Cs <- out[[i]]$Cs
    ws <- rep(1,(ncol(X)+ncol(Z)))
    dbi <- c(dbi, DBI(cbind(X,Z),K,Cs,ws))
  }
  return(list(num=which.max(dbi),Table=rbind(Tau, dbi), tau=Tau[which.max(dbi)], DBI=max(dbi)))
}

#' This function calculates the true error rate of a clustering result.
#' @value err is the true error rate of a clustering result.
#'
#' @param X is a clustering result in matrix format.
ErrorRate <- function(X)
{
  n=nrow(X)
  Error <- matrix(1,n,n)
  Error[1:n1,1:n1]<-0
  Error[(1+n1):(n1+n2),(1+n1):(n1+n2)]<-0
  Error[(1+n2+n1):n,(1+n2+n1):n] <- 0
  Denum <- sum(Error)
  err <- 0
  for(i in 1:ncol(X)){
    f <- matrix(X[,i],n,1)
    err <- err+sum(f%*%t(f)*Error)/Denum
  }
  return(err)
}

#' This function transfers a clustering result into a matrix format from a vector.
#' @value Kcs is a clustering result in the matrix format.
#' 
#' @param x is a clustering result in the vector format.
Kcs <- function(x)
{
  n <- length(x)
  K <- length(unique(x))
  Kcs <- matrix(0,n,K)
  for(i in 1:K){
    Kcs[which(x==i),i] <- 1
  }
  return(Kcs)
}

#' This function calculates the stablility of the simulation result
#' @value The output if the stability of the simulation result
#'
#' @param x is a 3 dimensional array. 
#'          The first dimension equals to sample size.
#'          The second dimension equals to number og clusters.
#'          The third dimension equals to the replication times. 
Stability <- function(x){
  N <- dim(x)[3]
  n <- dim(x)[1]
  sta <- 0
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      sta <- sta+sum(abs(x[,,i]-x[,,j]))/(n^2)
    }
  } 
  return(sta/choose(N,2))
}


