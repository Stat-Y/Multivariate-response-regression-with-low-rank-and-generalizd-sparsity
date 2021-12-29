# Data
load("Example_Data_Z_neq_I.Rdata")

# Setting
library(matrixcalc)
library(progress)
library(Matrix)
library(caret)
library(e1071)
library(dplyr)
library(lattice)
library(MASS)
library(foreach)
library(doParallel)
library(gplots)

# Sylvester Equation Function
sylvester=function(A,B,C){
  X=matrix(NA,dim(C)[1],dim(C)[2])
  XX=matrix(NA,dim(C)[1],dim(C)[2])
  TA=Schur(A)$T
  TB=Schur(B)$T
  UA=Schur(A)$Q
  UB=Schur(B)$Q
  CC=t(UA)%*%C%*%UB
  XX[,1]=solve(TA+TB[1,1]*diag(dim(C)[1]))%*%CC[,1]
  XX[,2]=solve(TA+TB[2,2]*diag(dim(C)[1]))%*%(CC[,2]-XX[,1:2-1]*TB[,2][1:2-1])
  if(dim(C)[2]>2){
    for(i in 3:dim(C)[2])
    {
      XX[,i]=solve(TA+TB[i,i]*diag(dim(C)[1]))%*%(CC[,i]-XX[,1:i-1]%*%TB[,i][1:i-1])
      i=i+1
    }
  }
  X=UA%*%XX%*%t(UB)
  return(X)
}

# ADMM Algorithm for Z=I
ADMM_3block <- function(initial_Beta_0,initial,X,Y,lambda_1,lambda_2,rho,maxiter,tole){
  pb <-  progress_bar$new(total = maxiter-1)
  abs <- 0.1^tole
  rel <- 0.1^tole
  b_0 <- matrix(0,dim(Y)[2],maxiter)
  B <-  array(0,c(dim(X)[2],dim(Y)[2],maxiter))
  Q <-  array(0,c(dim(X)[2],dim(Y)[2],maxiter))
  U <- array(0,c(dim(X)[2],dim(Y)[2],maxiter))
  S <- array(0,c(dim(X)[2],dim(Y)[2],maxiter))
  V <- array(0,c(dim(X)[2],dim(Y)[2],maxiter))
  b_0[,1] <- initial_Beta_0
  B[,,1] <- initial
  Q[,,1] <- initial
  U[,,1] <- initial
  S[,,1] <- initial
  V[,,1] <- initial
  error <- data.frame(matrix(0,maxiter,4))
  error_b_0 <- rep(0,dim(Y)[2])
  stopcriteria <- data.frame(matrix(0,maxiter,4))
  error[1,] <- Inf
  error_b_0[1] <- Inf
  names(error) <- c("B-Q","B-S","Q","S")
  i <- 2
  while ( (i>=2) && ( (i<=maxiter) && ((error[i-1,1] > sqrt(dim(X)[2]*dim(Y)[2])*abs + max( sqrt(sum(B[,,i-1]^2)), sqrt(sum(Q[,,i-1]^2)) )*rel) ||
                                       (error[i-1,2] > sqrt(dim(X)[2]*dim(Y)[2])*abs + max( sqrt(sum(B[,,i-1]^2)), sqrt(sum(S[,,i-1]^2)) )*rel) ||
                                       (rho*error[i-1,3] > sqrt(dim(X)[2]*dim(Y)[2])*abs + sqrt(sum(U[,,i-1]^2))*rel) ||
                                       (rho*error[i-1,4] > sqrt(dim(X)[2]*dim(Y)[2])*abs + sqrt(sum(V[,,i-1]^2))*rel) ||
                                       (rho*error_b_0[i-1] > sqrt(dim(Y)[2])*abs)   ) ) )  
  {
    pb$tick()        
    B[,,i] <- solve(t(X)%*%X/(dim(X)[1])+2*rho*diag(dim(X)[2]))%*%(t(X)%*%Y/(dim(X)[1])-t(X)%*%rep(1,dim(X)[1])%*%t(b_0[,i-1])/(dim(X)[1])-U[,,i-1]+rho*Q[,,i-1]-V[,,i-1]+rho*S[,,i-1])
    b_0[,i] <- (t(Y)%*%rep(1,dim(X)[1])-t(B[,,i])%*%t(X)%*%rep(1,dim(X)[1]))/(dim(X)[1])
    Q[,,i] <- svd(B[,,i]+(1/rho)*U[,,i-1])$u %*% pmax(diag(svd(B[,,i]+(1/rho)*U[,,i-1])$d)-(lambda_1)/rho,0) %*% t(svd(B[,,i]+(1/rho)*U[,,i-1])$v)
    S[,,i] <- hadamard.prod(sign(B[,,i]+(1/rho)*V[,,i-1]), pmax(abs(B[,,i]+(1/rho)*V[,,i-1])-(lambda_2)/rho,0))
    U[,,i] <- U[,,i-1]+rho*(B[,,i]-Q[,,i])
    V[,,i] <- V[,,i-1]+rho*(B[,,i]-S[,,i])
    error[i,1] <- sqrt(sum((B[,,i]-Q[,,i])^2))
    error[i,2] <- sqrt(sum((B[,,i]-S[,,i])^2))
    error[i,3] <- sqrt(sum((Q[,,i]-Q[,,i-1])^2))
    error[i,4] <- sqrt(sum((S[,,i]-S[,,i-1])^2))
    error_b_0[i] <- sqrt(sum((b_0[,i]-b_0[,i-1])^2))
    
    stopcriteria[i,1] <- sqrt(dim(X)[2]*dim(Y)[2])*abs + max( sqrt(sum(B[,,i]^2)), sqrt(sum(Q[,,i]^2)) )*rel
    stopcriteria[i,2] <- sqrt(dim(X)[2]*dim(Y)[2])*abs + max( sqrt(sum(B[,,i]^2)), sqrt(sum(S[,,i]^2)) )*rel
    stopcriteria[i,3] <- sqrt(dim(X)[2]*dim(Y)[2])*abs + sqrt(sum(U[,,i]^2))*rel
    stopcriteria[i,4] <- sqrt(dim(X)[2]*dim(Y)[2])*abs + sqrt(sum(V[,,i]^2))*rel
    
    errorZ_B_Q<<- c(error[i,1],stopcriteria[i,1],stopcriteria[i,1]/error[i,1] )
    errorZ_BZ_S<<-c(error[i,2],stopcriteria[i,2],stopcriteria[i,2]/error[i,2] )
    errorZ_Q_Q<<-c(rho*error[i,3],stopcriteria[i,3],stopcriteria[i,3]/(rho*error[i,3]) )
    errorZ_S_S<<-c(rho*error[i,4],stopcriteria[i,4],stopcriteria[i,4]/(rho*error[i,4]) )
    errorZ_b_0_b_0<<-c(rho*error_b_0[i],sqrt(dim(Y)[2])*abs,(sqrt(dim(Y)[2])*abs)/(rho*error_b_0[i]) )
    iteration<<-i
    i <- i+1
  }
  return(list(b_0[,i-1],B[,,i-1],Q[,,i-1],S[,,i-1]))
}

# ADMM Algorithm for Z!=I
ADMM_3block_Z <- function(initial_Beta_0,initial,X,Y,Z,lambda_1,lambda_2,rho,maxiter,tole){
  pb <- progress_bar$new(total = maxiter-1)
  abs <- 0.1^tole
  rel <- 0.1^tole
  b_0 <- matrix(0,dim(Y)[2],maxiter)
  B <- array(0,c(dim(X)[2],dim(Y)[2],maxiter))
  Q <- array(0,c(dim(X)[2],dim(Y)[2],maxiter))
  U <- array(0,c(dim(X)[2],dim(Y)[2],maxiter))
  S <- array(0,c(dim(X)[2],dim(Z)[2],maxiter))
  V <- array(0,c(dim(X)[2],dim(Z)[2],maxiter))
  b_0[,1] <- initial_Beta_0
  B[,,1] <- initial
  Q[,,1] <- initial
  U[,,1] <- initial
  S[,,1] <- initial%*%Z
  V[,,1] <- initial%*%Z
  error <- data.frame(matrix(0,maxiter,4))
  error_b_0 <- rep(0,dim(Y)[2])
  stopcriteria <- data.frame(matrix(0,maxiter,4))
  error[1,] <- Inf
  error_b_0[1] <- Inf
  names(error) <- c("B-Q","BZ-S","Q","S")
  i <- 2
  while ( (i>=2) && ( (i<=maxiter) && ((error[i-1,1] > sqrt(dim(X)[2]*dim(Y)[2])*abs + max( sqrt(sum(B[,,i-1]^2)), sqrt(sum(Q[,,i-1]^2)) )*rel) ||
                                       (error[i-1,2] > sqrt(dim(X)[2]*dim(Z)[2])*abs + max( sqrt(sum((B[,,i-1]%*%Z)^2)), sqrt(sum(S[,,i-1]^2)) )*rel) ||
                                       (rho*error[i-1,3] > sqrt(dim(X)[2]*dim(Y)[2])*abs + sqrt(sum(U[,,i-1]^2))*rel) ||
                                       (rho*error[i-1,4] > sqrt(dim(X)[2]*dim(Z)[2])*abs + sqrt(sum(V[,,i-1]^2))*rel) ||
                                       (rho*error_b_0[i-1] > sqrt(dim(Y)[2])*abs)  ) ) )  
  {
    pb$tick()        
    B[,,i] <- sylvester(t(X)%*%X/(dim(X)[1])+rho*diag(dim(X)[2]),rho*Z%*%t(Z),t(X)%*%Y/(dim(X)[1])-t(X)%*%rep(1,dim(X)[1])%*%t(b_0[,i-1])/(dim(X)[1])-U[,,i-1]+rho*Q[,,i-1]-V[,,i-1]%*%t(Z)+rho*S[,,i-1]%*%t(Z))
    b_0[,i] <- (t(Y)%*%rep(1,dim(X)[1])-t(B[,,i])%*%t(X)%*%rep(1,dim(X)[1]))/(dim(X)[1])
    Q[,,i] <- svd(B[,,i]+(1/rho)*U[,,i-1])$u %*% pmax(diag(svd(B[,,i]+(1/rho)*U[,,i-1])$d)-(lambda_1)/rho,0) %*% t(svd(B[,,i]+(1/rho)*U[,,i-1])$v)
    S[,,i] <- hadamard.prod(sign(B[,,i]%*%Z+(1/rho)*V[,,i-1]), pmax(abs(B[,,i]%*%Z+(1/rho)*V[,,i-1])-(lambda_2)/rho,0))
    U[,,i] <- U[,,i-1]+rho*(B[,,i]-Q[,,i])
    V[,,i] <- V[,,i-1]+rho*(B[,,i]%*%Z-S[,,i])
    error[i,1] <- sqrt(sum((B[,,i]-Q[,,i])^2))
    error[i,2] <- sqrt(sum((B[,,i]%*%Z-S[,,i])^2))
    error[i,3] <- sqrt(sum((Q[,,i]-Q[,,i-1])^2))
    error[i,4] <- sqrt(sum((S[,,i]-S[,,i-1])^2))
    error_b_0[i] <- sqrt(sum((b_0[,i]-b_0[,i-1])^2))
    
    stopcriteria[i,1] <- sqrt(dim(X)[2]*dim(Y)[2])*abs + max( sqrt(sum(B[,,i]^2)), sqrt(sum(Q[,,i]^2)) )*rel
    stopcriteria[i,2] <- sqrt(dim(X)[2]*dim(Z)[2])*abs + max( sqrt(sum((B[,,i]%*%Z)^2)), sqrt(sum(S[,,i]^2)) )*rel
    stopcriteria[i,3] <- sqrt(dim(X)[2]*dim(Y)[2])*abs + sqrt(sum(U[,,i]^2))*rel
    stopcriteria[i,4] <- sqrt(dim(X)[2]*dim(Z)[2])*abs + sqrt(sum(V[,,i]^2))*rel
    
    errorZ_B_Q<<- c(error[i,1],stopcriteria[i,1],stopcriteria[i,1]/error[i,1] )
    errorZ_BZ_S<<-c(error[i,2],stopcriteria[i,2],stopcriteria[i,2]/error[i,2] )
    errorZ_Q_Q<<-c(rho*error[i,3],stopcriteria[i,3],stopcriteria[i,3]/(rho*error[i,3]) )
    errorZ_S_S<<-c(rho*error[i,4],stopcriteria[i,4],stopcriteria[i,4]/(rho*error[i,4]) )
    errorZ_b_0_b_0<<-c(rho*error_b_0[i],sqrt(dim(Y)[2])*abs,(sqrt(dim(Y)[2])*abs)/(rho*error_b_0[i]) )
    iteration<<-i
    i <- i+1
  }
  return(list(b_0[,i-1],B[,,i-1],Q[,,i-1],S[,,i-1]))
}

# Initial Value Generator for Cross Validation
CV_5fold_MTR_initial <- function(Xgroup,Ygroup,rho,maxiter,tole){
  initial_beta_0 <- array(0,c(dim(Ygroup)[2]-1,5))
  initial_beta <- array(0,c(dim(Xgroup)[2]-1,dim(Ygroup)[2]-1,5))
  for(i in 1:5){
    CV_i<<-i
    initial_cv <- ADMM_3block(0,0,as.matrix(filter(Xgroup,group!=i)[,-1]),as.matrix(filter(Ygroup,group!=i)[,-1]),0,0.00001,rho,maxiter,tole)
    initialbeta <- as.matrix(as.data.frame(initial_cv[4]))
    initialbeta_0 <- as.matrix(as.data.frame(initial_cv[1]))
    initial_beta[,,i] <- initialbeta
    initial_beta_0[,i] <- initialbeta_0
  }
  return(list(initial_beta_0,initial_beta))
}

# CV Fit Function for LGSMTR
CVfit_LGSMTR <- function(filtered_indep,filtered_resp,group,maxgrid,gridlength=50,ite_ADMM=100){
  k <- dim(filtered_resp)[2]
  p <- dim(filtered_indep)[2]
  X <- as.matrix(filtered_indep)
  Y <- as.matrix(filtered_resp)
  
  grid <- seq(0.00001,maxgrid,length=gridlength)
  record_1 <- matrix(0,gridlength,2)
  record_1[,1] <- grid
  Xgroup <- data.frame(group,X)
  Ygroup <- data.frame(group,Y)
  X_GRP <- as.matrix(Xgroup)
  Y_GRP <- as.matrix(Ygroup)
  
  init_cv <- CV_5fold_MTR_initial(Xgroup,Ygroup,0.1,ite_ADMM,3)
  initial_beta_0 <- array(as.numeric(unlist(init_cv[1])), dim=c(k,5))
  initial_beta <- array(as.numeric(unlist(init_cv[2])), dim=c(p,k,5))  
    
  indexcorr <- combn(1:(dim(Ygroup)[2]-1),2)
  
  ZZZ <- array(0,c(dim(Ygroup)[2]-1,dim(indexcorr)[2],5))
  for (m in 1:5){
    for(u in 1:dim(indexcorr)[2]){
      ZZZ[indexcorr[1,u],u,m] <- cor(Y_GRP[Y_GRP[,1]!=m,][,-1][,indexcorr[1,u]],Y_GRP[Y_GRP[,1]!=m,][,-1][,indexcorr[2,u]])
      ZZZ[indexcorr[2,u],u,m] <- -sign(ZZZ[indexcorr[1,u],u,m])*ZZZ[indexcorr[1,u],u,m]
      u <- u+1
    }
    m <- m+1
  }  
    
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  clusterExport(cl, c("ADMM_3block", "ADMM_3block_Z", "sylvester" , "Schur", "progress_bar", "hadamard.prod"))
  registerDoParallel(cl)  
  
  chunk <- foreach(i = 1:length(grid), .combine= 'cbind') %dopar% {
    CV_error_GSLMTR_2 <- rep(0,5)
    for(m in 1:5){
      temp <- ADMM_3block_Z(initial_beta_0[,m],initial_beta[,,m],X_GRP[X_GRP[,1]!=m,][,-1],Y_GRP[Y_GRP[,1]!=m,][,-1],ZZZ[,,m],record_1[i,1],record_1[i,1]/sqrt(k),0.1,ite_ADMM,3)
      tempbeta <- as.matrix(as.data.frame(temp[3]))
      tempbeta_0 <- as.matrix(as.data.frame(temp[1]))
      CV_error_GSLMTR_2[m] <- mean((Y_GRP[Y_GRP[,1]==m,][,-1]-rep(1,dim(Y_GRP[Y_GRP[,1]==m,][,-1])[1])%*%t(tempbeta_0)-X_GRP[X_GRP[,1]==m,][,-1]%*%tempbeta)^2)
      m <- m+1 
    }
    c(mean(CV_error_GSLMTR_2))
  }
  stopCluster(cl)
  record_1[,2] <- as.vector(chunk)
  rec <- as.data.frame(record_1)
  colnames(rec) <- c("lambda_1","CV_Error")
  
  opt_lambda_1 <- rec[order(rec$CV_Error),][1,1]
  
  SS <- ADMM_3block(0,0,X,Y,0,0.00001,0.1,ite_ADMM*10,3)
  initial_0 <- as.matrix(as.data.frame(SS[1]))
  initialSol <- as.matrix(as.data.frame(SS[4]))
  
  ZZ <- matrix(0,k,dim(indexcorr)[2])
  for(i in 1:dim(indexcorr)[2]){
    ZZ[indexcorr[1,i],i] <- cor(filtered_resp[,indexcorr[1,i]],filtered_resp[,indexcorr[2,i]])
    ZZ[indexcorr[2,i],i] <- -sign(cor(filtered_resp[,indexcorr[1,i]],filtered_resp[,indexcorr[2,i]]))*cor(filtered_resp[,indexcorr[1,i]],filtered_resp[,indexcorr[2,i]])
    i <- i+1
  }
  
  Sol_GSLMTR_2=ADMM_3block_Z(initial_0,initialSol,X,Y,ZZ,opt_lambda_1,opt_lambda_1/sqrt(k),0.1,ite_ADMM*10,3)
  SS_GSLMTR_2 <- as.matrix(as.data.frame(Sol_GSLMTR_2[3]))
  drugname <- c("AEW541","Nilotinib","17-AAG","PHA-665752","Lapatinib","Nutlin-3","AZD0530","PF2341066","L-685458","ZD-6474","Panobinostat","Sorafenib","Irinotecan","Topotecan","LBW242","PD-0325901","PD-0332991","Paclitaxel","AZD6244","PLX4720","RAF265","TAE684","TKI258","Erlotinib")
  colnames(SS_GSLMTR_2) <- drugname[-11]
  row.names(SS_GSLMTR_2) <- as.numeric(gsub("V",0,colnames(filtered_indep)))
  
  output <- list(rec,SS_GSLMTR_2)
  names(output) <- c("CV Error", "Estimated Coefficient")
  return(output)
}
