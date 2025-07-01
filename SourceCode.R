## Functions for ZIPM data generation, density calculation, model fitting, 
## and confidence interval calculations.
library(matrixStats)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gridExtra)

# generate data from a ZIPM model
rZIPM <- function(I,J,pi,epsilon,mu,nu,t,seed=NULL){
  if(pi<=0 | pi > 0.5 | mu <=0 | nu <=0 | epsilon <=0 | epsilon >= 1 | length(t) != I){
    stop("Problem with parameters. Check conditions.")
  }
  if(!is.null(seed)){set.seed(seed)}
  
  Y <- rbinom(J,1,pi)
  Z <- matrix(rbinom(I*J,1,epsilon),nrow=I,ncol=J)
  M <- matrix(NA,nrow=I,ncol=J)
  for(j in 1:J){
    if(Y[j]==0){M[,j] <- rpois(I,t*nu)
    }else{M[,j] <- rpois(I,t*mu)}
  }
  N <- Z*M
  
  return(N)
}

# calculate density of data from a ZIPM model
dZIPM <- function(N,I,J,pi,epsilon,mu,nu,t,log=FALSE){
  
  dN <- matrix(NA,nrow=I,ncol=J)
  for(i in 1:I){for(j in 1:J){
    dN[i,j] <- pi*epsilon*dpois(N[i,j],t[i]*mu)+(1-pi)*epsilon*dpois(N[i,j],t[i]*nu)+(1-epsilon)*(0^(N[i,j]))
  }}
  logD <- sum(log(dN))
  if(log){return(logD)}else{return(exp(logD))}
  
}

# fit frequentist model via EM (three helper functions, main "EM" function,
# and a function to try multiple random EM initializers)
get_yhat <- function(N,I,J,pi,epsilon,mu,nu,t){
  
  one_eq_ij <- 0^N
  one_uneq_ij <- 1-(0^N)
  t_uneq_j <- apply(t*one_uneq_ij,2,sum)
  n_j <- apply(N,2,sum)
  
  yhat <- numeric(J)
  for(j in 1:J){
    term1<-exp(-t_uneq_j[j]*(nu-mu)+n_j[j]*log(nu/mu)+sum((one_eq_ij[,j])*log((epsilon*exp(-t*nu)+(1-epsilon))/(epsilon*exp(-t*mu)+(1-epsilon)))))
    yhat[j] <- pi/(pi+(1-pi)*term1)
  }
  
  
  # yhat <- numeric(J)
  # for(j in 1:J){
  #   dens_A <- dpois(N[,j],t*mu)
  #   dens_B <- dpois(N[,j],t*nu)
  #   dens_0 <- (0^N[,j])
  #   yhat[j] <- pi/(pi+(1-pi)*prod((epsilon*dens_B+(1-epsilon)*dens_0)/(epsilon*dens_A+(1-epsilon)*dens_0)))
  # }

  
  return(yhat)
}
get_yzhat <- function(N,I,J,pi,epsilon,mu,nu,t,yhat){
  
  one_eq_ij <- 0^N
  yzhat <- matrix(0,nrow=I,ncol=J)
  for(i in 1:I){for(j in 1:J){
    term1 <- (1+((1-epsilon)/epsilon)*exp(t[i]*mu))^(-one_eq_ij[i,j])
    yzhat[i,j] <- yhat[j]*term1
    
  }}
  
  # yzhat <- matrix(0,nrow=I,ncol=J)
  # for(i in 1:I){for(j in 1:J){
  #   dens_A <- dpois(N[i,j],t[i]*mu)
  #   dens_0 <- 0^N[i,j]
  #   
  #   term1 <- epsilon*dens_A
  #   term2 <- (1-epsilon)*dens_0
  #   yzhat[i,j] <- yhat[j]*term1/(term1+term2)
  # }}
  
  
  return(yzhat)
}
get_onelessyzhat <- function(N,I,J,pi,epsilon,mu,nu,t,yhat){
  
  one_eq_ij <- 0^N
  onelessyzhat <- matrix(0,nrow=I,ncol=J)
  for(i in 1:I){for(j in 1:J){
    term1 <- (1+((1-epsilon)/epsilon)*exp(t[i]*nu))^(-one_eq_ij[i,j])
    onelessyzhat[i,j] <- (1-yhat[j])*term1
  }}
  
  # yzhat <- matrix(0,nrow=I,ncol=J)
  # for(i in 1:I){for(j in 1:J){
  #   dens_A <- dpois(N[i,j],t[i]*mu)
  #   dens_0 <- 0^N[i,j]
  #   
  #   term1 <- epsilon*dens_A
  #   term2 <- (1-epsilon)*dens_0
  #   yzhat[i,j] <- yhat[j]*term1/(term1+term2)
  # }}
  
  
  return(onelessyzhat)
}
get_zhat <- function(yzhat,onelessyzhat){
  zhat <- yzhat + onelessyzhat
  
  # zhat <- matrix(0,nrow=I,ncol=J)
  # for(i in 1:I){for(j in 1:J){
  #   
  #   dens_A <- dpois(N[i,j],t[i]*mu)
  #   dens_B <- dpois(N[i,j],t[i]*nu)
  #   dens_0 <- 0^N[i,j]
  #   
  #   term1 <- epsilon*(pi*dens_A+(1-pi)*dens_B)
  #   term2 <- (1-epsilon)*dens_0
  #   
  #   zhat[i,j] <- term1/(term1+term2)
  # }}
  
  return(zhat)
}
EM <- function(N,t,seed=NULL,omega0=NULL,tol=0.0001,update_t=FALSE){
  
  if(!is.null(seed)){set.seed(seed)}
  if(is.null(omega0)){
    pi <- runif(1,0,.5)
    epsilon <- runif(1)
    mu <- max(.1,mean(N/t)+runif(1,-2,2))
    nu <- max(.1,mean(N/t)+runif(1,-2,2))
  }else{
    pi <- omega0[1]
    epsilon <- omega0[2]
    mu <- omega0[3]
    nu <-omega0[4]
  }
  I <- nrow(N)
  J <- ncol(N)
  
  diff <- Inf
  iter <- 1
  while(diff > tol & iter < 5000){
    
    yhat <- get_yhat(N,I,J,pi,epsilon,mu,nu,t)
    yzhat <- get_yzhat(N,I,J,pi,epsilon,mu,nu,t,yhat)
    onelessyzhat <- get_onelessyzhat(N,I,J,pi,epsilon,mu,nu,t,yhat)
    zhat <- get_zhat(yzhat,onelessyzhat)
    
    (pi_new <- mean(yhat))
    (epsilon_new <- mean(zhat))
    (mu_new <- sum(apply(N,2,sum)*yhat)/sum(t*yzhat))
    (nu_new <- sum(apply(N,2,sum)*(1-yhat))/sum(t*(onelessyzhat)))
    
    if(update_t){
      for(i in 1:I){
        t[i] <- optim(par=t[i],fn=function(t){
          value <- -sum((yhat*zhat[i,])*(-t*mu + N[i,]*log(t*mu))+
                          ((1-yhat)*zhat[i,])*(-t*nu+N[i,]*log(t*nu))+
                          log(0^(N[i,]*(1-zhat[i,]))))
          if(is.infinite(value)){return(-100000)}
          return(value)
        },method="L-BFGS-B",lower=.01)$par
      }
    }
    
    if(pi_new == 1 | pi_new == 0 | is.nan(pi_new)){pi_new <- 0.5}
    if(epsilon_new == 1 | epsilon_new == 0 | is.nan(epsilon_new)){epsilon_new <- 0.5}
    if(is.nan(mu_new)){mu_new <- mean(N/t)+runif(1,0,1)}
    if(is.nan(nu_new)){nu_new <- mean(N/t)+runif(1,0,1)}
    
    diff <- sum(abs(c(pi_new,epsilon_new,mu_new,nu_new)-c(pi,epsilon,mu,nu)))
    pi <- pi_new
    epsilon <- epsilon_new
    mu <- mu_new
    nu <- nu_new
    iter <- iter + 1
  }
  
  if(pi <= 0.5){
    return(list(pi=pi,epsilon=epsilon,mu=mu,nu=nu,theta=mu/nu,zhat=zhat,yhat=yhat,t=t))
  }else{
    return(list(pi=1-pi,epsilon=epsilon,mu=nu,nu=mu,theta=nu/mu,zhat=zhat,yhat=yhat,t=t))
  }
}
EM_multistart <- function(N,t,tol=0.0001,starts=2,update_t=FALSE){
  I <- nrow(N)
  J <- ncol(N)
  res <- replicate(starts,{
    EM(N,t,tol=tol,update_t=update_t)})
  d <- apply(res,2,function(res){dZIPM(N,I,J,res$pi,res$epsilon,res$mu,res$nu,res$t,log=TRUE)})
  return(res[,which.max(d)])
}

# calculate standard error and 100x(1-alpha) confidence intervals for theta
get_se <- function(N,I,J,pi,epsilon,mu,nu,t){
  
  # calculate notational constants
  one_uneq_ij <- 1-(0^N)
  one_uneq_j <- apply(one_uneq_ij,2,sum)
  one_uneq <- sum(one_uneq_j)
  t_uneq_j <- apply(t*one_uneq_ij,2,sum)
  one_eq_ij <- 0^N
  one_eq_j <- apply(one_eq_ij,2,sum)
  one_eq <- sum(one_eq_j)
  t_eq_j <- apply(t*one_eq_ij,2,sum)
  n_j <- apply(N,2,sum)
  K <- I*J
  
  # get q_j (good)
  logA <- -t_uneq_j*mu + n_j*log(mu) + apply(one_eq_ij,2,function(one_eq_ij_j){
    sum(one_eq_ij_j*log(epsilon*exp(-t*mu)+(1-epsilon)))
  })
  logB <- -t_uneq_j*nu + n_j*log(nu) + apply(one_eq_ij,2,function(one_eq_ij_j){
    sum(one_eq_ij_j*log(epsilon*exp(-t*nu)+(1-epsilon)))
  })
  q_j <- pi/(pi+(1-pi)*(exp(logB - logA)))
  rm(logA,logB)
  
  # get rho (good)
  sum <- 0
  for(i in 1:I){for(j in 1:J){
    if(one_eq_ij[i,j] == 1){
      sum <- sum + exp(log(q_j[j])-log(epsilon+(1-epsilon)*exp(t[i]*mu)))+
        exp(log(1-q_j[j])-log(epsilon+(1-epsilon)*exp(t[i]*nu)))
    }
  }}
  # sum <- 0
  # for(i in 1:I){for(j in 1:J){
  #   sum <- sum + one_eq_ij[i,j]*((q_j[j])/(epsilon+(1-epsilon)*exp(t[i]*mu))+
  #                                  (1-q_j[j])/(epsilon+(1-epsilon)*exp(t[i]*nu)))
  # }}
  rho <- one_uneq/K + (epsilon/K)*sum
  rm(sum)
  
  # get (log_)psi (good!)
  #term1 <- pi*exp(-t_uneq_j*mu)*mu^n_j*unlist(lapply(1:J,function(j){prod((epsilon*(exp(-t*mu)-1)+1)^one_eq_ij[,j])}))
  log_term1 <- log(pi)-(t_uneq_j*mu)+n_j*log(mu)+unlist(lapply(1:J,function(j){sum(one_eq_ij[,j]*log(epsilon*(exp(-t*mu)-1)+1))}))
  #term2 <- (1-pi)*exp(-t_uneq_j*nu)*nu^n_j*unlist(lapply(1:J,function(j){prod((epsilon*(exp(-t*nu)-1)+1)^one_eq_ij[,j])}))
  log_term2 <- log(1-pi)-(t_uneq_j*nu)+n_j*log(nu)+unlist(lapply(1:J,function(j){sum(one_eq_ij[,j]*log(epsilon*(exp(-t*nu)-1)+1))}))
  log_psi_j <- unlist(lapply(1:J,function(j){logSumExp(c(log_term1[j],log_term2[j]))}))
  rm(log_term1,log_term2)
  
  # get phi (good)
  phi_j <- array(NA,c(4,1,J))
  phi_j[1,1,] <- 1/(pi*(1-pi))
  phi_j[2,1,] <- unlist(lapply(1:J,function(j){
    term1 <- one_eq_ij[,j]*(exp(-t*mu)-exp(-t*nu))
    term2 <- (epsilon*(exp(-t*mu)-1)+1)*(epsilon*(exp(-t*nu)-1)+1)
    sum(term1/term2)
  }))
  phi_j[3,1,] <- n_j/mu-t_uneq_j-epsilon*unlist(lapply(1:J,function(j){
    sum(one_eq_ij[,j]*exp(log(t)-t*mu-log(epsilon*(exp(-t*mu)-1)+1)))
  }))
  phi_j[4,1,] <- n_j/nu-t_uneq_j-epsilon*unlist(lapply(1:J,function(j){
    sum(one_eq_ij[,j]*exp(log(t)-t*nu-log(epsilon*(exp(-t*nu)-1)+1)))
  }))
  
  # get chi (good)
  chi_i_1 <- chi_i_0 <- array(NA,c(4,1,I))
  chi_i_1[1,1,] <- chi_i_0[1,1,] <- 0
  chi_i_1[2,1,] <- chi_i_0[2,1,] <- 1/(epsilon*(1-epsilon))
  chi_i_1[3,1,] <- chi_i_0[4,1,] <- t
  chi_i_1[4,1,] <- chi_i_0[3,1,] <- 0
  
  ## term 1 matrix (good)
  term1_mat <- matrix(0,nrow=4,ncol=4)
  term1_mat[1,1] <- (J*((1-2*pi)*mean(q_j)+pi^2))/((pi^2)*(1-pi)^2)
  term1_mat[2,2] <- (K*((1-2*epsilon)*rho+epsilon^2))/(epsilon^2*(1-epsilon)^2)
  term1_mat[3,3] <- K*(mean(t(N)*q_j)/(mu^2))
  term1_mat[4,4] <- K*(mean(t(N)*(1-q_j))/(nu^2))
  
  ## term 2 matrix
  log_num <- -t_uneq_j*(mu+nu)+n_j*log(mu*nu)+unlist(lapply(1:J,function(j){sum(one_eq_ij[,j]*(log(epsilon*(exp(-t*mu)-1)+1)+log(epsilon*(exp(-t*nu)-1)+1)))}))
  log_denom <- 2*log_psi_j
  constant <- exp(log_num-log_denom)
  mat_j <- lapply(1:J,function(j){(constant[j])*(t(t(phi_j[,,j])) %*% t(phi_j[,,j]))})
  mat_j <- array(unlist(mat_j),c(4,4,J))
  term2_mat_a <- -pi*(1-pi)*apply(mat_j,c(1,2),sum)
  
  term2_mat_b <- matrix(0,nrow=4,ncol=4)
  for(i in 1:I){for(j in 1:J){
    term1 <- exp(-t[i]*mu-2*log(epsilon*(exp(-t[i]*mu)-1)+1)+log(q_j[j]))*(t(t(chi_i_1[,,i]))%*%t(chi_i_1[,,i]))
    #term1 <- (exp(-t[i]*mu)/((epsilon*(exp(-t[i]*mu)-1)+1)^2))*(q_j[j])*(chi_i_1[,,i,drop=FALSE] %*% t(chi_i_1[,,i]))
    term2 <- exp(-t[i]*nu-2*log(epsilon*(exp(-t[i]*nu)-1)+1)+log(1-q_j[j]))*(t(t(chi_i_0[,,i]))%*%t(chi_i_0[,,i]))
    #term2 <- (exp(-t[i]*nu)/((epsilon*(exp(-t[i]*nu)-1)+1)^2))*(chi_i_0[,,i,drop=FALSE] %*% t(chi_i_0[,,i]))*(1-q_j[j])
    term2_mat_b <- term2_mat_b + one_eq_ij[i,j]*(term1+term2)
  }}
  term2_mat_b <- -epsilon*(1-epsilon)*term2_mat_b
  term2_mat <- term2_mat_a+term2_mat_b
  
  ### calculate observed information matrix
  I_N <- term1_mat+term2_mat
  
  ### calculate standard error, tau
  I_11 <- I_N[1:2,1:2]
  I_21 <- I_N[3:4,1:2]
  I_12 <- I_N[1:2,3:4]
  I_22 <- I_N[3:4,3:4]
  mat1 <- matrix(c(1/nu,-mu/(nu^2)),nrow=1)
  mat2 <- solve(I_22 - (I_21 %*% solve(I_11) %*% I_12))
  tau_squared <- c(K*(mat1 %*% mat2 %*% t(mat1)))
  tau <- sqrt(tau_squared)
  
  return(tau)
}
get_conf <- function(theta,tau,I,J,alpha=0.05){
  K <- I*J
  return(theta+c(tau/sqrt(K)*qnorm(alpha/2),-tau/sqrt(K)*qnorm(alpha/2)))
}
