######################
# Manifold EL code iterate
######################
library(rgl) #3D plot
library(rotasym) #Generate von Mises-Fisher dist.
library(RiemBase) #Compute Frechet mean for spherical data
library(progress)
######################
#C-matrix
C_mat <- function(th,ph){
  CM <- array(0,dim=c(3,3))
  CM[1,1] <- cos(th)*cos(ph);CM[1,2] <- -sin(ph);CM[1,3] <- sin(th)*cos(ph)
  CM[2,1] <- cos(th)*sin(ph);CM[2,2] <-  cos(ph);CM[2,3] <- sin(th)*sin(ph)
  CM[3,1] <- -sin(th);       CM[3,2] <- 0;       CM[3,3] <- cos(th)
  return(CM)
}
#g function
g_func_med <- function(x,th,ph){#x: vector, th: theta, ph: phi
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  m <- c(m1,m2,m3)
  A <- x - sum(m*x)*m
  B <- sqrt(sum(A^2))
  C <- A/B
  D <- t(C_mat(th,ph))%*%C
  E <- c(D[1],D[2])
  return(E)
}

######################
#g function theta = pi, phi = 0
g_func3 <- function(x){#x: vector, th: theta, ph: phi
  g1 <- x[1]*cos(th)*cos(ph) + x[2]*cos(th)*sin(ph) - x[3]*sin(th)
  g2 <- -x[1]*sin(ph) + x[2]*cos(ph)
  A <- acos(-x[3])
  B <- sqrt(1-(-x[3])^2)
  C <- -2*sin(1)*(A/B)
  D <- c(g1,g2)*C
  return(D)
}
#######################
#Compute dual form iteratively
MC <- 1000
n <- 50 #sample size
boot <-199 #number of bootstrap replications
######################
#Frechet mean
#first sample
dim <- 2
theta1 <- 0
phi1 <- 0
mu1 <- sin(theta1)*cos(phi1)
mu2 <- sin(theta1)*sin(phi1)
mu3 <- cos(theta1)
mu_fir <- c(mu1,mu2,mu3)
#second sample
theta2 <- 0
phi2 <- 0
mu12 <- sin(theta2)*cos(phi2)
mu22 <- sin(theta2)*sin(phi2)
mu32 <- cos(theta2)
mu_sec <- c(mu12,mu22,mu32)
#################################
lam_save_1 <- array(0,dim=c(MC,2))
lam_save_2 <- array(0,dim=c(MC,2))
lam_save_1b <- array(0,dim=c(MC,boot,2))
lam_save_2b <- array(0,dim=c(MC,boot,2))
EL_dual_opt_1 <- rep(0,MC)
EL_dual_opt_2 <- rep(0,MC)
EL_dual_opt_3 <- rep(0,MC)
EL_dual_opt_1b <- array(0,dim=c(MC,boot))
EL_dual_opt_2b <- array(0,dim=c(MC,boot))
EL_dual_opt_3b <- array(0,dim=c(MC,boot))
quant_TS <- array(0,dim=c(MC,3))
####################################
pb <- progress_bar$new(total = MC,
                       format = "[:bar] :percent : :eta",
                       clear = TRUE)
kappa1 <- 15 #parameter of concentration of von Mises-Fisher dist. 
kappa2 <- 25 #parameter of concentration of von Mises-Fisher dist. 
for(j in 1:MC){
  pb$tick()
  ##########################################
  #Generate data from von Mises-Fisher dist.
  X1 <- r_vMF(n,mu_fir,kappa1)
  X2 <- r_vMF(n,mu_sec,kappa2)
  X_merge <- rbind(X1,X2)
  n_merge <- 2*n
  ##########################################
  #Compute Frechet median/mean of marged sample
  th_max <- 64
  ph_max <- 128
  Fre_med <- rep(0,th_max*ph_max)
  Fre_mean <- rep(0,th_max*ph_max)
  th_cand <- (seq(1:th_max)-1)/(th_max)*(pi/2)
  ph_cand <- (seq(1:ph_max)-1)/ph_max*(2*pi)
  for(t in 1:ph_max){
    for(k in 1:th_max){
      Fmu1 <- sin(th_cand[k])*cos(ph_cand[t])
      Fmu2 <- sin(th_cand[k])*sin(ph_cand[t])
      Fmu3 <- cos(th_cand[k])
      Fmu <- c(Fmu1,Fmu2,Fmu3)
      A_m <- 0
      B_m <- 0
      for(r in 1:n_merge){
        A_m <- A_m + acos(sum(Fmu*X_merge[r,]))
        B_m <- B_m + acos(sum(Fmu*X_merge[r,]))^2
      }
      Fre_med[(k-1)*ph_max+t] <- A_m/n_merge
      Fre_mean[(k-1)*ph_max+t] <- B_m/n_merge
    }
  }
  ind_Fmed <- which(Fre_med == min(Fre_med))
  ind_Fmed2 <- floor(ind_Fmed/ph_max)
  ind_Fmed3 <- ind_Fmed-ph_max*ind_Fmed2
  ind_Fmed_th <- ind_Fmed2+1
  ind_Fmed_ph <- ind_Fmed3+1
  theta_m <- max(th_cand[ind_Fmed_th])  #theta of the Frechet median computed from the marged sample
  phi_m <- max(ph_cand[ind_Fmed_ph])   #phi of the Frechet median computed from the marged sample
  #######################
  #Compute dual form of EL
  #first sample
  EL_dual_1 <- function(lam){
    g_vec2 <- rep(0,n)
    #when null theta neq pi
    for(w1 in 1:n){
      g_vec2[w1] <- -2*log(1+sum(lam*g_func_med(X1[w1,],theta_m,phi_m)))
    }
    #when null theta = pi
    #for(i in 1:n){
    #  g_vec2[i] <- -2*log(1+sum(lam*g_func3(X[i,])))
    #}
    return(sum(g_vec2))
  }
  #second sample
  EL_dual_2 <- function(lam){
    g_vec22 <- rep(0,n)
    #when null theta neq pi
    for(w2 in 1:n){
      g_vec22[w2] <- -2*log(1+sum(lam*g_func_med(X2[w2,],theta_m,phi_m)))
    }
    #when null theta = pi
    #for(i in 1:n){
    #  g_vec2[i] <- -2*log(1+sum(lam*g_func3(X[i,])))
    #}
    return(sum(g_vec22))
  }
  lambda01 <- optim(c(0,0),EL_dual_1)$par
  lambda02 <- optim(c(0,0),EL_dual_2)$par
  EL_dual_opt_1[j] <- -EL_dual_1(lambda01) 
  EL_dual_opt_2[j] <- -EL_dual_2(lambda02) 
  EL_dual_opt_3[j] <- EL_dual_opt_1[j] + EL_dual_opt_2[j]
  lam_save_1[j,] <- lambda01
  lam_save_2[j,] <- lambda02
  ######################################
  #sample index for bootstrap
  ind_boot_1 <- sample(c(1:n_merge),size=n*boot,replace=TRUE)
  ind_boot_2 <- sample(c(1:n_merge),size=n*boot,replace=TRUE)
  m_boot_1 <- matrix(ind_boot_1,ncol=boot)
  m_boot_2 <- matrix(ind_boot_2,ncol=boot)
  for(bt in 1:boot){
    X1b <- X_merge[m_boot_1[,bt],]
    X2b <- X_merge[m_boot_2[,bt],]
    #######################
    #Compute dual form of EL
    #first sample
    EL_dual_1b <- function(lam){
      g_vec2_1b <- rep(0,n)
      #when null theta neq pi
      for(i1 in 1:n){
        g_vec2_1b[i1] <- -2*log(1+sum(lam*g_func_med(X1b[i1,],theta_m,phi_m)))
        }
      #when null theta = pi
      #for(i1 in 1:n){
      #  g_vec2_1b[i1] <- -2*log(1+sum(lam*g_func3(Xb[i1,])))
      #}
      return(sum(g_vec2_1b))
      }
    #second sample
    EL_dual_2b <- function(lam){
      g_vec22_2b <- rep(0,n)
      #when null theta neq pi
      for(i2 in 1:n){
        g_vec22_2b[i2] <- -2*log(1+sum(lam*g_func_med(X2b[i2,],theta_m,phi_m)))
        }
      #when null theta = pi
      #for(i2 in 1:n){
      #  g_vec22_2b[i2] <- -2*log(1+sum(lam*g_func3(X2b[i2,])))
      #}
      return(sum(g_vec22_2b))
      }
    lambda_1b <- optim(c(0,0),EL_dual_1b)$par
    lambda_2b <- optim(c(0,0),EL_dual_2b)$par
    EL_dual_opt_1b[j,bt] <- -EL_dual_1b(lambda_1b) 
    EL_dual_opt_2b[j,bt] <- -EL_dual_2b(lambda_2b) 
    EL_dual_opt_3b[j,bt] <- EL_dual_opt_1b[j,bt] + EL_dual_opt_2b[j,bt]
    lam_save_1b[j,bt,] <- lambda_1b
    lam_save_2b[j,bt,] <- lambda_2b
    }
  ind_quant <- which(EL_dual_opt_3b[j,]<=5000)
  quant1 <- as.numeric(quantile(EL_dual_opt_3b[j,ind_quant],c(0.90,0.95,0.99)))
  quant_TS[j,] <- quant1
  Sys.sleep(1 / MC)
  }

#####################################
#compute size
ind_quant_NA <- which(is.na(quant_TS[,1])==T)
ind_quant <- which(is.na(quant_TS[,1])==F)
n_NA <- length(ind_quant_NA)
n_TRUE <- length(ind_quant)
EL_stat <- EL_dual_opt_3[ind_quant]
quant_TS1 <- quant_TS[ind_quant,1]
quant_TS2 <- quant_TS[ind_quant,2]
quant_TS3 <- quant_TS[ind_quant,3]
size_90 <- (length(which(EL_stat>quant_TS1))+n_NA)/MC
size_95 <- (length(which(EL_stat>quant_TS2))+n_NA)/MC
size_99 <- (length(which(EL_stat>quant_TS3))+n_NA)/MC
size_90
size_95
size_99


