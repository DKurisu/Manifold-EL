######################
# Manifold EL code iterate
######################
library(rgl) #3D plot
library(rotasym) #Generate von Mises-Fisher dist.
library(RiemBase) #Compute Frechet mean for spherical data
library(progress)
######################
#Lambda/C-function
#Lambda matrix function
Lam_11_f <- function(x,th,ph){#x: vector, th, ph: angle of the null Frechet mean
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  m <- c(m1,m2,m3)
  AL1 <- sum(m*x)
  BL1 <- sqrt(1-AL1^2)
  gL1 <- x[1]*cos(th)*cos(ph) + x[2]*cos(th)*sin(ph) - x[3]*sin(th)
  CL1 <- 2*sin(1)^2*(1/(BL1^2) - acos(AL1)*AL1/(BL1^3))*(gL1^2)
  CL2 <- -2*(acos(AL1)/BL1)*(cos(1)-2*sin(1))*AL1
  DL1 <- CL1 + CL2
  return(DL1)
}
Lam_12_f <- function(x,th,ph){#x: vector, th, ph: angle of the null Frechet mean
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  m <- c(m1,m2,m3)
  AL2 <- sum(m*x)
  BL2 <- sqrt(1-AL2^2)
  gL2 <- x[1]*cos(th)*cos(ph) + x[2]*cos(th)*sin(ph) - x[3]*sin(th)
  gL3 <- -x[1]*sin(ph) + x[2]*cos(ph)
  DL2 <- 2*sin(1)^2*(1/(BL2^2) - acos(AL2)*AL2/(BL2^3))*gL2*gL3
  return(DL2)
}
Lam_22_f <- function(x,th,ph){#x: vector, th, ph: angle of the null Frechet mean
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  m <- c(m1,m2,m3)
  AL3 <- sum(m*x)
  BL3 <- sqrt(1-AL3^2)
  gL4 <- -x[1]*sin(ph) + x[2]*cos(ph)
  CL3 <- 2*sin(1)^2*(1/(BL3^2) - acos(AL3)*AL3/(BL3^3))*(gL4^2)
  CL4 <- -2*(acos(AL3)/BL3)*(cos(1)-2*sin(1))*AL3
  DL3 <- CL3 + CL4
  return(DL3)
}
Lam_f <- function(x,th,ph){
  L <- array(0,dim=c(2,2))
  L11 <- Lam_11_f(x,th,ph)
  L12 <- Lam_12_f(x,th,ph)
  L22 <- Lam_22_f(x,th,ph)
  L[1,1] <- L11
  L[1,2] <- L12
  L[2,1] <- L12
  L[2,2] <- L22
  return(L)
}
######################
#C matrix function
#g function
g_func <- function(x,th,ph){#x: vector, th: theta, ph: phi
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  m <- c(m1,m2,m3)
  gf1 <- x[1]*cos(th)*cos(ph) + x[2]*cos(th)*sin(ph) - x[3]*sin(th)
  gf2 <- -x[1]*sin(ph) + x[2]*cos(ph)
  Ag <- acos(sum(m*x))
  Bg <- sqrt(1-(sum(m*x))^2)
  Cg <- -2*sin(1)*(Ag/Bg)
  Dg <- c(gf1,gf2)*Cg
  return(Dg)
}
C_f <- function(x,gm,th,ph){
  g0 <- g_func(x,th,ph)
  c0<- (g0-gm)%*%t(g0-gm)
  return(c0)
}
C_f1 <- function(x,th,ph){
  g0 <- g_func(x,th,ph)
  c0<- g0%*%t(g0)
  return(c0)
}
########################
#logarithmic map
#lotate north pole matrix
lotate_north <- function(th,ph){
  MA <- array(0,dim=c(3,3))
  MB <- array(0,dim=c(3,3))
  MC <- array(0,dim=c(3,3))
  MA[1,1] <- cos(th);MA[1,3] <- sin(th)
  MA[2,2] <- 1
  MA[3,1] <- -sin(th);MA[3,3] <- cos(th) 
  MB[1,1] <- cos(ph);MB[1,2] <- -sin(ph)
  MB[2,1] <- sin(ph);MB[2,2] <- cos(ph)
  MB[3,3] <- 1
  MC <- t(MB%*%MA)
  return(MC)
}
log_map <- function(Fm,th,ph){
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  m <- c(m1,m2,m3)
  Alm <- sum(Fm*m)
  Blm <- Fm - Alm*m
  Clm <- acos(Alm)*Blm/sqrt(sum(Blm^2))
  Dlm <- lotate_north(th,ph)
  Elm <- Dlm%*%Clm
  Flm <- c(Elm[1,],Elm[2,])
  return(Flm)
}
#######################
#true Frechet mean
theta <- 0
phi <- 0
mu1 <- sin(theta)*cos(phi)
mu2 <- sin(theta)*sin(phi)
mu3 <- cos(theta)
mu <- c(mu1,mu2,mu3)
#######################
#Compute Wald statistics iteratively
MC <- 1000
kappa <- 1/2 #parameter of concentration of von Mises-Fisher dist.
n <- 200 #sample size
dim <- 2
Wald_save <- numeric(MC)
pb <- progress_bar$new(total = MC,
                       #format = "[:bar] :percent: :eta",
                       format = "[:bar] :percent: :elapsed",
                       clear = TRUE)
########################
#null Frechet mean
th_n <- 0 #null theta
ph_n <- 0 #null theta
########################
for(j in 1:MC){
  pb$tick()
  ##########################################
  #Generate data from von Mises-Fisher dist.
  X <- r_vMF(n,mu,kappa)
  #######################
  #Compute Frechet mean
  data  = list()
  for (i in 1:n){
    data[[i]] <- c(X[i,1],X[i,2],X[i,3])
  }
  data <- riemfactory(data, name="sphere")
  out1 <- rbase.mean(data)
  out1_r <- sqrt(out1$x[1,1]^2 + out1$x[2,1]^2 + out1$x[3,1]^2)
  FM <- out1$x/out1_r
  #th_FM <- acos(FM[3])
  #ph_FM <- sign(FM[2])*acos(FM[1]/sqrt(FM[1]^2+FM[2]^2))
  th_FM <- th_n
  ph_FM <- ph_n
  #######################
  #Compute Wald statistics
  Lam_hat <- array(0,dim=c(2,2))
  C_hat <- array(0,dim=c(2,2))
  g_hat_mean <- array(0,dim=c(2,1))
  for(i in 1:n){
    g_hat_mean <- g_hat_mean + g_func(X[i,],th_FM,ph_FM)
  }
  g_hat_mean <- g_hat_mean/n
  for(i in 1:n){
    Lam_hat <- Lam_hat + Lam_f(X[i,],th_FM,ph_FM) 
    C_hat <- C_hat + C_f(X[i,],g_hat_mean, th_FM,ph_FM)
    #C_hat <- C_hat + C_f1(X[i,],th_FM,ph_FM)
  }
  Lam_hat <- Lam_hat/n
  C_hat_inv <- solve(C_hat/n)
  V_Wald <- Lam_hat%*%C_hat_inv%*%Lam_hat
  log_map_FM <- log_map(FM,th_n,ph_n)
  #Wald_stat <- n*t(FM - mu_0)%*%V_Wald%*%(FM - mu_0)
  Wald_stat <- n*t(log_map_FM)%*%V_Wald%*%t(t(log_map_FM))
  Wald_save[j] <- Wald_stat
  Sys.sleep(1 / MC)
}

ind_Wald <- which(Wald_save<=5000)
#####################################
#compute size
q1 <- qchisq(0.95,df=dim)
rate <- length(which(Wald_save>q1))/MC
size <- rate
power <- rate
size
power

