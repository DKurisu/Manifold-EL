######################
# Manifold EL code iterate
####################################
library(rgl) #3D plot
library(rotasym) #Generate von Mises-Fisher dist.
library(RiemBase) #Compute Frechet mean for spherical data
library(progress)
library(shapes)
###################################
load("T1mice.rda")
T1mice$x #print data
T1mice$group #print group label
which(T1mice$group=="c")
nc <- length(which(T1mice$group=="c"))
nl <- length(which(T1mice$group=="l"))
ns <- length(which(T1mice$group=="s"))
vec_mice <- c(1,2,3)
lg <- T1mice$x[vec_mice,,(1:nl)+nc]
sg <- T1mice$x[vec_mice,,(1:ns)+nc+nl]
k <- length(lg[,1,1]) #number of landmarks
m <- length(lg[1,,1]) #dimension of the space where landmarks live


#g function
g_func <- function(x,th,ph){#x: vector, th: theta, ph: phi
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  thx <- acos(x[3])
  phx <- sign(x[2])*acos(x[1]/sqrt(x[1]^2+x[2]^2))
  x1 <- sin(thx)*cos(phx)
  x2 <- sin(thx)*sin(phx)
  x3 <- cos(thx)
  m <- c(m1,m2,m3)
  g1 <- x1*cos(th)*cos(ph) + x2*cos(th)*sin(ph) - x3*sin(th)
  g2 <- -x1*sin(ph) + x2*cos(ph)
  A <- -2*sin(1)*acos(sum(m*x))
  B <- 1/sqrt(1-(sum(m*x))^2)
  C <- B*c(g1,g2)
  D <- A*C
  return(D)
}
#g function-A
g_func_A <- function(x,th,ph){#x: vector, th: theta, ph: phi
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  m <- c(m1,m2,m3)
  A <- acos(sum(m*x))
  return(A)
}
g_func_C <- function(x,th,ph){#x: vector, th: theta, ph: phi
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  thx <- acos(x[3])
  phx <- sign(x[2])*acos(x[1]/sqrt(x[1]^2+x[2]^2))
  x1 <- sin(thx)*cos(phx)
  x2 <- sin(thx)*sin(phx)
  x3 <- cos(thx)
  m <- c(m1,m2,m3)
  g1 <- x1*cos(th)*cos(ph) + x2*cos(th)*sin(ph) - x3*sin(th)
  g2 <- -x1*sin(ph) + x2*cos(ph)
  A <- -2*sin(1)*acos(sum(m*x))
  B <- 1/sqrt(1-(sum(m*x))^2)
  C <- B*c(g1,g2)
  return(C)
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


###########################################
#mapping oliginal data to pre-shape space
#centering
#large group
Lgc <- array(0,dim=c(k,2,nl))
Lg <- array(0,dim=c(k,2,nl))
for(i in 1:nl){
  lgc <- apply(lg[,,i],2,mean.default) #apply mean w.r.t row (1) or column (2)
  Lgc[1:3,1,i] <- lgc[1]
  Lgc[1:3,2,i] <- lgc[2]
  Lg[,,i] <- lg[,,i] - Lgc[,,i] #centering large group
}
#small group
Sgc <- array(0,dim=c(k,2,ns))
Sg <- array(0,dim=c(k,2,ns))
for(i in 1:ns){
  sgc <- apply(sg[,,i],2,mean.default) #apply mean w.r.t row (1) or column (2)
  Sgc[1:3,1,i] <- sgc[1]
  Sgc[1:3,2,i] <- sgc[2]
  Sg[,,i] <- sg[,,i] - Sgc[,,i] #centering large group
}
###############################
# mapping pre-shape space
k1 <- k-1 
vec_pre <- c(2,k)
#large group
Lgps1 <- Lg[vec_pre,,]
Lgd <- sqrt(apply(Lgps1[,1,]^2,2,sum)+apply(Lgps1[,2,]^2,2,sum)) #compute Euclidean distance for each data
Lgps2 <- array(0,dim=c(k1,m,nl))
for(i in 1:nl){
  Lgps2[,,i] <- Lgps1[,,i]/Lgd[i]
}
#small group
Sgps1 <- Sg[vec_pre,,]
Sgd <- sqrt(apply(Sgps1[,1,]^2,2,sum)+apply(Sgps1[,2,]^2,2,sum)) #compute Euclidean distance for each data
Sgps2 <- array(0,dim=c(k1,m,ns))
for(i in 1:ns){
  Sgps2[,,i] <- Sgps1[,,i]/Sgd[i]
}
####################################
#rotatation
#large group
Lgps4 <- array(0,dim=c(k1,m,nl))
ind_Lgps2_n<-which(Lgps2[1,2,]<0)
Lgps2_2 <- Lgps2
Lgps2_2[,2,ind_Lgps2_n] <- -Lgps2[,2,ind_Lgps2_n]
for(i in 1:nl){
  Lgps3 <- Lgps2_2[1,,i]/sqrt(sum(Lgps2_2[1,,i]^2))
  th_ps1 <- acos(Lgps3[1])
  Rot_L <- matrix(c(cos(-th_ps1),sin(-th_ps1),-sin(-th_ps1),cos(-th_ps1)),ncol=2)
  for(j in 1:k1){
    Lgps4[j,,i] <- t(Rot_L%*%t(t(Lgps2_2[j,,i])))
  }
}
Lgps4_2 <- Lgps4
Lgps4_2[,2,ind_Lgps2_n] <- -Lgps4[,2,ind_Lgps2_n]
#small group
Sgps4 <- array(0,dim=c(k1,m,ns))
ind_Sgps2_n<-which(Sgps2[1,2,]<0)
Sgps2_2 <- Sgps2
Sgps2_2[,2,ind_Sgps2_n] <- -Sgps2[,2,ind_Sgps2_n]
for(i in 1:ns){
  Sgps3 <- Sgps2_2[1,,i]/sqrt(sum(Sgps2_2[1,,i]^2))
  th_ps2 <- acos(Sgps3[1])
  Rot_S <- matrix(c(cos(-th_ps2),sin(-th_ps2),-sin(-th_ps2),cos(-th_ps2)),ncol=2)
  for(j in 1:k1){
    Sgps4[j,,i] <- t(Rot_S%*%t(t(Sgps2_2[j,,i])))
  }
}
Sgps4_2 <- Sgps4
Sgps4_2[,2,ind_Sgps2_n] <- -Sgps4[,2,ind_Sgps2_n]
####################################
# compute map from Sigma_2^k to S^{2k-4}
k2 <- 2*k-3
#large group
Lgps5 <- array(0,dim=c(k2,nl))
for(i in 1:nl){
  vec_L1 <- c(Lgps4_2[2,,i],Lgps4_2[1,1,i])
  Lgps5[,i] <- vec_L1/sqrt(sum(vec_L1^2))
}
#small group
Sgps5 <- array(0,dim=c(k2,ns))
for(i in 1:ns){
  vec_S1 <- c(Sgps4_2[2,,i],Sgps4_2[1,1,i])
  Sgps5[,i] <- vec_S1/sqrt(sum(vec_S1^2))
}


#######################
boot <-199 #number of bootstrap replications
X1 <- t(Lgps5);n1 <- length(X1[,1]) #Lgps5
X2 <- t(Sgps5);n2 <- length(X2[,1])
X_merge <- rbind(X1,X2);n_merge <- n1+n2
######################
lam_save_1b <- array(0,dim=c(boot,2))
lam_save_2b <- array(0,dim=c(boot,2))
EL_dual_opt_1b <- rep(0,boot)
EL_dual_opt_2b <- rep(0,boot)
EL_dual_opt_3b <- rep(0,boot)
quant_TS <- rep(0,3)
  ##########################################
  #Compute Frechet median/mean of marged sample
  th_max <- 128
  ph_max <- 256
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
  ind_Fmean <- which(Fre_mean == min(Fre_mean))
  ind_Fmean2 <- floor(ind_Fmean/ph_max)
  ind_Fmean3 <- ind_Fmean-ph_max*ind_Fmean2
  ind_Fmean_th <- ind_Fmean2+1
  ind_Fmean_ph <- ind_Fmean3+1
  theta_m <- min(th_cand[ind_Fmean_th])  #theta of the Frechet mean computed from the marged sample
  phi_m <- min(ph_cand[ind_Fmean_ph])   #phi of the Frechet mean computed from the marged sample
  #######################
  #Compute dual form of EL
  #first sample
  EL_dual_1 <- function(lam){
    g_vec2 <- rep(0,n1)
    #when null theta neq pi
    for(w1 in 1:n1){
      g_vec2[w1] <- -2*log(1+sum(lam*g_func(X1[w1,],theta_m,phi_m)))
    }
    #when null theta = pi
    #for(i in 1:n1){
    #  g_vec2[i] <- -2*log(1+sum(lam*g_func3(X1[w1,])))
    #}
    return(sum(g_vec2))
  }
  #second sample
  EL_dual_2 <- function(lam){
    g_vec22 <- rep(0,n2)
    #when null theta neq pi
    for(w2 in 1:n2){
      g_vec22[w2] <- -2*log(1+sum(lam*g_func(X2[w2,],theta_m,phi_m)))
    }
    #when null theta = pi
    #for(w2 in 1:n2){
    #  g_vec2e[w2] <- -2*log(1+sum(lam*g_func3(X2[w2,])))
    #}
    return(sum(g_vec22))
  }
  lambda01 <- optim(c(0,0),EL_dual_1)$par
  lambda02 <- optim(c(0,0),EL_dual_2)$par
  EL_dual_opt_1 <- -EL_dual_1(lambda01) 
  EL_dual_opt_2 <- -EL_dual_2(lambda02) 
  EL_dual_opt_3 <- EL_dual_opt_1 + EL_dual_opt_2
  ######################################
  #sample index for bootstrap
  ind_boot_1 <- sample(c(1:n_merge),size=n1*boot,replace=TRUE)
  ind_boot_2 <- sample(c(1:n_merge),size=n2*boot,replace=TRUE)
  m_boot_1 <- matrix(ind_boot_1,ncol=boot)
  m_boot_2 <- matrix(ind_boot_2,ncol=boot)
  ####################################
  pb <- progress_bar$new(total = boot,
                         format = "[:bar] :percent : :eta",
                         clear = TRUE)
  for(bt in 1:boot){
    pb$tick()
    ##########################################
    X1b <- X_merge[m_boot_1[,bt],]
    X2b <- X_merge[m_boot_2[,bt],]
    #######################
    #Compute dual form of EL
    #first sample
    EL_dual_1b <- function(lam){
      g_vec2_1b <- rep(0,n1)
      #when null theta neq pi
      for(i1 in 1:n1){
        g_vec2_1b[i1] <- -2*log(1+sum(lam*g_func(X1b[i1,],theta_m,phi_m)))
      }
      #when null theta = pi
      #for(i1 in 1:n){
      #  g_vec2_1b[i1] <- -2*log(1+sum(lam*g_func3(Xb[i1,])))
      #}
      return(sum(g_vec2_1b))
    }
    #second sample
    EL_dual_2b <- function(lam){
      g_vec22_2b <- rep(0,n2)
      #when null theta neq pi
      for(i2 in 1:n2){
        g_vec22_2b[i2] <- -2*log(1+sum(lam*g_func(X2b[i2,],theta_m,phi_m)))
      }
      #when null theta = pi
      #for(i2 in 1:n){
      #  g_vec22_2b[i2] <- -2*log(1+sum(lam*g_func3(X2b[i2,])))
      #}
      return(sum(g_vec22_2b))
    }
    lambda_1b <- optim(c(0,0),EL_dual_1b)$par
    lambda_2b <- optim(c(0,0),EL_dual_2b)$par
    EL_dual_opt_1b[bt] <- -EL_dual_1b(lambda_1b) 
    EL_dual_opt_2b[bt] <- -EL_dual_2b(lambda_2b) 
    EL_dual_opt_3b[bt] <- EL_dual_opt_1b[bt] + EL_dual_opt_2b[bt]
    lam_save_1b[bt,] <- lambda_1b
    lam_save_2b[bt,] <- lambda_2b
    Sys.sleep(1 / boot)
  }
  ind_quant <- which(EL_dual_opt_3b<=5000)
  quant1 <- as.numeric(quantile(EL_dual_opt_3b[ind_quant],c(0.90,0.95,0.99)))
  quant1
  EL_dual_opt_3
  