######################
# Manifold EL code iterate
######################
library(rgl) #3D plot
library(rotasym) #Generate von Mises-Fisher dist.
library(RiemBase) #Compute Frechet mean for spherical data
library(progress)
######################
#Frechet median
dim <- 2
theta <- 0
phi <- 0
mu1 <- sin(theta)*cos(phi)
mu2 <- sin(theta)*sin(phi)
mu3 <- cos(theta)
mu <- c(mu1,mu2,mu3)
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
#######################
#Compute dual form iteratively
MC <- 1000
kappa <- 1 #parameter of concentration of von Mises-Fisher dist.
n <- 200 #sample size
EL_dual_opt <- rep(0,MC)
pb <- progress_bar$new(total = MC,
                       format = "[:bar] :percent : :eta",
                       clear = TRUE)
theta <- 0 #null theta
phi <- 0   #null phi
for(j in 1:MC){
  pb$tick()
  ##########################################
  #Generate data from von Mises-Fisher dist.
  X <- r_vMF(n,mu,kappa)
  #######################
  #Compute dual form of EL
  EL_dual <- function(lam){
    g_vec2 <- rep(0,n)
    for(i in 1:n){
      g_vec2[i] <- -2*log(1+sum(lam*g_func_med(X[i,],theta,phi)))
    }
    return(sum(g_vec2))
    }
  lambda <- optim(c(0,0),EL_dual)$par
  EL_dual_opt[j] <- -EL_dual(lambda) 
  Sys.sleep(1 / MC)
}

#####################################
#compute size
q1 <- qchisq(0.95,df=2)
rate <- length(which(EL_dual_opt>q1))/MC
size <- rate
power <- rate
size
power

