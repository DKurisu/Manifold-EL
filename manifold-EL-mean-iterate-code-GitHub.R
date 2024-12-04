######################
# Manifold EL code iterate
######################
library(rgl) #3D plot
library(rotasym) #Generate von Mises-Fisher dist.
library(RiemBase) #Compute Frechet mean for spherical data
library(progress)
######################
#Frechet mean
dim <- 2
theta <- 0
phi <- 0
mu1 <- sin(theta)*cos(phi)
mu2 <- sin(theta)*sin(phi)
mu3 <- cos(theta)
mu <- c(mu1,mu2,mu3)
######################
#g function
g_func <- function(x,th,ph){#x: vector, th: theta, ph: phi
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  m <- c(m1,m2,m3)
  g1 <- x[1]*cos(th)*cos(ph) + x[2]*cos(th)*sin(ph) - x[3]*sin(th)
  g2 <- -x[1]*sin(ph) + x[2]*cos(ph)
  A <- acos(sum(m*x))
  B <- sqrt(1-(sum(m*x))^2)
  C <- -2*sin(1)*(A/B)
  D <- c(g1,g2)*C
  return(D)
}
######################
#g function theta = pi, phi = 0
g_func3 <- function(x,th,ph){#x: vector, th: theta, ph: phi
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
kappa <- 1 #parameter of concentration of von Mises-Fisher dist.
n <- 200 #sample size
EL_dual_opt <- rep(0,MC)
pb <- progress_bar$new(total = MC,
                       #format = "[:bar] :percent : :eta",
                       format = "[:bar] :percent : :elapsed",
                       clear = TRUE)
theta <- 0@#null theta
phi <- 0    #null phi
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
      g_vec2[i] <- -2*log(1+sum(lam*g_func(X[i,],theta,phi)))
    }
    #when null theta = pi
    #for(i in 1:n){
    #  g_vec2[i] <- -2*log(1+sum(lam*g_func3(X[i,],theta,phi)))
    #}
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

