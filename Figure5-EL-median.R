##########################################
# Manifold EL code iterate for smeary case
##########################################
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

g_func_med2 <- function(x,th,ph){#x: vector, th: theta, ph: phi
  m1 <- sin(th)*cos(ph)
  m2 <- sin(th)*sin(ph)
  m3 <- cos(th)
  m <- c(m1,m2,m3)
  g1 <- x[1]*cos(th)*cos(ph) + x[2]*cos(th)*sin(ph) - x[3]*sin(th)
  g2 <- -x[1]*sin(ph) + x[2]*cos(ph)
  A <- sum(m*x)
  B <- sqrt(1-(sum(m*x))^2)
  C <- 2*sin(1)*A/B
  D <- c(g1,g2)*C
  return(D)
}
#######################
#Compute dual form iteratively
MC <-1000
n <- 2000 # sample size
v2 <- 2*pi^{3/2}/gamma(3/2)
v3 <- 2*pi^{4/2}/gamma(4/2)
ga <- v3/(2*v2)
p <- 1-1/(1+ga) # success rate of Bernoulli r.v.s
EL_dual_opt <- rep(0,MC)
pb <- progress_bar$new(total = MC,
                       format = "[:bar] :percent : :eta",
                       clear = TRUE)
theta <- 0 #null theta
phi <- 0   #null phi
for(j in 1:MC){
  pb$tick()
  ##########################################
  #Generate data 
  ###################################
  #Generate Bernolli r.v.s
  be <- rbinom(n,1,p)
  ind_lhs <- which(be==0)
  ind_north <- which(be==1)
  n_lhs <- length(ind_lhs)
  n_north <- length(ind_north)
  ####################################
  #uniform dist on lower-half sphere
  rt <- acos(runif(n_lhs,min=-1,max=0)) #angle for z-axis
  rp <- runif(n_lhs,min=0,max=(2*pi)) #angle for xy-plane
  r_sphere <- array(dim=c(n_lhs,3))
  for(i in 1:n_lhs){
    r_sphere[i,1] <- sin(rt[i])*sin(rp[i])
    r_sphere[i,2] <- sin(rt[i])*cos(rp[i])
    r_sphere[i,3] <- cos(rt[i])
  }
  ###################################
  #generate two smeary data
  X <- array(0,dim=c(n,3))
  X[ind_north,1] <- 0
  X[ind_north,2] <- 0
  X[ind_north,3] <- 1
  for(i in 1:n_lhs){
    X[ind_lhs[i],] <- r_sphere[i,]
  }
  #######################
  #Compute dual form of EL
  EL_dual <- function(lam){
    g_vec2 <- rep(0,n)
    for(i in 1:n){
      g_vec2[i] <- -2*log(1+sum(lam*g_func_med(X[i,],theta,phi)))
    }
    #when null theta = 0  
    #for(i in 1:n_lhs){
    # g_vec2[ind_lhs[i]] <- -2*log(1+sum(lam*g_func_med(X[ind_lhs[i],],0,phi)))
    #}
    #for(i in 1:n_north){
    # #g_vec2[ind_north[i]] <- -2*log(1+sum(lam*g_func_med(X[ind_north[i],],1e-16,phi)))
    #  g_vec2[ind_north[i]] <- -2*log(1+sum(lam*c(0,0)))
    #}
    #when null theta = pi/2 
    #for(i in 1:n_lhs){
    #  g_vec2[ind_lhs[i]] <- -2*log(1+sum(lam*g_func_med(X[ind_lhs[i],],0,phi)))
    #}
    #for(i in 1:n_north){
    #  g_vec2[ind_north[i]] <- -2*log(1+sum(lam*c(1,0)))
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

