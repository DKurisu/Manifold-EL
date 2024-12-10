library(MASS)
library(progress)
library(circular)

turtle <- c(8,9,13,13,14,18,22,27,30,34,
            38,38,40,44,45,47,48,48,48,48,
            50,53,56,57,58,58,61,63,64,64,
            64,65,65,68,70,73,78,78,78,83,
            83,88,88,88,90,92,92,93,95,96,
            98,100,103,106,113,118,138,153,153,155,
            204,215,223,226,237,238,243,244,250,251,
            257,268,285,319,343,350)
turtle <- turtle*pi/180
n <- length(turtle)
t_direction <- array(0,dim=c(n,2))
for(i in 1:n){
  t_direction[i,1] <- cos(turtle[i])
  t_direction[i,2] <- sin(turtle[i])
}
########################################
#mean direction
cm <- sum(t_direction[,1])
sm <- sum(t_direction[,2])
t_mean_d <- coord2rad(cm,sm)
t_mean_d <- c(cos(1.120001),sin(1.120001))
t_mean_d
##########################################
#Compute Frechet mean
Fre_cand <- 128
Fre <- rep(0,Fre_cand)
ph_cand <- (seq(1:Fre_cand)-1)/Fre_cand*(2*pi)
for(i in 1:Fre_cand){
  m1 <- cos(ph_cand[i])
  m2 <- sin(ph_cand[i])
  m <- c(m1,m2)
  A <- 0
  for(k in 1:n){
    A <- A + acos(sum(m*t_direction[k,]))^2
  }
  Fre[i] <- A/n
}
##################################################
xlabels <- c(expression(paste(0,pi,"/",128)),
             expression(paste(40,pi,"/",128)),
             expression(paste(80,pi,"/",128)),
             expression(paste(120,pi,"/",128)),
             expression(paste(160,pi,"/",128)),
             expression(paste(200,pi,"/",128)),
             expression(paste(240,pi,"/",128)))
xvalues <- c(1:13-1)*10
xvalues2 <- c(1:7-1)*20
plot(Fre,type="l",xaxt="n",xlab=expression(theta),ylab="Frechet function")
axis(side=1,at=xvalues,labels=F,tcl=-0.25) 
axis(side=1,at=xvalues2,labels=xlabels) 
ind_FM <- which(Fre == min(Fre))
FM <- c(cos(ph_cand[ind_FM]),sin(ph_cand[ind_FM]))
FM
##########################################
#Compute Frechet median
F_med_cand <- 128
F_med <- rep(0,F_med_cand)
ph_cand <- (seq(1:F_med_cand)-1)/F_med_cand*(2*pi)
for(i in 1:F_med_cand){
  m1 <- cos(ph_cand[i])
  m2 <- sin(ph_cand[i])
  m <- c(m1,m2)
  A <- 0
  for(k in 1:n){
    A <- A + acos(sum(m*t_direction[k,]))
  }
  F_med[i] <- A/n
}

plot(F_med,type="l",xaxt="n",xlab=expression(theta),ylab="Frechet function")
axis(side=1,at=xvalues,labels=F,tcl=-0.25) 
axis(side=1,at=xvalues2,labels=xlabels) 
ind_F_med <- which(F_med == min(F_med))
Fmed <- c(cos(ph_cand[ind_F_med]),sin(ph_cand[ind_F_med]))
Fmed
######################
#g function
g_func <- function(x,ph){#x: vector, th: theta, ph: phi
  m1 <- cos(ph)
  m2 <- sin(ph)
  m <- c(m1,m2)
  #g1 <- x[1]*cos(ph) + x[2]*sin(ph)
  g1 <- -x[1]*sin(ph) + x[2]*cos(ph)
  A <- acos(sum(m*x))
  B <- sqrt(1-(sum(m*x))^2)
  C <- -2*sin(1)*(A/B)
  D <- C*g1
  return(D)
}
g_func2 <- function(x,ph){#x: vector, th: theta, ph: phi
  m1 <- cos(ph)
  m2 <- sin(ph)
  m <- c(m1,m2)
  g1 <- x[1]*cos(ph) + x[2]*sin(ph)
  C <- -2*sin(1)
  D <- C*g1
  return(D)
}
#C-matrix
C_mat <- function(ph){
  CM <- array(0,dim=c(2,2))
  CM[1,1] <- cos(ph);CM[1,2] <- -sin(ph)
  CM[2,1] <- sin(ph);CM[2,2] <-  cos(ph)
  return(CM)
}
######################
#g function for median
g_func_med <- function(x,ph){#x: vector, ph: phi
  m1 <- cos(ph)
  m2 <- sin(ph)
  m <- c(m1,m2)
  A <- x - sum(m*x)*m
  B <- sqrt(sum(A^2))
  C <- A/B
  D <- t(C_mat(ph))%*%C
  E <- D[2]
  return(E)
}
#######################
ph_vec <- (seq(1:128)-1)/128*(2*pi)
ph_N <- length(ph_vec)
EL_dual_opt <- array(0,dim=c(ph_N))
EL_dual_opt_med <- array(0,dim=c(ph_N))
#################################
g_mat <- array(0,dim=c(n,ph_N))
g_mat_med <- array(0,dim=c(n,ph_N))
for(i in 1:n){
  for(k in 1:ph_N){
    g_mat[i,k] <- g_func(t_direction[i,],ph_vec[k])
    g_mat_med[i,k] <- g_func_med(t_direction[i,],ph_vec[k])
  }
}

g_mat[15,17] <- g_func2(t_direction[15,],ph_vec[17])
g_mat[45,33] <- g_func2(t_direction[45,],ph_vec[33])
g_mat[15,81] <- 10000
g_mat[45,97] <- 10000
g_mat_med[15,17] <- 0
g_mat_med[45,33] <- 0
pb <- progress_bar$new(total = ph_N,
                       format = "[:bar] :percent : :eta",
                       clear = TRUE)
############################################
for(j in 1:ph_N){
  pb$tick()
  #Compute dual form of EL
  EL_dual <- function(lam1){
    g_vec2 <- rep(1,n)
    for(i in 1:n){
      g_vec2[i] <- -2*log(1+lam1*g_mat[i,j])
      }
    return(sum(g_vec2))
  }
  EL_dual_med <- function(lam2){
    g_vec2_med <- rep(1,n)
    for(i in 1:n){
      g_vec2_med[i] <- -2*log(1+lam2*g_mat_med[i,j])
    }
    return(sum(g_vec2_med))
  }
  lambda <- optim(0,EL_dual)$par
  lambda_med <- optim(0,EL_dual_med)$par
  #lambda <- nlm(EL_dual,0)
  EL_dual_opt[j] <- -EL_dual(lambda) 
  EL_dual_opt_med[j] <- -EL_dual_med(lambda_med) 
  Sys.sleep(1/ph_N)
  }
########################################
#compute confidence region
sp_point <- function(ph){
  s1 <- cos(ph)
  s2 <- sin(ph)
  s <- c(s1,s2)
  return(s)
}
dim <- 1
q_EL <- qchisq(0.95,df=dim)
#q_med <- 2*log(n)
q_med <- qchisq(0.95,df=dim)
ind_EL <- which(EL_dual_opt<=q_EL)
ind_EL_med <- which(EL_dual_opt_med<=q_med)
cr_point <- length(ind_EL)
CR_point <- array(0,dim=c(cr_point,2))
cr_point_med <- length(ind_EL_med)
CR_point_med <- array(0,dim=c(cr_point_med,2))
for(i in 1:cr_point){
  CR_point[i,1] <- cos(ph_vec[ind_EL[i]])
  CR_point[i,2] <- sin(ph_vec[ind_EL[i]])
}
for(i in 1:cr_point_med){
  CR_point_med[i,1] <- cos(ph_vec[ind_EL_med[i]])
  CR_point_med[i,2] <- sin(ph_vec[ind_EL_med[i]])
}

##################################################
xlabels <- c(expression(paste(0,pi,"/",128)),
             expression(paste(40,pi,"/",128)),
             expression(paste(80,pi,"/",128)),
             expression(paste(120,pi,"/",128)),
             expression(paste(160,pi,"/",128)),
             expression(paste(200,pi,"/",128)),
             expression(paste(240,pi,"/",128)))
xvalues <- c(1:13-1)*10
xvalues2 <- c(1:7-1)*20
##################################################
plot(EL_dual_opt,type="l",xlim=c(0,ph_N),ylim=c(0,50),xaxt="n",xlab=expression(theta),ylab="EL")
axis(side=1,at=xvalues,labels=F,tcl=-0.25) 
axis(side=1,at=xvalues2,labels=xlabels) 
abline(h=q_EL,col="red")
abline(v=ind_FM,col="red")
##################################################
plot(EL_dual_opt_med,type="l",xlim=c(0,ph_N),ylim=c(0,50),xaxt="n",xlab=expression(theta),ylab="EL")
axis(side=1,at=xvalues,labels=F,tcl=-0.25) 
axis(side=1,at=xvalues2,labels=xlabels) 
abline(h=q_med,col="purple")
abline(v=ind_F_med,col="purple")