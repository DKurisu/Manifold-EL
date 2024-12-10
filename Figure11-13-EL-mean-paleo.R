####################################
# Manifold EL confidence region code
####################################
library(rgl) #3D plot
library(scatterplot3d) #3D plot
library(rotasym) #Generate von Mises-Fisher dist.
library(RiemBase) #Compute Frechet mean for spherical data
library(progress)
library(lattice)
library(latticeExtra)
######################
VGP <- read.csv("vgps_filtered_v1.csv",header = T)
VGP_la <- VGP$VGP_lat
VGP_lo <- VGP$VGP_lon
ind_VGP_la1 <- which(VGP_la <= 0)
ind_VGP_la2 <- which(VGP_la > 0)
VGP_la1 <- VGP_la
VGP_la1[ind_VGP_la1] <- 90-VGP_la[ind_VGP_la1] 
VGP_la1[ind_VGP_la2] <- 90-VGP_la[ind_VGP_la2] 
ind_VGP_lo <- which(VGP_lo <= 0)
VGP_lo1<- VGP_lo
VGP_lo1[ind_VGP_lo] <- 360+VGP_lo[ind_VGP_lo] 
VGP_th <- VGP_la1/180*pi
VGP_ph <- VGP_lo1/180*pi
n <- length(VGP_th)
##############################
X <- array(0,dim=c(n,3))
for(i in 1:n){
  X[i,1] <- sin(VGP_th[i])*cos(VGP_ph[i])
  X[i,2] <- sin(VGP_th[i])*sin(VGP_ph[i])
  X[i,3] <- cos(VGP_th[i])
}
##########################################
#Compute Frechet median/mean
th_max <- 128
ph_max <- 128
Fre_med <- rep(0,th_max*ph_max)
Fre_mean <- rep(0,th_max*ph_max)
th_cand <- (seq(1:th_max)-1)/(th_max)*pi 
ph_cand <- (seq(1:ph_max)-1)/ph_max*(2*pi)
pb <- progress_bar$new(total = ph_max,
                       format = "[:bar] :percent : :eta",
                       clear = TRUE)
for(i in 1:ph_max){
  pb$tick()
  for(k in 1:th_max){
    mu1 <- sin(th_cand[k])*cos(ph_cand[i])
    mu2 <- sin(th_cand[k])*sin(ph_cand[i])
    mu3 <- cos(th_cand[k])
    mu <- c(mu1,mu2,mu3)
    A <- 0
    B <- 0
    for(j in 1:n){
      A <- A + acos(sum(mu*X[j,]))
      B <- B + acos(sum(mu*X[j,]))^2
    }
    Fre_med[(k-1)*ph_max+i] <- A/n
    Fre_mean[(k-1)*ph_max+i] <- B/n
  }
  Sys.sleep(1/th_max)
}
ind_Fmed <- which(Fre_med == min(Fre_med))
ind_Fmean <- which(Fre_mean == min(Fre_mean))
ind_Fmed2 <- floor(ind_Fmed/ph_max)
ind_Fmed3 <- ind_Fmed-ph_max*ind_Fmed2
ind_Fmean2 <- floor(ind_Fmean/ph_max)
ind_Fmean3 <- ind_Fmean-ph_max*ind_Fmean2
ind_Fmed_th <- ind_Fmed2+1
ind_Fmed_ph <- ind_Fmed3
ind_Fmean_th <- ind_Fmean2+1
ind_Fmean_ph <- ind_Fmean3
F_med <- c(sin(th_cand[ind_Fmed_th])*cos(ph_cand[ind_Fmed_ph]),
           sin(th_cand[ind_Fmed_th])*sin(ph_cand[ind_Fmed_ph]),
           cos(th_cand[ind_Fmed_th]))
F_med
F_mean <- c(sin(th_cand[ind_Fmean_th])*cos(ph_cand[ind_Fmean_ph]),
            sin(th_cand[ind_Fmean_th])*sin(ph_cand[ind_Fmean_ph]),
            cos(th_cand[ind_Fmean_th]))
F_mean
##################################################
#heat map of sample Frechet function
#################################################
th_plot <- c(1:128-1)/128*pi
ph_plot <- (seq(1:128)-1)/128*(2*pi)
Fmean_grid <- expand.grid(x=ph_plot,y=th_plot)
Fmed_grid <- expand.grid(x=ph_plot,y=th_plot)
Fmean_surface <- rep(0,128*128)
Fmed_surface <- rep(0,128*128)
for(i in 1:128){
  for(j in 1:128){
    Fmean_surface[(j-1)*128+i] <- Fre_mean[(j-1)*128+i]
    Fmed_surface[(j-1)*128+i] <- Fre_med[(j-1)*128+i]
  }
}
Fmean_grid$z <- Fmean_surface
Fmed_grid$z <- Fmed_surface
brks <- seq(min(Fre_mean),max(Fre_mean),by=0.05)
fig <- levelplot(z~x*y,data=Fmean_grid, asp=1,col.regions=terrain.colors,xlab=expression(xi),ylab=expression(theta))#,at=brks)
fig + as.layer(contourplot(z~x*y, data=Fmean_grid, asp=1, at=brks, labels=F))
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
##########################################
#Compute EL statistics to construct confidence region
th_vec <- c(1:64-33)/64*pi
ph_vec <- (seq(1:128)-1)/128*(2*pi)
th_N <- length(th_vec)
ph_N <- length(ph_vec)
EL_dual_opt <- array(0,dim=c(th_N,ph_N))
pb <- progress_bar$new(total = ph_N,
                       format = "[:bar] :percent : :eta",
                       clear = TRUE)
for(j in 1:ph_N){
  pb$tick()
  for(k in 1:th_N){
    #Compute dual form of EL
    EL_dual <- function(lam){
      g_vec2 <- rep(0,n)
      for(i in 1:n){
        g_vec2[i] <- -2*log(1+sum(lam*g_func(X[i,],th_vec[k],ph_vec[j])))
      }
      return(sum(g_vec2))
    }
    lambda <- optim(c(0,0),EL_dual)$par
    EL_dual_opt[k,j] <- -EL_dual(lambda) 
  }
  Sys.sleep(1/ph_N)
}
######################
th_plot <- c(1:64-1)/64*pi
ph_plot <- (seq(1:128)-1)/128*(2*pi)
EL_grid <- expand.grid(x=ph_plot,y=th_plot)
EL_surface <- rep(0,64*128)
for(i in 1:128){
  for(j in 1:64){
    EL_surface[(j-1)*128+i] <- EL_dual_opt[j,i]
  }
}
EL_grid$z <- EL_surface
th_mean <- th_cand[ind_Fmean_th]
ph_mean <- ph_cand[ind_Fmean_ph]
##################################################
#heat map of EL statistics
#################################################
brks <- seq(min(EL_grid$z),max(EL_grid$z),by=20)
fig <- levelplot(z~x*y,data=EL_grid, asp=1,col.regions=terrain.colors,xlab=expression(xi),ylab=expression(theta))#,xlim=c(0,2*pi),ylim=c(0,pi))#,at=brks)
fig + as.layer(contourplot(z~x*y, data=EL_grid, asp=1, at=brks, labels=F))