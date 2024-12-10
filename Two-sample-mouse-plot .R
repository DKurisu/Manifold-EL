####################################
# 3 landmark data
####################################
library(shapes)
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

####################################
par(mfrow=c(1,1))
#Large
plot(c(lg[1:3,1,1],lg[1:3,1,1]),c(lg[1:3,2,1],lg[1:3,2,1]),type="o",pch=19,cex=2,lty=1,lwd=1.5,col=rgb(0,0,0,0.3),
     xlim=c(0,250),ylim=c(0,250),xlab="",ylab="")
title("Large")
par(new=T)
for(itp in 1:nl){
  plot(c(lg[1:3,1,itp],lg[1:3,1,itp]),c(lg[1:3,2,itp],lg[1:3,2,itp]),type="o",pch=19,cex=2,lty=1,lwd=1.5,col=rgb(0,0,0,0.3),
       xlim=c(0,250),ylim=c(0,250),xlab="",ylab="",xaxt="n",yaxt="n")
  par(new=T)
}
par(new=F)
#Small
plot(c(sg[1:3,1,1],sg[1:3,1,1]),c(sg[1:3,2,1],sg[1:3,2,1]),type="o",pch=19,cex=2,lty=1,lwd=1.5,col=rgb(0,0,0,0.3),
     xlim=c(0,250),ylim=c(0,250),xlab="",ylab="")
title("Small")
par(new=T)
for(itp in 1:ns){
  plot(c(sg[1:3,1,itp],sg[1:3,1,itp]),c(sg[1:3,2,itp],sg[1:3,2,itp]),type="o",pch=19,cex=2,lty=1,lwd=1.5,col=rgb(0,0,0,0.3),
       xlim=c(0,250),ylim=c(0,250),xlab="",ylab="",xaxt="n",yaxt="n")
  par(new=T)
}
par(new=F)
