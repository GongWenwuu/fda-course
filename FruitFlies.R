rm(list = ls())
library(fda)
library(ggplot2)
load('E://R Files/Functional Data Analysis/Data Sets/fly.Rdata')# medfly(row,col): the number of eggs, the different fruit fly
flyname = colnames(medfly$eggcount)
day = 1:26
y = medfly$eggcount[,flyname[1]]
####################################################################    
plot(day,y,col=2,cex=1.5,xlab='day',ylab='eggcount',ylim=c(0,140),
     main=flyname[1],cex.lab=1.5, cex.axis=1.5)

par(mfrow=c(1,2),ask=T)

for (i in 1:50){
  y = medfly$eggcount[,flyname[i]]
  plot(day,y,col=2,cex=1.5,xlab='day',ylab='y',cex.lab=1.5,
       main=flyname[i],cex.axis=1.5)
}

#################################################################### 
# Find out the order of basis expansion by SSE(OCV)
# Polynomial basis
num1 = 6
P_bvals = matrix(day,26,num1) 
for(i in 1:num1){ P_bvals[,i] = ((P_bvals[,i]-mean(y))/var(y))^(i-1)} 
P_bvals[,1] = 0.01

smsse = rep(0,num1) 
for(i in 1:num1){ 
  S = P_bvals[,1:i]%*%solve( t(P_bvals[,1:i])%*%P_bvals[,1:i])%*%t(P_bvals[,1:i])
  h = (1-diag(S))^2
  errs = y - S%*%y
  smsse[i] = sum(errs^2/h )/26 # leave-one-out CV for linear smoother
}
print(smsse)
which.min(smsse)

ggplot(data=NULL,aes(1:num1,smsse[1:num1]))+geom_point(color = "darkred",size=5)+geom_line(linetype="dotted",size=2)+
  theme_bw()+theme(panel.grid=element_blank())+xlab("No Basis Functions") +ylab("Ordinary Cross-Validated Error")+
  annotate("text",x =4, y = 150,parse = T,label="min(CV) == 133.08")  

# Fourier basis
daybasis = create.fourier.basis(c(0, 26), 26)
bvals = eval.basis(day,daybasis) 
num2 = 12

smsse = rep(0,num2)
for(i in 1:num2){ 
  F_bvals = bvals[,1:(2*i+1)]
  S = F_bvals%*%solve( t(F_bvals)%*%F_bvals)%*%t(F_bvals)
  h = (1-diag(S))^2
  errs = y - S%*%y
  smsse[i] = sum(errs^2/h )/26
}
print(smsse)
which.min(smsse)
ggplot(data=NULL,aes(1:num2,smsse[1:num2]))+geom_point(color = "darkred",size=5)+geom_line(linetype="dotted",size=2)+
  theme_bw()+theme(panel.grid=element_blank())+xlab("No Basis Functions") +ylab("Ordinary Cross-Validated Error")+
  annotate("text",x =11, y = 180,parse = T,label="min(CV) == 154.29") 

# Spline basis
# knots = cumsum(c(0,5,5,5,5,6)) 
bbasis = list()
for (i in 1:8){
  bbasis[[i]] = create.bspline.basis(rangeval=c(0,26), nbasis = 6+i, norder=4) # nbasis=nbreaks+norder-2
}

smsse = rep(0,8)
for(i in 1:8){ 
  bvals = eval.basis(day,bbasis[[i]])
  S = bvals%*%solve( t(bvals)%*%bvals)%*%t(bvals)
  h = (1-diag(S))^2
  errs = y - S%*%y
  smsse[i] = sum(errs^2/h )/26
}
print(smsse)
which.min(smsse)
ggplot(data=NULL,aes(7:14,smsse[1:8]))+geom_point(color = "darkred",size=5)+geom_line(linetype="dotted",size=2)+
  theme_bw()+theme(panel.grid=element_blank())+xlab("No Basis Functions") +ylab("Ordinary Cross-Validated Error")+
  annotate("text",x =12, y = 280,parse = T,label="min(CV) == 154.43") 

# Least-square estimate for fixed num basis functions
par(mfrow=c(1,1), xpd=NA, bty="n")

pbasis = create.power.basis(c(0,26),4,seq(0,3,1)) 
P_dayEggfd = smooth.basis(day, y, pbasis)
plotfit.fd(y, day, P_dayEggfd$fd, residual = FALSE, titles = "512",
           xlab="Day", ylab= "eggcount", cex.lab=1.5, cex.axis=1.5)
plotfit.fd(y, day, P_dayEggfd$fd, residual = TRUE,titles = "RMSE = 9.85",
           xlab="Day", ylab= "eggcount", cex.lab=1.5, cex.axis=1.5)
RMSE1 = sqrt(mean((eval.fd(day, P_dayEggfd$fd) - y)^2))

fbasis = create.fourier.basis(c(0, 26), 11) # 11 basis functions
F_dayEggfd = smooth.basis(day, y, fbasis)
plotfit.fd(y, day, F_dayEggfd$fd, residual = FALSE, titles = '512',
           xlab="Day", ylab= "eggcount", cex.lab=1.5, cex.axis=1.5)
plotfit.fd(y, day, F_dayEggfd$fd, residual = TRUE, titles = 'RMSE = 8.38',
           xlab="Day", ylab= "eggcount", cex.lab=1.5, cex.axis=1.5)
RMSE2 = sqrt(mean((eval.fd(day, F_dayEggfd$fd) - y)^2))

B_dayEggfd = smooth.basis(day, y, bbasis[[6]]) # 8 knots, 12 basis functions
plotfit.fd(y, day, B_dayEggfd$fd, residual = FALSE, titles = '512',
           xlab="Day", ylab= "eggcount", cex.lab=1.5, cex.axis=1.5)
plotfit.fd(y, day, B_dayEggfd$fd, residual = TRUE, titles = 'RMSE = 7.25',
           xlab="Day", ylab= "eggcount", cex.lab=1.5, cex.axis=1.5)
RMSE3 = sqrt(mean((eval.fd(day, B_dayEggfd$fd) - y)^2))

# Penalty estimate
bbasis_F = create.bspline.basis(rangeval=c(0,26), nbasis = 28, norder=4)
gcv = rep(0,7)
df = rep(0,7)
sse = rep(0,7)
Lfd = int2Lfd(2) # We can use different Lfd(nderiv,a list of fd object)
for(i in 1:7){
  lambda=10^{i-4}
  D2fdPar = fdPar(bbasis_F,Lfdobj=Lfd,lambda=lambda) 
  P_dayeggfd = smooth.basis(day,y,D2fdPar)
  gcv[i] = P_dayeggfd$gcv
  df[i] = P_dayeggfd$df
  sse[i] = P_dayeggfd$SSE
  mainstr = paste('lambda = ',lambda,' df = ',df[i],' gcv = ',gcv[i],sep='')
  plotfit.fd(y,day,P_dayeggfd$fd,title=mainstr, 
             xlab="Day", ylab= "eggcount",cex.lab=1.5, cex.axis=1.5)
}
par(mfrow=c(1,1),ask=F)
lambda = 10^(-3:3)
res = data.frame(lambda,gcv,df,sse)
# res = melt(res,id="lambda")
# ggplot(data=res,aes(log(lambda),value,colour=variable))+geom_line(linetype="dotted")+
#   theme_bw()+theme(panel.grid=element_blank())

ggplot(data=res,aes(log(lambda),gcv))+geom_point(color = "darkred",size=5)+geom_line(linetype="dotted",size=2)+
  theme_bw()+theme(panel.grid=element_blank())+
  annotate("text",x =2.4, y = 150,parse = T,label="lamda == 10")
ggplot(data=res,aes(log(lambda),df))+geom_point(color = "darkred",size=5)+geom_line(linetype="dotted",size=2)+
  theme_bw()+theme(panel.grid=element_blank())
ggplot(data=res,aes(log(lambda),sse))+geom_point(color = "darkred",size=5)+geom_line(linetype="dotted",size=2)+
  theme_bw()+theme(panel.grid=element_blank())
which.min(gcv)

par(mfrow=c(1,1), xpd=NA, bty="n")
D2fdPar = fdPar(bbasis_F,Lfdobj=Lfd,lambda=lambda[5]) # functional parameter par
P_dayEggfd = smooth.basis(day,y,D2fdPar)
plotfit.fd(y,day,P_dayEggfd$fd,title="512,lambda=10", 
           xlab="Day", ylab= "eggcount",cex.lab=1.5, cex.axis=1.5)
plotfit.fd(y, day, P_dayEggfd$fd, residual = TRUE, titles = 'RMSE = 9.09',
           xlab="Day", ylab= "eggcount", cex.lab=1.5, cex.axis=1.5)
RMSE = sqrt(mean((eval.fd(day, P_dayEggfd$fd) - y)^2))


 