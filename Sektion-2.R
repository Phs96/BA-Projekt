###############OPSÆTNING#################
library(quantmod) # check
###############OPSÆTNING#################

#### Figur 1 - sandsynlighed
###Sandsynlighed plot
p_funktion <- function(t,sigma,r) return(pnorm(-(1/(sqrt(t)*sigma))-((r*sqrt(t))/sigma),0,1))
p_vektor <- Vectorize(function(t, sigma) {p_funktion(t,sigma,0.1)})
sigma_ssh <- seq(0.01, 0.99, 0.005)
deltat_ssh <- seq(0.01, 0.99, 0.005)
ssh <- matrix(NA, nrow = length(deltat_ssh), ncol = length(sigma_ssh))
for (i in 1:length(sigma_ssh)){
  ssh[i,] <- p_vektor(deltat_ssh[i],sigma_ssh)
}
teorip <- persp(sigma_ssh, deltat_ssh, ssh, theta = 35, xlab = "Sigma", ylab ="T", zlab = "p",expand = 0.5 ,col="blue",
                ticktype = "detailed", phi = 0, zlim = c(0,0.15), border = NA)
points(trans3d(0.1,1/10,p_funktion(1/10,0.1,0.03),pmat = teorip), col = "red", pch = 16)

#### Test af fordeling
sp100 <- new.env()
getSymbols("^OEX", env = sp100, src = "yahoo", from = as.Date("1982-08-02"), to = as.Date("2017-08-02"))
sp100 <- sp100$OEX
S <- as.vector(sp100$OEX.Close)

Afkast <- function(S){
  Safkast <- c()
  for (i in 2:length(S)){
    Safkast[i-1] <- log(S[i]/S[i-1])
  }
  return(Safkast)
}
#Figur 2
par(mfrow = c(1, 2))
Saf <- Afkast(S)
h <- hist(Saf, breaks = 50, xlab = "Log Afkast", main = NA)
points(min(Saf),0, pch = 16, col = "red", cex = 1.5)
points(min(Saf[Saf != min(Saf)]),0, pch = 16, col = "blue", cex = 1.5)
h <- hist(Saf, breaks = 300, xlim = c(-0.05,0.05), xlab = "Log Afkast", main = NA)
xfit<-seq(min(Saf),max(Saf),length=length(Saf))
yfit<-dnorm(xfit,mean=mean(Saf),sd=sd(Saf))
yfit <- yfit*diff(h$mids[1:2])*length(Saf)
lines(xfit, yfit, col="chartreuse2", lwd=2)
#Datoer
sp100[which.min(Saf)];sp100[which.min(Saf[Saf != min(Saf)])]
par(mfrow = c(1, 1))
####Mu sigma
a <- function(n,S) {return ((log(S[n])-log(S[1]))/n)}
b <- function(n,S){
  b <- 0
  for (i in 2:n){
    b <- b + ((log(S[i]/S[i-1]) - a(n,S))^2)/n
  }
  return (b)
}
t <- 35/length(S) # 
n <- length(S)
sigmahat <- sqrt(b(n,S)/t); sigmahat
muhat <- 1/2*sigmahat^2+a(n,S)/(t); muhat

#### optim loglikelihood
si <- numeric(n-1)
si1 <- numeric(n-1)
for (i in 1:(n-1)){
  si[i] <- S[i+1]
  si1[i] <- S[i] 
}
loglike <- function( pars, si. = si, si1. = si1) {
  mu <- pars[1]
  sigma2 <- pars[2]^2
  val.log.like <- n/2*log(2*pi)+n/2*log(sigma2*t)+1/(2*sigma2*t)*sum((log(si./si1.)-(mu - 1/2*sigma2)*t)^2)
  return(val.log.like )
}
optimization <- optim( c( 0.1, 0.5), loglike)
answer <- matrix( optimization$par, 2, 1) 
answer #estimatorne

###mu sigma 2010
S2010 <- as.vector(sp100$OEX.Close['2010'])
n2010 <- length(S2010)
t2010 <- 1/n2010
sigmahat2010 <- b(n2010,S2010)/t2010; sqrt(sigmahat2010)
mu2010 <- 1/2*sigmahat2010+a(n2010,S2010)/(t2010); mu2010
