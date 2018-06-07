###############OPSÆRNING#################
library(grid) #c
library(gridExtra) #c
library(gridGraphics) # c
library(reshape) # c
library(ggplot2) #c
###############OPSÆRNING#################
##############Nødvendige ting fra sektion-2
sigmahat = 0.1827096
muhat <- 0.1021119
a <- function(n,S) {return ((log(S[n])-log(S[1]))/n)}
b <- function(n,S){
  b <- 0
  for (i in 2:n){
    b <- b + ((log(S[i]/S[i-1]) - a(n,S))^2)/n
  }
  return (b)
}
##############Nødvendige ting fra sektion-2
####simulation plot
##Method A plot
simulateAplot <- function(sigma,mu,n_time,n_sim,T,S0){
  Sa <- matrix(NA, n_time+1,n_sim) #nrow, ncol
  Sa[1,] <- S0
  t <- T/n_time
  for (i in 2:(n_time+1)){
    Sa[i,] <- Sa[i-1,] + mu*Sa[i-1,]*t+Sa[i-1,]*sigma*sqrt(t)*rnorm(n_sim,0,1)
  }
  data_plot <- data.frame('time' =  seq(0, T, by = t), "GBM" = Sa)
  data_plot <- melt(data_plot,  id = c('time'))
  plot_gbm <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable)) +
    theme(legend.position = "none")+ 
    ggtitle("A") + theme(plot.title = element_text(hjust = 0.5))
  return(plot_gbm)
}
## Method B plot
simulateBplot <- function(sigma,mu,n_time,n_sim,T,S0){
  Sb <- matrix(NA, n_time+1,n_sim) #nrow, ncol
  Sb[1,] <- S0
  t <- T/n_time
  for (i in 2:(n_time+1)){
    Sb[i,] <- Sb[i-1,]*exp((mu-0.5*sigma^2)*t+sigma*sqrt(t)*rnorm(n_sim,0,1))
  }
  data_plot <- data.frame('time' =  seq(0, T, by = t), "GBM" = Sb)
  data_plot <- melt(data_plot,  id = c('time'))
  plot_gbm <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable)) +
    theme(legend.position = "none") + 
    ggtitle("B") + theme(plot.title = element_text(hjust = 0.5))
  return(plot_gbm)
}
## Method C plot
simulateCplot <- function(sigma,mu,n_time,n_sim,T,S0){
  Sc <- matrix(NA, n_time+1,n_sim) #nrow, ncol
  Sc[1,] <- S0
  dt <- T/n_time
  t <- 0
  for (i in 2:(n_time+1)){
    t <- t + dt
    Sc[i,] <- Sc[1,]*exp((mu-0.5*sigma^2)*t+sigma*sqrt(t)*rnorm(n_sim,0,1))
  }
  data_plot <- data.frame('time' =  seq(0, T, by = dt), "GBM" = Sc)
  data_plot <- melt(data_plot,  id = c('time'))
  #colnames(data_plot)[which(names(data_plot) == "value")] <- "S-værdi"
  plot_gbm <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable)) +
    theme(legend.position = "none") + 
    ggtitle("C") + theme(plot.title = element_text(hjust = 0.5))
  return(plot_gbm)
}
grid.arrange(
  simulateAplot(0.1,0.03,200,10,1,100),
  simulateBplot(0.1,0.03,200,10,1,100), 
  simulateCplot(0.1,0.03,200,10,1,100), ncol = 3
)

##Diskretionsfejl
simulateA <- function(sigma,mu,n_time,n_sim,T,S0,k){
  Sa <- matrix(NA, n_time+1,n_sim) #nrow, ncol
  Sa[1,] <- S0
  t <- T/n_time
  for (i in 2:(n_time+1)){
    set.seed(k+i)
    Sa[i,] <- Sa[i-1,] + mu*Sa[i-1,]*t+Sa[i-1,]*sigma*sqrt(t)*rnorm(n_sim,0,1)
  }
  return(Sa)
}
simulateB <- function(sigma,mu,n_time,n_sim,T,S0,k){
  Sb <- matrix(NA, n_time+1,n_sim) #nrow, ncol
  Sb[1,] <- S0
  t <- T/n_time
  for (i in 2:(n_time+1)){
    set.seed(k+i)
    Sb[i,] <- Sb[i-1,]*exp((mu-0.5*sigma^2)*t+sigma*sqrt(t)*rnorm(n_sim,0,1))
  }
  return(Sb)
}
f11 <- mean(abs(simulateA(0.2,0.1,100,10000,1,100,1)[101,]-simulateB(0.2,0.1,100,10000,1,100,1)[101,]))
f12 <- mean(abs(simulateA(0.2,0.1,200,10000,1,100,1)[201,]-simulateB(0.2,0.1,200,10000,1,100,1)[201,]))
f13 <- mean(abs(simulateA(0.2,0.1,500,10000,1,100,1)[501,]-simulateB(0.2,0.1,500,10000,1,100,1)[501,]))
f14 <- mean(abs(simulateA(0.2,0.1,1000,10000,1,100,1)[1001,]-simulateB(0.2,0.1,1000,10000,1,100,1)[1001,]))
f15 <- mean(abs(simulateA(0.2,0.1,2000,10000,1,100,1)[2001,]-simulateB(0.2,0.1,2000,10000,1,100,1)[2001,]))
f21 <- mean(abs(simulateA(0.2,0.1,100,10000,1,100,2001)[101,]-simulateB(0.2,0.1,100,10000,1,100,2001)[101,]))
f22 <- mean(abs(simulateA(0.2,0.1,200,10000,1,100,2001)[201,]-simulateB(0.2,0.1,200,10000,1,100,2001)[201,]))
f23 <- mean(abs(simulateA(0.2,0.1,500,10000,1,100,2001)[501,]-simulateB(0.2,0.1,500,10000,1,100,2001)[501,]))
f24 <- mean(abs(simulateA(0.2,0.1,1000,10000,1,100,2001)[1001,]-simulateB(0.2,0.1,1000,10000,1,100,2001)[1001,]))
f25 <- mean(abs(simulateA(0.2,0.1,2000,10000,1,100,2001)[2001,]-simulateB(0.2,0.1,2000,10000,1,100,2001)[2001,]))
f31 <- mean(abs(simulateA(0.2,0.1,100,10000,1,100,3001)[101,]-simulateB(0.2,0.1,100,10000,1,100,3001)[101,]))
f32 <- mean(abs(simulateA(0.2,0.1,200,10000,1,100,3001)[201,]-simulateB(0.2,0.1,200,10000,1,100,3001)[201,]))
f33 <- mean(abs(simulateA(0.2,0.1,500,10000,1,100,3001)[501,]-simulateB(0.2,0.1,500,10000,1,100,3001)[501,]))
f34 <- mean(abs(simulateA(0.2,0.1,1000,10000,1,100,3001)[1001,]-simulateB(0.2,0.1,1000,10000,1,100,3001)[1001,]))
f35 <- mean(abs(simulateA(0.2,0.1,2000,10000,1,100,3001)[2001,]-simulateB(0.2,0.1,2000,10000,1,100,3001)[2001,]))
konvergens <- c(f11,f12,f13,f14,f15)/sqrt(c(1/100,1/200,1/500,1/1000,1/2000))
#konvergens = 2.479564 2.507576 2.511029 2.465562 2.491357

##Asymptotisk fordeling
fordeling <- function(sigma, mu, n_time, n_sim, T, S0){
  Ssim <- simulateB(sigma,mu,n_time,n_sim,T,S0,k=1000) 
  t <- T/n_time
  n <- n_sim
  estimator <- matrix(NA,nrow = 2, ncol =n_sim)
  for (i in 1:n){
    estimator[1,i] <- sqrt(b(n_time,Ssim[,i])/t)                     #sigma
    estimator[2,i] <- 1/2*estimator[1,i]^2+a(n_time,Ssim[,i])/(t)    #mu
  }
  return(estimator)
}
fordeling_sigma_mu <- fordeling(sigma = sigmahat,
                                mu = muhat, n_time = 100, n_sim = 10000, T=1, S0 = 100)

par(mfrow = c(1, 2))

h_mu <- hist(fordeling_sigma_mu[2,], breaks = 30, xlab = "mu")
xfit<-seq(min(fordeling_sigma_mu[2,]),max(fordeling_sigma_mu[2,]),length=40)
yfit<-dnorm(xfit,mean=mean(fordeling_sigma_mu[2,]),sd=sd(fordeling_sigma_mu[2,]))
yfit <- yfit*diff(h_mu$mids[1:2])*length(fordeling_sigma_mu[2,])
lines(xfit, yfit, col="blue", lwd=2)

h_sig <- hist(fordeling_sigma_mu[1,], breaks = 30, xlab = "sigma")
xfit<-seq(min(fordeling_sigma_mu[1,]),max(fordeling_sigma_mu[1,]),length=40)
yfit<-dnorm(xfit,mean=mean(fordeling_sigma_mu[1,]),sd=sd(fordeling_sigma_mu[1,]))
yfit <- yfit*diff(h_sig$mids[1:2])*length(fordeling_sigma_mu[1,])
lines(xfit, yfit, col="blue", lwd=2)

par(mfrow = c(1, 1))

#### Figur 5,6,7
MetodeA <- function(sigma,r,n_time,n_sim,T,S0,K,t){
  EhatV <- matrix(NA, n_time,n_sim) 
  for (k in 1:n_sim){
    print(k)
    z <- matrix(rnorm(n_time*50,0,1),n_time,50)
    S <- matrix(NA,51, n_time)
    S[1,] <- S0 
    for (i in 2:51){
      S[i,] <- S[i-1,]+ r*S[i-1,]*t + S[i-1,]*sigma*sqrt(t)*z[,i-1]
    }
    V <- pmax(K-S[51,],0)
    for (j in 1:n_time){
      EhatV[j,k] <- mean(V[1:j])
    }
  }
  Vhat <- EhatV*exp(-r)
  data_plot <- data.frame('time' =  seq(1, n_time, by = 1), "GBM" = Vhat)
  data_plot <- melt(data_plot,  id = c('time'))
  plot_gbm <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable)) +
    ylim(4.2,4.7) +
    theme(legend.position = "none") +
    xlab("n")  + ggtitle("Metode A")
  return(plot_gbm)
}
##Figur 3.4 (EU put) 
MetodeA_AV <- function(sigma,r,n_time,n_sim,T,S0,K,t){
  EhatVplus <- matrix(NA, n_time,n_sim) 
  EhatVminus <- matrix(NA, n_time,n_sim) 
  for (k in 1:n_sim){
    print(k)
    z <- matrix(rnorm(n_time*50,0,1),n_time,50)
    Splus <- matrix(NA,51, n_time)
    Sminus <- matrix(NA,51, n_time)
    Splus[1,] <- S0 
    Sminus[1,] <- S0 
    for (i in 2:51){
      Splus[i,] <- Splus[i-1,]+ r*Splus[i-1,]*t + Splus[i-1,]*sigma*sqrt(t)*z[,i-1]
      Sminus[i,] <- Sminus[i-1,]+ r*Sminus[i-1,]*t + Sminus[i-1,]*sigma*sqrt(t)*(-z[,i-1])
    }
    Vplus <- pmax(K-Splus[51,],0)
    Vminus <- pmax(K-Sminus[51,],0)
    for (j in 1:n_time){
      EhatVplus[j,k] <- mean(Vplus[1:j]) 
      EhatVminus[j,k] <- mean(Vminus[1:j]) 
    }
  }
  Vhat <- (EhatVplus*exp(-r) + EhatVminus*exp(-r))/2
  data_plot <- data.frame('time' =  seq(1, n_time, by = 1), "GBM" = Vhat)
  data_plot <- melt(data_plot,  id = c('time'))
  plot_gbm <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable)) +
    ylim(4.2,4.7) +
    theme(legend.position = "none") + xlab("n") + ggtitle("Metode A med AV") + 
    theme(plot.title = element_text(hjust = 0.5))
  return(plot_gbm)
}
#Figue 3.5
MetodeB <- function(sigma,r,n_time,n_sim,T,S0,K,t){
  EhatV <- matrix(NA, n_time,n_sim) 
  for (k in 1:n_sim){
    print(k)
    z <- rnorm(n_time,0,1)
    S <- S0*exp((r-0.5*sigma^2)+sigma*z)
    V <- pmax(K-S,0)
    for (j in 1:n_time){
      EhatV[j,k] <- mean(V[1:j]) 
    }
  }
  Vhat <- EhatV*exp(-r)
  data_plot <- data.frame('time' =  seq(1, n_time, by = 1), "GBM" = Vhat)
  data_plot <- melt(data_plot,  id = c('time'))
  plot_gbm <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable)) +
    ylim(4.2,4.7) +
    theme(legend.position = "none")  + xlab("n") + ggtitle("Metode B")
  return(plot_gbm)
}
MetodeB_AV <- function(sigma,r,n_time,n_sim,T,S0,K,t){
  EhatVplus <- matrix(NA, n_time,n_sim)
  EhatVminus <- matrix(NA, n_time,n_sim)
  for (k in 1:n_sim){
    print(k)
    z <- rnorm(n_time,0,1)
    Splus <- S0*exp((r-0.5*sigma^2)+sigma*z)
    Sminus <- S0*exp((r-0.5*sigma^2)+sigma*(-z))
    Vplus <- pmax(K-Splus,0)
    Vminus <- pmax(K-Sminus,0)
    for (j in 1:n_time){
      EhatVplus[j,k] <- mean(Vplus[1:j]) 
      EhatVminus[j,k] <- mean(Vminus[1:j])
    }
  }
  Vhat <- (EhatVplus*exp(-r)+EhatVminus*exp(-r))/2
  data_plot <- data.frame('time' =  seq(1, n_time, by = 1), "GBM" = Vhat)
  data_plot <- melt(data_plot,  id = c('time'))
  plot_gbm <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable)) +
    ylim(4.2,4.7) +
    theme(legend.position = "none") + xlab("n") + ggtitle("Metode B med AV") +
    theme(plot.title = element_text(hjust = 0.5))
  return(plot_gbm)
}
grid.arrange(MetodeA(0.3,0.06,10000,10,1,5,10,0.02),
             MetodeA_AV(0.3,0.06,10000,10,1,5,10,0.02),
             MetodeB(0.3,0.06,10000,5,1,5,10,0.02),
             MetodeB_AV(0.3,0.06,10000,5,1,5,10,0.02),
             nrow = 2)
### Asian
Asiancall <- function(sigma,mu,n_time,n_sim,T,S0,K,q) {
  EhatV <- matrix(NA, n_time,n_sim)
  t <- T/n_time
  Splus <- matrix(NA,n_time + 1,n_sim)
  Sminus <- matrix(NA,n_time + 1,n_sim)
  Splus[1,] <- Sminus[1,] <- S0
  z <- matrix(rnorm(n_time*n_sim,0,1),n_time,n_sim)
  Payoff <- c()
  for (i in 2:(n_time+1)){
    Splus[i,] <- Splus[i-1,]*exp((mu-q-1/2*sigma^2)*t+sigma*sqrt(t)*z[i-1,])
    Sminus[i,] <- Sminus[i-1,]*exp((mu-q-1/2*sigma^2)*t+sigma*sqrt(t)*(-z[i-1,]))
  }
  for (i in 1:n_sim){
    Payoff[i] <- 1/2*(max(mean(Splus[,i])-K,0)+max(mean(Sminus[,i])-K,0))
  }
  Vhat <- exp(-mu*T)*mean(Payoff)
  return(Vhat)
}
Asianput <- function(sigma,mu,n_time,n_sim,T,S0,K) {
  t <- T/n_time
  Splus <- matrix(NA,n_time + 1,n_sim)
  Sminus <- matrix(NA,n_time + 1,n_sim)
  Splus[1,] <- Sminus[1,] <- S0
  z <- matrix(rnorm(n_time*n_sim,0,1),n_time,n_sim)
  Payoff <- c()
  for (i in 2:(n_time+1)){
    Splus[i,] <- Splus[i-1,]*exp((mu-1/2*sigma^2)*t+sigma*sqrt(t)*z[i-1,])
    Sminus[i,] <- Sminus[i-1,]*exp((mu-1/2*sigma^2)*t+sigma*sqrt(t)*(-z[i-1,]))
  }
  Payoff <- c()
  for (i in 1:n_sim){
    Payoff[i] <- exp(-mu*T)*1/2*(max(K-mean(Splus[,i]),0)+max(K-mean(Sminus[,i]),0))
  }
  EhatV <- mean(Payoff)
  Vhat <- EhatV
  return(Vhat)
}
blackcall <- function(sigma,r,T,S0,K,q){
  d1 <- 1/(sigma*sqrt(T)) * (log(S0/K)+(r-q+1/2*sigma^2)*(T))
  d2 <- d1 - sigma*sqrt(T)
  C <- exp(-r*T)*(exp((r-q)*T)*S0*pnorm(d1)-K*pnorm(d2))
  return(C)
}
blackput <- function(sigma,r,T,S0,K,q){
  P <- blackcall(sigma,r,T,S0,K,q) - S0 * exp(-q*T) + K* exp(-r*T)
  return(P)
}
EUvsA <- function(sigma,r,K) { return(c(Asiancall(sigma, r, 1000, 10000, 1, 10, K,0),
                                        blackcall(sigma,r,1,10,K,0),
                                        Asianput(sigma, r, 1000, 10000, 1, 10, K),
                                        blackput(sigma,r,1,10,K,0)))}
TabelAsian <- matrix(c(EUvsA(0.1,0.1,5),EUvsA(0.1,0.1,10),EUvsA(0.1,0.1,15),EUvsA(0.2,0.1,5),EUvsA(0.2,0.1,10),
               EUvsA(0.2,0.1,15),EUvsA(0.3,0.1,5),EUvsA(0.3,0.1,10),EUvsA(0.3,0.1,15)), nrow = 9,byrow = TRUE)
TabelAsian
