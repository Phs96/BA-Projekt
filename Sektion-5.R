###############OPSÆRNING#################
library(grid)  #c
library(gridExtra) #c
library(gridGraphics) #c
library(reshape) #c
library(ggplot2) #c
library(ggpubr) #c
library(pracma) #c
###############OPSÆRNING#################
##########################
#A_t simulationen er uoptimeret
##########################
##### Sektion 6
r0 = 0.01
a = 0.15
b = 0.042
lambda = -0.23
sigmar = 0.01
rho = -0.15
sigmahat = 0.1827096
muhat <- 0.1021119
q = 0.025

A_udvikling <- function(n_time, T, x, n_sim){
  A_t <- matrix(NA, n_time + 1, n_sim)
  A_t[1,] <- 1000
  q <- 0.025
  dt <- T/n_time
  B <- function(t,T) {return(1/a*(1-exp(-a*(T-t))))}
  A <- function(t,T) {return((sigmar^2/(2*a^2)-b + lambda*sigmar/a)*((T-t)-B(t,T))-sigmar^2/(4*a)*B(t,T)^2)}
  for (j in 1:n_sim){
    r_t <- c(r0)
    t <- 0
    w1 <- rnorm(n_time,0,1)
    w2 <- rnorm(n_time,0,1)
    w3 <- rho * sqrt(dt) * w1 + sqrt(1-rho^2) * sqrt(dt) * w2
    for (i in 2:(n_time+1)){
      A_t[i,j] <-  A_t[i-1,j] + x[1]*A_t[i-1,j]*r_t[i-1]*dt + 
        x[2]*A_t[i-1,j]*(muhat*dt + sigmahat*w3[i-1]) +
        x[3]*A_t[i-1,j]*sum(1/10*((r_t[i-1]-lambda*sigmar*
                                     B(t,floor(t)+seq(1:10)))*dt-
                                    sigmar*B(t,floor(t)+seq(1:10))*sqrt(dt)* w1[i-1]))
      t <- t + dt
      r_t[i]  <- r_t[i-1] + a*(b-r_t[i-1])*dt + sigmar * sqrt(dt) * w1[i-1] 
    }
  }
  q <- matrix(NA, nrow = 121, ncol = 4)
  gennemsnit <- c()
  for (i in 1:121){
    q[i,1:3] <- quantile(A_t[i,], probs = c(0.005,0.5,0.995))
    q[i,4] <- mean(A_t[i,])
  }
  data_plot <- data.frame('time' =  seq(0, 10, by = 1/12), "GBM" = q)
  data_plot <- melt(data_plot,  id = c('time'))
  levels(data_plot[,2]) <- c("0.5% quantile", "50% quantile", "99.5% quantile", "Gennemsnit")
  colnames(data_plot)[which(names(data_plot) == "variable")] <- "Linjer"
  data_plot2 <- data_plot[,!(names(data_plot) %in% c("Linjer"))][1:121,]
  data_plot2[,3] <- data_plot[243:363,3]
  data_plot2[,4] <- gennemsnit
  
  plot_gbm <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = Linjer,group = Linjer),size = 1.2) +
    geom_ribbon(data=data_plot2, 
                aes(ymin=value,ymax=V3), fill="blue", alpha="0.2") +
    theme(legend.title=element_blank())
  return(plot_gbm)
}
q1 <- A_udvikling(120,10,c(1,0,0),10000) + ggtitle("X = (1,0,0)")
q2 <- A_udvikling(120,10,c(0,1,0),10000) + ggtitle("X = (0,1,0)")
q3 <- A_udvikling(120,10,c(0,0,1),10000) + ggtitle("X = (0,0,1)")
q4 <- A_udvikling(120,10,c(1/3,1/3,1/3),10000) + ggtitle("X = (1/3,1/3,1/3)")
ggarrange(q1,q2,q3,q4,common.legend = TRUE, legend="bottom" ,nrow = 2, ncol = 2)

## Sandsynlighed
Sandsynlighed <- function(n_time, T, x, n_sim){
  A_t <- matrix(NA, n_time + 1, n_sim)
  A_t[1,] <- 1000
  dt <- T/n_time
  t <- 0
  B <- function(t,T) return(1/a*(1-exp(-a*(T-t))))
  A <- function(t,T) {return((sigmar^2/(2*a^2)-b + lambda*sigmar/a)*((T-t)-B(t,T))-sigmar^2/(4*a)*B(t,T)^2)}
  for (j in 1:n_sim){
    r_t <- r0
    w1 <- rnorm(n_time,0,1)
    w2 <- rnorm(n_time,0,1)
    w3 <- rho * sqrt(dt) * w1 + sqrt(1-rho^2) * sqrt(dt) * w2
    for (i in 2:(n_time+1)){
      A_t[i,j] <-  A_t[i-1,j] + x[1]*A_t[i-1,j]*r_t*dt + 
        x[2]*A_t[i-1,j]*(muhat*dt + sigmahat*w3[i-1]) +
        x[3]*A_t[i-1,j]*sum(1/10*((r_t-lambda*sigmar*B(t,floor(t)+seq(1:10)))*dt-
                                    sigmar*B(t,floor(t)+seq(1:10))*sqrt(dt)* w1[i-1]))
      r_t  <- r_t + a*(b-r_t)*dt + sigmar * sqrt(dt) * w1[i-1] 
      t <- t + dt
    }
  }
  LT <- 1000*(1+q)^T
  BT <- max(A_t[121,] - LT)
  return(sum(A_t[121,] < 1000*(1+q)^T)/n_sim)
}
#Det tager 4 timer at køre
Xs <- rev(seq(0, 1, 0.01))
Xb <- seq(0, 1, 0.01)
p <- matrix(NA, ncol = length(Xb), nrow = length(Xs))
for (i in 1:length(Xs)){
  for (j in 1:i){
    print(c(i,j))
    Xbeta <- rev(Xb[1:i])
    p[length(Xs)-i+1,j] <-   Sandsynlighed(120,10,c(Xbeta[j],Xs[i],Xb[j]), 2000)
  }
}
col.pal<-colorRampPalette(c("dark green","green","yellow", "red","dark red"))
colors<-col.pal(100)
z.facet.center <- (p[-1, -1] + p[-1, -ncol(p)] + p[-nrow(p), -1] + p[-nrow(p), -ncol(p)])/4
z.facet.range<-cut(z.facet.center, 100)
persp(rev(Xs),Xb,p, theta = 120, col=colors[z.facet.range],expand = 0.5,ticktype = "detailed",
      phi = 13, xlab= "Xb",ylab = "Xs",border = NA)


### Quantiler
Quantiles <- function(n_time, T, x, n_sim){
  A_t <- matrix(NA, n_time + 1, n_sim)
  A_t[1,] <- 1000
  q <- 0.025
  dt <- T/n_time
  Bt <- matrix(NA, ncol = n_sim, nrow = (n_time+1))
  loss <- matrix(NA, ncol = n_sim, nrow = (n_time-11))
  B <- function(t,T) {return(1/a*(1-exp(-a*(T-t))))}
  A <- function(t,T) {return((sigmar^2/(2*a^2)-b + lambda*sigmar/a)*((T-t)-B(t,T))-sigmar^2/(4*a)*B(t,T)^2)}
  Bt[1,] <- A_t[1,]-1000*(1+q)^(T)*exp(A(0,T)-B(0,T)*r0)
  for (j in 1:n_sim){
    r_t <- c(r0)
    t <- 0
    w1 <- rnorm(n_time,0,1)
    w2 <- rnorm(n_time,0,1)
    w3 <- rho * sqrt(dt) * w1 + sqrt(1-rho^2) * sqrt(dt) * w2
    for (i in 2:(n_time+1)){
      A_t[i,j] <-  A_t[i-1,j] + x[1]*A_t[i-1,j]*r_t[i-1]*dt + 
        x[2]*A_t[i-1,j]*(muhat*dt + sigmahat*w3[i-1]) +
        x[3]*A_t[i-1,j]*sum(1/10*((r_t[i-1]-lambda*sigmar*B(t,floor(t)+seq(1:10)))*dt-
                                    sigmar*B(t,floor(t)+seq(1:10))*sqrt(dt)* w1[i-1]))
      t <- t + dt
      r_t[i]  <- r_t[i-1] + a*(b-r_t[i-1])*dt + sigmar * sqrt(dt) * w1[i-1] 
      Bt[i,j] <- max(A_t[i,j] - 1000*(1+q)^(T)*exp(A(t,T)-B(t,T)*r_t[i]),0)
    }
    t <- 0
    for(i in 1:(n_time-11)){
      loss[i,j] <- Bt[i,j] - Bt[i + 12,j]*exp(A(t,dt*12+t)-B(t,dt*12+t)*r_t[i])
      t <- t + dt
    }
  }
  qb <- matrix(NA, ncol=3, nrow= (n_time+1))    #Quantile for B()
  for(i in 1:(n_time+1)){
    qb[i,] <- quantile(Bt[i,], probs = c(0.005,0.5,0.995))
  }
  ql <- c()                                     #Quantile for loss
  for(i in 1:(n_time-11)){
    ql[i] <- quantile(loss[i,],probs = c(0.995))
  }
  Cr <- matrix(NA, ncol = n_sim, nrow = (n_time-11))  #Cr
  qCr <- matrix(NA, ncol = 3, nrow = (n_time-11))     #Quantile for Cr
  cRp <- c()                                          #Punktvise sandsynlighed for insolvens
  for(i in 1:(n_time-11)){
    Cr[i,] <- Bt[i,]/ql[i]
    qCr[i,] <- quantile(Cr[i,], probs = c(0.005,0.5,0.995)) 
  }
  for(i in 1:(n_time - 11)){
    længde <- length(Cr[i,])
    cRp[i] <- sum(Cr[i,] < 1)/længde
    if (length(which(Cr[i,] < 1) > 0)){
      Cr <- Cr[,-which(Cr[i,] < 1)]
    }
  }
  cp <- c()                                           #Akkumuleret sandsynlighed
  for(i in 1:(n_time-11)){
    cp[i] = prod(1-cRp[1:(i-1)])*cRp[i]+(1-prod(1-cRp[1:(i-1)]))  
  }
  
  data_plot <- data.frame('time' =  seq(0, T-1, by = 1/12), "GBM" = qCr)
  data_plot <- melt(data_plot,  id = c('time'))
  levels(data_plot[,2]) <- c("0.5% quantile", "50% quantile", "99.5% quantile")
  qCrplot <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable,group = variable)) + ggtitle("Quantile af Cr") +theme(legend.position = c(0.2, 0.9)) +theme(legend.title=element_blank())
  
  sshplot <- qplot(seq(0, T-1, by = 1/12),cp, main = "Insolvent Sandsynlighed") 
  
  
  data_plot <- data.frame('time' =  seq(0, T, by = 1/12), "GBM" = qb)
  data_plot <- melt(data_plot,  id = c('time'))
  levels(data_plot[,2]) <- c("0.5% quantile", "50% quantile", "99.5% quantile")
  Bplot <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable,group = variable)) + ggtitle("B quantiles") + xlim(c(0,9)) + theme(legend.position="bottom") +theme(legend.title=element_blank())
  
  data_plot <- data.frame('time' =  seq(0, T-1, by = 1/12), "GBM" = ql)
  data_plot <- melt(data_plot,  id = c('time'))
  levels(data_plot[,2]) <- c("Loss")
  qLossplot <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable,group = variable)) + ggtitle("Loss") + theme(legend.position="bottom") +theme(legend.title=element_blank())
  
  return(grid.arrange(Bplot, qLossplot, sshplot , qCrplot,ncol = 2))
}
Quantiles(120,10,c(0,0.05,0.95), 10000)

### Eta og lt
Lt_eta <- function(n_time, T, x, n_sim){
  A_t <- matrix(NA, n_time + 1, n_sim)
  A_t[1,] <- 1000
  dt <- T/n_time
  t <- 0
  B <- function(t,T) return(1/a*(1-exp(-a*(T-t))))
  A <- function(t,T) {return((sigmar^2/(2*a^2)-b + lambda*sigmar/a)*((T-t)-B(t,T))-sigmar^2/(4*a)*B(t,T)^2)}
  for (j in 1:n_sim){
    r_t <- r0
    w1 <- rnorm(n_time,0,1)
    w2 <- rnorm(n_time,0,1)
    w3 <- rho * sqrt(dt) * w1 + sqrt(1-rho^2) * sqrt(dt) * w2
    for (i in 2:(n_time+1)){
      A_t[i,j] <-  A_t[i-1,j] + x[1]*A_t[i-1,j]*r_t*dt + 
        x[2]*A_t[i-1,j]*(r_t*dt + sigmahat*w3[i-1]) +
        A_t[i-1,j]*x[3]*sum(1/10*(r_t*dt-
                                    sigmar*B(t,floor(t)+seq(1:10))*sqrt(dt)* w1[i-1]))
      r_t  <- r_t + a*(b-lambda*sigmar/a-r_t)*dt + sigmar * sqrt(dt) * w1[i-1]
      t <- t + dt
    }
  }
  LT <- 1000*(1+q)^T
  bonus_dis <- c()
  bonus <- c()
  for (i in 1:length(seq(13,121,12))){
    bonus_dis[i] <- mean(exp(A(0,T)-B(0,T)*r0)*(pmax(A_t[seq(13,121,12)[i],]-1000*(1+q)^i,0))) #Prisen af call optioner - Garenteret rente
    bonus[i] <-   mean(pmax(A_t[seq(13,121,12)[i],]-1000*(1+q)^i,0)) 
  }
  bonus_dis[11] <- sum(bonus_dis[1:10]) 
  bonus[11] <- sum(bonus[1:10])
  BT <- mean(pmax(A_t[121,] - LT,0)) 
  
  eta1 <- bisect(function(eta){ return(exp(A(0,T)-B(0,T)*r0)*(eta*BT+LT)-1000)}, 0,1)
  eta2 <- bisect(function(eta){ return(exp(A(0,T)-B(0,T)*r0)*(eta*bonus[11]+LT)-1000)}, 0,1)
  eta <- c(eta1$root,eta2$root)
  return(eta)
}
Lt_eta(120, 10, c(0,0.05,0.95), 10000)


## Nulkupon hedge
CRZCBhedge <- function(n_time, T, x, n_sim){
  A_ti <- matrix(NA, n_time + 1, n_sim)
  A_t <- matrix(NA, n_time + 1, n_sim)
  q <- 0.025
  dt <- T/n_time
  Bt <- matrix(NA, ncol = n_sim, nrow = (n_time+1))
  qb <- matrix(NA, ncol = 3, nrow = (n_time+1))
  loss <- matrix(NA, ncol = n_sim, nrow = (n_time-11))
  B <- function(t,T) {return(1/a*(1-exp(-a*(T-t))))}
  A <- function(t,T) {return((sigmar^2/(2*a^2)-b + lambda*sigmar/a)*((T-t)-B(t,T))-sigmar^2/(4*a)*B(t,T)^2)}
  A_ti[1,] <- 1000 - 1000*(1+q)^(T)*exp( A(0,T) - B(0,T)*r0) #Det der bliver investeret uden ZCB
  A_t[1,] <- 1000   #Værdien af alle aktiver med ZCB
  Bt[1,] <- A_t[1,]-1000*(1+q)^(T)*exp(A(0,T)-B(0,T)*r0)
  for (j in 1:n_sim){
    r_t <- c(r0)
    t <- 0
    w1 <- rnorm(n_time,0,1)
    w2 <- rnorm(n_time,0,1)
    w3 <- rho * sqrt(dt) * w1 + sqrt(1-rho^2) * sqrt(dt) * w2
    for (i in 2:(n_time+1)){
      A_ti[i,j] <-  A_ti[i-1,j] + x[1]*A_ti[i-1,j]*r_t[i-1]*dt + 
        x[2]*A_ti[i-1,j]*(muhat*dt + sigmahat*w3[i-1]) +
        A_ti[i-1,j]*x[3]*sum(1/10*((r_t[i-1]-lambda*sigmar*B(t,floor(t)+seq(1:10)))*dt-
                                     sigmar*B(t,floor(t)+seq(1:10))*sqrt(dt)* w1[i-1]))
      t <- t + dt
      r_t[i]  <- r_t[i-1] + a*(b-r_t[i-1])*dt + sigmar * sqrt(dt) * w1[i-1] 
      A_t[i,j] <- max(A_ti[i,j] + 1000*(1+q)^(T)*exp(A(t,T)-B(t,T)*r_t[i]),0)
      Bt[i,j] <- max(A_t[i,j] - 1000*(1+q)^(T)*exp(A(t,T)-B(t,T)*r_t[i]),0)
    }
    t <- 0
    for(i in 1:(n_time-11)){
      loss[i,j] <- Bt[i,j] - Bt[i + 12,j]*exp(A(t,dt*12+t)-B(t,dt*12+t)*r_t[i])
      t <- t + dt
    }
  }
  
  qb <- matrix(NA, ncol=3, nrow= (n_time+1))
  for(i in 1:(n_time+1)){
    qb[i,] <- quantile(Bt[i,], probs = c(0.005,0.5,0.995))
  }
  ql <- matrix(NA, ncol=1, nrow= (n_time-11))
  for(i in 1:(n_time-11)){
    ql[i] <- quantile(loss[i,],probs = c(0.995))
  }
  Cr <- matrix(NA, ncol = n_sim, nrow = (n_time-11))
  qCr <- matrix(NA, ncol = 3, nrow = (n_time-11))
  cRp <- c()
  for(i in 1:(n_time-11)){
    Cr[i,] <- Bt[i,]/ql[i]
    qCr[i,] <- quantile(Cr[i,], probs = c(0.005,0.5,0.995)) 
  }
  for(i in 1:(n_time - 11)){
    længde <- length(Cr[i,])
    cRp[i] <- sum(Cr[i,] < 1)/længde
    if (length(which(Cr[i,] < 1) > 0)){
      Cr <- Cr[,-which(Cr[i,] < 1)]
    }
  }
  cp <- c()
  for(i in 1:(n_time-11)){
    cp[i] = prod(1-cRp[1:(i-1)])*cRp[i]+(1-prod(1-cRp[1:(i-1)]))
  }
  
  data_plot <- data.frame('time' =  seq(0, T-1, by = 1/12), "GBM" = qCr)
  data_plot <- melt(data_plot,  id = c('time'))
  levels(data_plot[,2]) <- c("0.5% quantile", "50% quantile", "99.5% quantile")
  qCrplot <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable,group = variable)) + ggtitle("Quantile af Cr") + theme(legend.title=element_blank()) +  theme(legend.position="bottom") 
  
  sshplot <- qplot(seq(0, T-1, by = 1/12),cp, main = "Insolvent Sandsynlighed") 
  
  data_plot <- data.frame('time' =  seq(0, T, by = 1/12), "GBM" = qb)
  data_plot <- melt(data_plot,  id = c('time'))
  levels(data_plot[,2]) <- c("0.5% quantile", "50% quantile", "99.5% quantile")
  Bplot <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable,group = variable)) + ggtitle("B quantiles") + xlim(c(0,9)) + theme(legend.title=element_blank()) +  theme(legend.position="bottom") 
  
  data_plot <- data.frame('time' =  seq(0, T-1, by = 1/12), "GBM" = ql)
  data_plot <- melt(data_plot,  id = c('time'))
  levels(data_plot[,2]) <- c("99.5% quantile")
  qLossplot <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable,group = variable)) + ggtitle("Loss quatile") + theme(legend.title=element_blank()) +  theme(legend.position="bottom") 
  return(grid.arrange(Bplot,qLossplot,qCrplot,sshplot,ncol = 2,nrow = 2)) 
  return(ggarrange(qLossplot,qCrplot, legend="bottom",ncol = 2,nrow = 1)) 
}
CRZCBhedge(120,10,c(0,1,0), 10000)

###Swap Hedge
CRSWAPhedge <- function(n_time, T, x, n_sim){
  A_ti <- matrix(NA, n_time + 1, n_sim)
  A_t <- matrix(NA, n_time + 1, n_sim)
  q <- 0.025
  dt <- T/n_time
  Bt <- matrix(NA, ncol = n_sim, nrow = (n_time+1))
  rt <- matrix(NA, ncol = n_sim, nrow = (n_time+1))
  qb <- matrix(NA, ncol = 3, nrow = (n_time+1))
  B <- function(t,T) {return(1/a*(1-exp(-a*(T-t))))}
  loss <- matrix(NA, ncol = n_sim, nrow = (n_time-11))
  A <- function(t,T) {return((sigmar^2/(2*a^2)-b + lambda*sigmar/a)*((T-t)-B(t,T))-sigmar^2/(4*a)*B(t,T)^2)}
  A_ti[1,] <- 1000 
  A_t[1,] <- A_ti[1,]+1000*(1+q)^(T)*exp(A(0,T)-B(0,T)*r0)
  Bt[1,] <- A_t[1,]-2*1000*(1+q)^(T)*exp(A(0,T)-B(0,T)*r0) 
  k <- function(t,r){ return(1000*(1+q)^(T)*exp(A(0,T)-B(0,T)*r0) / (exp(A(t,t+dt)-B(t,t+dt)*r)) - 1000*(1+q)^(T)*exp(A(0,T)-B(0,T)*r0)) }
  for (j in 1:n_sim){
    r_t <- c(r0)
    t <- 0
    w1 <- rnorm(n_time,0,1)
    w2 <- rnorm(n_time,0,1)
    w3 <- rho * sqrt(dt) * w1 + sqrt(1-rho^2) * sqrt(dt) * w2
    for (i in 2:(n_time+1)){
      A_ti[i,j] <-  A_ti[i-1,j] + x[1]*A_ti[i-1,j]*r_t[i-1]*dt + 
        x[2]*A_ti[i-1,j]*(muhat*dt + sigmahat*w3[i-1]) +
        A_ti[i-1,j]*x[3]*sum(1/10*((r_t[i-1]-lambda*sigmar*B(t,floor(t)+seq(1:10)))*dt-
                                     sigmar*B(floor(t),t+seq(1:10))*sqrt(dt)* w1[i-1]))
      t <- t + 1/12
      A_ti[i,j] <- A_ti[i,j] - k(t,r_t[i-1])    #k er floting leg
      r_t[i]  <- r_t[i-1] + a*(b-r_t[i-1])*dt + sigmar * sqrt(dt) * w1[i-1] 
      A_t[i,j] <- A_ti[i,j] + 1000*(1+q)^(T)*exp(A(t,T)-B(t,T)*r_t[i])
      Bt[i,j] <- max(A_t[i,j] - 1000*(1+q)^(T)*exp(A(t,T)-B(t,T)*r_t[i]) - 1000*(1+q)^(T)*exp(A(0,T)-B(0,T)*r_t[1]) ,0)
    }
    rt[,j] <- r_t
    t <- 0
    for(i in 1:(n_time-11)){
      loss[i,j] <- Bt[i,j] - Bt[i + 12,j]*exp(A(t,dt*12+t)-B(t,dt*12+t)*r_t[i])
      t <- t + dt
    }
  }
  qb <- matrix(NA, ncol=3, nrow= (n_time+1))
  for(i in 1:(n_time+1)){
    qb[i,] <- quantile(Bt[i,], probs = c(0.005,0.5,0.995))
  }
  
  ql <- matrix(NA, ncol=1, nrow= (n_time-11))
  for(i in 1:(n_time-11)){
    ql[i] <- quantile(loss[i,],probs = c(0.995))
  }
  
  Cr <- matrix(NA, ncol = n_sim, nrow = (n_time-11))
  qCr <- matrix(NA, ncol = 3, nrow = (n_time-11))
  cRp <- c()
  for(i in 1:(n_time-11)){
    Cr[i,] <- Bt[i,]/ql[i]
    qCr[i,] <- quantile(Cr[i,], probs = c(0.005,0.5,0.995)) 
  }
  
  for(i in 1:(n_time - 11)){
    længde <- length(Cr[i,])
    cRp[i] <- sum(Cr[i,] < 1)/længde
    if (length(which(Cr[i,] < 1) > 0)){
      Cr <- Cr[,-which(Cr[i,] < 1)]
    }
  }
  cp <- c()
  for(i in 1:(n_time-11)){
    cp[i] = prod(1-cRp[1:(i-1)])*cRp[i]+(1-prod(1-cRp[1:(i-1)]))
  }
  
  data_plot <- data.frame('time' =  seq(0, T-1, by = 1/12), "GBM" = qCr)
  data_plot <- melt(data_plot,  id = c('time'))
  levels(data_plot[,2]) <- c("0.5% quantile", "50% quantile", "99.5% quantile")
  qCrplot <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable,group = variable)) + ggtitle("Quantile af Cr") + xlim(c(0,9))+ theme(legend.position = c(0.14, 0.85)) +theme(legend.title=element_blank())
  
  sshplot <- qplot(seq(0, T-1, by = 1/12),cp, main = "Insolvent Sandsynlighed") 
  
  data_plot <- data.frame('time' =  seq(0, T, by = 1/12), "GBM" = qb)
  data_plot <- melt(data_plot,  id = c('time'))
  levels(data_plot[,2]) <- c("0.5% quantile", "50% quantile", "99.5% quantile")
  Bplot <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable,group = variable)) + ggtitle("B quantiles") + xlim(c(0,9))+ theme(legend.position = c(0.14, 0.85)) +theme(legend.title=element_blank())
  
  data_plot <- data.frame('time' =  seq(0, T-1, by = 1/12), "GBM" = ql)
  data_plot <- melt(data_plot,  id = c('time'))
  levels(data_plot[,2]) <- c("99.5% quantile")
  qLossplot <- ggplot(data_plot, aes(time, value)) +
    geom_line(aes(colour = variable,group = variable)) + ggtitle("Loss quatile") + xlim(c(0,9))+ theme(legend.position = c(0.14, 0.95)) +theme(legend.title=element_blank())
  
  return(grid.arrange(Bplot,qLossplot,qCrplot,sshplot + xlab("time"), nrow = 2))
  return(grid.arrange(Bplot, qLossplot, sshplot, qCrplot,ncol = 2))
}
CRSWAPhedge(120,10,c(0,0.05,0.95), 10000)
CRSWAPhedge(120,10,c(1,0,0), 10000)
