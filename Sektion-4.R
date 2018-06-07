###############OPSÆRNING#################
library(grid) #c
library(gridExtra) #c
library(gridGraphics) #c
library(reshape) #c
library(ggplot2) #c
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_171') #Java lokation til xlsx pakken
                                                               #er måske ikke nødvendig
library(xlsx) #c
library(pracma) #c
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
###############OPSÆRNING#################

## Black scholes sektor 4
#### Figur 8
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
blackcallplotr <- Vectorize(function(x) {blackcall(sigma = 0.1, x, T = 1, S0 = 10, K = 5, q=0)})
blackcallplotT <- Vectorize(function(x) {blackcall(sigma = 0.1, r = 0.06, x, S0 = 10, K = 5, q=0)})
blackcallplotS0 <- Vectorize(function(x) {blackcall(sigma = 0.1, r = 0.06, T=1, x, K = 5, q=0)})

grid.arrange(
  ggplot(data.frame(x=c(0,1)),aes(x)) + stat_function(fun = blackcallplotr) + xlab("r") + ylab("c_t"),
  ggplot(data.frame(x=c(1,10)),aes(x)) + stat_function(fun = blackcallplotT) + xlab("T") + ylab("c_t"),
  ggplot(data.frame(x=c(1,10)),aes(x)) + stat_function(fun = blackcallplotS0) + xlab("S0") + ylab("c_t")
  , ncol = 3)

##sigma og mu 500
# setwd("nuværende placerings af R fil")
SP500EU <- read.xlsx("./Data/SP500EU.xlsx", sheetIndex = 1)
SEU <- SP500EU$GSPC.Close

a <- function(n,S) {return ((log(S[n])-log(S[1]))/n)}
b <- function(n,S){
  b <- 0
  for (i in 2:n){
    b <- b + ((log(S[i]/S[i-1]) - a(n,S))^2)/n
  }
  return (b)
}
t <- 38/length(SEU) # 
n <- length(SEU)
sigmahat500 <- sqrt(b(n,SEU)/t); sigmahat500
muhat500 <- 1/2*sigmahat500^2+a(n,SEU)/(t); muhat500

##Figur 10
OptionEUcall1805 <-read.xlsx("./Data/CallsEU.xlsx", sheetIndex = 9)
OptionEUput1805  <-read.xlsx("./Data/PutsEU.xlsx", sheetIndex = 9)

blackputpris1805 <- Vectorize(function(K) { return(blackput(sigmahat500, 
                                                            0.0224, T = 15/252, SEU[length(SEU)],K, q = 0))})
blackcallpris1805 <- Vectorize(function(K) { return(blackcall(sigmahat500, 
                                                              0.0224, T = 15/252, SEU[length(SEU)],K, q = 0))})
xeucall <- blackcallpris1805(OptionEUcall1805$Strike)
yeucall <- OptionEUcall1805$Last
xeuput <- blackputpris1805(OptionEUput1805$Strike)
yeuput <- OptionEUput1805$Last

grid.arrange(qplot(xeucall,yeucall) + geom_abline(intercept = 0, slope = 1) + xlab("Teoretiske Call") + ylab("Observeret Call"),
             qplot(xeuput,yeuput) + geom_abline(intercept = 0, slope = 1)+ xlab("Teoretiske Put") + ylab("Observeret Put"),
             nrow = 1)

##AM sigma og mu
SP500AM <- read.xlsx("./Data/SP500AM.xlsx", sheetIndex = 1)
SAM <- SP500AM$SPY.Close

t <- 11/length(SAM) 
n <- length(SAM)
sigmahat500AM <- sqrt(b(n,SAM)/t); sigmahat500AM
muhat500AM <- 1/2*sigmahat500^2+a(n,SAM)/(t); muhat500AM

## Figur 11
OptionAMcall1805 <-read.xlsx("./Data/CallsAM.xlsx", sheetIndex = 9)
OptionAMput1805  <-read.xlsx("./Data/PutsAM.xlsx", sheetIndex = 9)

OptionAMcall3112 <-read.xlsx("./Data/CallsAM.xlsx", sheetIndex = 25)
OptionAMput3112  <-read.xlsx("./Data/PutsAM.xlsx", sheetIndex = 25)

blackputprisAM1805 <- Vectorize(function(K) { return(blackput(sigmahat500AM, 
                                                              0.0224, T = 15/252, SAM[length(SAM)],K, q= 0))})
blackcallprisAM1805 <- Vectorize(function(K) { return(blackcall(sigmahat500AM, 
                                                                0.0224, T = 15/252, SAM[length(SAM)],K, q= 0))})
blackputprisAM3112 <- Vectorize(function(K) { return(blackput(sigmahat500, 
                                                              0.0224, T = 247/365, SAM[length(SAM)],K, q = 0))})
blackcallprisAM3112 <- Vectorize(function(K) { return(blackcall(sigmahat500, 
                                                                0.0224, T = 247/365, SAM[length(SAM)],K, q = 0))})
xamcall1 <- blackcallprisAM1805(OptionAMcall1805$Strike)
yamcall1 <- OptionAMcall1805$Last
xamput1 <- blackputprisAM1805(OptionAMput1805$Strike)
yamput1 <- OptionAMput1805$Last

xamcall2 <- blackcallprisAM3112(OptionAMcall3112$Strike)
yamcall2 <- OptionAMcall3112$Last
xamput2 <- blackputprisAM3112(OptionAMput3112$Strike)
yamput2 <- OptionAMput3112$Last

grid.arrange(arrangeGrob(qplot(xamcall1,yamcall1) + geom_abline(intercept = 0, slope = 1)+ xlab("Teoretiske EU Call") + ylab("Observeret AM Call"),
                         qplot(xamput1,yamput1) + geom_abline(intercept = 0, slope = 1) + xlab("Teoretiske EU Put") + ylab("Observeret AM Put"), top = "18-05-2018", nrow = 1),
             arrangeGrob(qplot(xamcall2,yamcall2) + geom_abline(intercept = 0, slope = 1) + xlab("Teoretiske EU Call") + ylab("Observeret AM Call"),
                         qplot(xamput2,yamput2) + geom_abline(intercept = 0, slope = 1)+ xlab("Teoretiske EU Put") + ylab("Observeret AM Put"), top = "31-12-2018", nrow = 1),
             nrow = 2)  

## Implicit vol
IMV <- Vectorize(function(K,l){bisect(function(sigma){
  blackcall(sigma = sigma, r = muhat500, T = 17/252, q = 0, K = K,S0 =  SEU[length(SEU)])- l
} , -1000, 100000)$root})
IMV1805 <- IMV(OptionEUcall1805$Strike,OptionEUcall1805$Last)
IMV1805P <- IMV1805[IMV1805 > 0]  #Fjerne de negaive vol
SK <- OptionEUcall1805$Strike[which(IMV1805 > 0)]/SEU[length(SEU)]
grid.arrange(
  qplot(SK,IMV1805P) + geom_line(linetype="dashed", color="blue")+
    geom_point(color="red", size=3) + xlab("S/K") + ylab("Implied Volatility"),
  qplot(SK, OptionEUcall1805$OI[which(IMV1805 > 0)]/SEU[length(SEU)],geom="linerange",ymin = 0, ymax = OptionEUcall1805$OI[which(IMV1805 > 0)]/SEU[length(SEU)]) 
  + xlab("S/K") + ylab("Open Interest"), ncol =2)

