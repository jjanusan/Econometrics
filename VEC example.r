library(urca)  # data
library(vars)  
library(tsDyn) 


# define data 1958Q1 - 1984Q2
data(finland)
fin <- finland
nr_fin <- nrow(fin)

qy <- expand.grid(1:4, 1958:1984)[1:nr_fin,]
colnames(qy) <- c("q", "yyyy")
rownames(qy) <- NULL

### ADF tests
# Money
mon<-c(rep(0,20))
for (j in 1:length(mon)){
  aa<-fUnitRoots::adfTest(fin[,1],lags=j,type="c")
  mon[j]<-aa@test$statistic
}
# real income
rinc<-c(rep(0,20))
for (j in 1:length(rinc)){
  aa<-fUnitRoots::adfTest(fin[,2],lags=j,type="c")
  rinc[j]<-aa@test$statistic
}

# interest rate
ir<-c(rep(0,20))
for (j in 1:length(ir)){
  aa<-fUnitRoots::adfTest(fin[,3],lags=j,type="c")
  ir[j]<-aa@test$statistic
}
# inflation
infl<-c(rep(0,20))
for (j in 1:length(infl)){
  aa<-fUnitRoots::adfTest(fin[,4],lags=j,type="c")
  infl[j]<-aa@test$statistic
}


#### Testing if all variables are cointegrated
# quarterly centered dummy variables
qy$Q1 <- (qy$q==1)-1/4
qy$Q2 <- (qy$q==2)-1/4
qy$Q3 <- (qy$q==3)-1/4
dum_season <- qy[,-c(1,2)]

fin=ts(fin)

colnames(fin)=c("money","real income","interest rate","infl")

  plots <- list()
for (i in 1:4){
  plots[[i]] <- autoplot(fin[,i], color = "blue")+labs(x = "Time", y = colnames(fin)[i])
}
gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2, nrow = 2)

# differenced data
dif <- as.data.frame(diff(as.matrix(fin), lag = 1))

# cointegration Test
coint_ca.jo <- ca.jo(
  fin, ecdet = "none", type  = "eigen", K = 2, 
  spec = "transitory", season = 4, dumvar = NULL)
summary(coint_ca.jo)

# Trace test
coint_ca.jo <- ca.jo(
  fin, ecdet = "none", type  = "trace", K = 2, 
  spec = "transitory", season = 4, dumvar = NULL)
summary(coint_ca.jo)

# VECM 

cajorls_ca.jo <- cajorls(coint_ca.jo, r=2)
vec2var_ca.jo <- vec2var(coint_ca.jo, r=2)

