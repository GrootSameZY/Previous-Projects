setwd("D:/rcode/qf5210/GroupPro/")
library(fGarch)
library(rugarch)
library(tseries)
library(TSA)
library(ggplot2)
library(vars)
library(zoo)
library(reshape2)

adj = read.csv("adjp.csv",header = T,row.names = 1) #Adjusted Prices
total_assets = read.csv("datelineta.csv",header = T,row.names = 1) #Daily Total Assets
fintech = read.csv("Fintechindex_adjusted.csv",header = T,row.names = 1) #Fintech Index
MIX = read.csv("adjindexprice.csv",header = T,row.names = 1) #Mixed Index from SS,SZ,HK Indices
log_MIX = log(MIX[2:nrow(MIX),]/MIX[1:(nrow(MIX)-1),]) #Daily Log Returns of Mixed Index
MIX1520 = log_MIX[which(rownames(as.data.frame(log_MIX))=="2015/1/6"):which(rownames(as.data.frame(log_MIX))=="2020/6/30"),4] #Mixed Index: 20150106--20200630

#Compute Bank Returns(Industry)
R = log(adj[2:nrow(adj),]/adj[1:(nrow(adj)-1),]) #Daily Log Returns
R[is.na(R)] <- 0
row.sums = apply(total_assets,1,sum) #Sum of Total Weights
row.sums = row.sums[-1]
w = total_assets[-1,]/row.sums #Weights for Each Stock
R_bk = apply(R*w,1,sum)#Returns of Bank Industry: Weighted Sum
#Split the Bank Returns into two parts: 2010--2015,2015--2020
R_bk_1015 = R_bk[1:which(rownames(as.data.frame(R_bk))=="2015/9/2")] #20100902---20150902
R_bk_1520 = R_bk[which(rownames(as.data.frame(R_bk))=="2015/1/6"):which(rownames(as.data.frame(R_bk))=="2020/6/30")] #20150106 --- 20200630 

#Fintech Index
Fin_Ind = fintech[2:which(rownames(as.data.frame(fintech))=="2020/6/30"),which(colnames(as.data.frame(fintech))=="Fintech.Index")] #Composite Fintech Index: 20150106--20200630
Fin_sep_log = as.data.frame(fintech[2:which(rownames(as.data.frame(fintech))=="2020/6/30"),c(2,4,6,8,10)]) #Separate Fintech Factors(log)
Fin_sep = as.data.frame(fintech[2:which(rownames(as.data.frame(fintech))=="2020/6/30"),c(1,3,5,7,9,11)]) #Separate Fintech Factors (fintech mktcap adj.)
SC_Fin_sep = scale(Fin_sep,center = TRUE,scale = TRUE)  #Normalized Scaling of Fintech Factors
Cor1 = cor(cbind(R_bk_1520,Fin_sep_log))
Cor2 = cor(cbind(R_bk_1520,SC_Fin_sep))

############################################################################################################
# Tests for Properties of FTS
basicStats(R_bk) #Skewness:0.026234(Right Skewed), Excess Kurtosis:2.940971
#Time Series Plots(Non-diff, Diff)
par(mfrow = c(2,1))
p <- plot.ts(R_bk,xlab ='Time',ylab ='Bank Returns')
p2 <- plot.ts(diff(R_bk),xlab ='Time',ylab ='Diff of Bank Returns')
#Stationarity: Stationary
adf.test(R_bk,k = 1) #lag:13 (Default)
summary(ur.df(R_bk,type = "trend",selectlags = "BIC")) #tau3<0: No Unit Roots 
#Autocorrelations
par(mfrow = c(2,1))
acf(R_bk)
pacf(R_bk)
par(mfrow = c(2,1))
acf(R_bk_1015)
pacf(R_bk_1015)
par(mfrow = c(2,1))
acf(R_bk_1520,lag = 90)
pacf(R_bk_1520,lag = 90)
#pairs(cbind(R_bk_1520,Fin_sep_log))
#pairs(cbind(R_bk_1520,Fin_sep))
#Normality: Reject H0
jbTest(R_bk) 
shapiroTest(R_bk)
#Heteroskedasticity: Reject H0. Heteroskedasticity exists
par(mfrow = c(1,1))
McLeod.Li.test(y = R_bk) 
library(forecast)
auto.arima(R_bk) #arima(0,0,0)
auto.arima(R_bk_1015) #arima(0,0,0)
armamodel = auto.arima(R_bk_1520,seasonal = TRUE) #arima(1,0,1)
plot(residuals(armamodel))
McLeod.Li.test(y = residuals(armamodel))  #ARCH Effect of arma residuals: Significant

#######################################################################################
#Fintech Factor Plugged in (2015--2020)
#Find the best linear model of R_bk_1520
#PART I
#1.According to acf of R_bk_1520, 
#try R_bk_1520_{t} ~ R_bk_1520_{t-1}+R_bk_1520_{t-2}+R_bk_1520_{t-3}+Fintech Factors(log)+Mixed Index
MA_Fin_sep_log = as.matrix(Fin_sep_log)
FIN_reg1 = MA_Fin_sep_log[-1,]
FIN_reg3 = MA_Fin_sep_log[-(1:3),]
Rbk_lags = matrix(0,nrow = (length(R_bk_1520)-3), ncol = 3)
for(j in 1:3){
  Rbk_lags[,j] = R_bk_1520[j:(length(R_bk_1520)+j-4)] #col 1: lag=3, col 2: lag=2, col 3: lag=1 
}
Rbk_lags <- Rbk_lags[,c(3,2,1)]
#FIT and FIT_ols are the same. To use "ols_step_best_subset" to check regressors, we create FIT_ols
FIT = lm(R_bk_1520[-(1:3)]~Rbk_lags+FIN_reg3+MIX1520[-(1:3)])
summary(FIT)
#Stepwise Regression: Rbk_lags1(lag = 3),MIX1520[-(1:3)]
FIT_step<-step(FIT)
summary(FIT_step)
#Forward Selection: 
#Best Results: R_bk_1520_{t} ~ R_bk_1520_{t-1}+R_bk_1520_{t-2}+R_bk_1520_{t-3}+ATM+POS+IT+Mixed Index
FIT_ols = lm(R_bk_1520[-(1:3)]~Rbk_lags[,1]+Rbk_lags[,2]+Rbk_lags[,3]+FIN_reg3[,1]+FIN_reg3[,2]+FIN_reg3[,3]+FIN_reg3[,4]+FIN_reg3[,5]+MIX1520[-(1:3)])
library(olsrr)
ols_step_best_subset(FIT_ols) 

#2.Acording to auto.arima(1,0,1) of R_bk_1520,
#try R_bk_1520_{t} ~ R_bk_1520_{t-1}+Fintech Factors(log)+Mixed Index
FIT2 = lm(R_bk_1520[-1]~R_bk_1520[-length(R_bk_1520)]+FIN_reg1[,1]+FIN_reg1[,2]+FIN_reg1[,3]+FIN_reg1[,4]+FIN_reg1[,5]+MIX1520[-1])
summary(FIT2)
#Stepwise Regression: FIN_reg1[, 4], MIX1520[-1]
FIT2_step<-step(FIT2)
summary(FIT2_step)
#Forward Selection: 
#Best Result: R_bk_1520_{t} ~ATM+POS+IT+Mixed Index
ols_step_best_subset(FIT2)
##################################################################################################################################3
#GARCH Families Fitting & VaR Computation
m1=cbind(FIN_reg3[,1],FIN_reg3[,2]) #variance$external.regressors
#m2=cbind(Rbk_lags[,3],FIN_reg3[,3],FIN_reg3[,4],FIN_reg3[,5],MIX1520[-(1:3)]) #mean$external.regressors
m3=cbind(FIN_reg3[,3],FIN_reg3[,4],FIN_reg3[,5],MIX1520[-(1:3)]) 
ugarch_sol = list("nlminb", "solnp", "lbfgs", "gosolnp", "nloptr", "hybrid")
ugarch_dist = list("norm","snorm","std","sstd")#,"ged","sged","nig","ghyp","jsu")
set.seed(1)
##############################################################################
#egarch fitting
egarch_spec=ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1), external.regressors = as.matrix(m1),variance.targeting = FALSE),
                       mean.model = list(armaOrder = c(1, 1), include.mean = TRUE, archm = FALSE, archpow = 1, arfima = FALSE, external.regressors = as.matrix(m3), archex = FALSE),
                       distribution.model = "sstd")#distribution.model = dist)  #as.matrix(FIN_reg3[,c(1,2)])   
egarch_fit=ugarchfit(egarch_spec,data=as.data.frame(R_bk_1520[-(1:3)]),solver="hybrid")
egarch_fit
plot(egarch_fit,which="all")
VaR_egarch = fitted(egarch_fit) + sigma(egarch_fit)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew  = coef(egarch_fit)["skew"], shape=coef(egarch_fit)["shape"])
VaRTest(0.05, as.numeric(R_bk_1520[-(1:3)]), as.numeric(VaR_egarch))

#garch fitting
garch_spec=ugarchspec(variance.model=list(model="fGARCH",garchOrder=c(1,0),submodel="GARCH", external.regressors=as.matrix(m1)),
                      mean.model=list(armaOrder=c(1,0),external.regressors=as.matrix(m3)),distribution.model = "sstd")
garch_fit=ugarchfit(garch_spec,data=as.data.frame(R_bk_1520[-(1:3)]),solver="hybrid")
garch_fit
plot(garch_fit,which="all")
VaR_garch = fitted(garch_fit) + sigma(garch_fit)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew  = coef(garch_fit)["skew"], shape=coef(garch_fit)["shape"])
VaRTest(0.05, as.numeric(R_bk_1520[-(1:3)]), as.numeric(VaR_garch))

#aparch fitting
aparch_spec=ugarchspec(variance.model = list(model = "fGARCH", submodel = "APARCH", garchOrder = c(1, 1), external.regressors = as.matrix(m1), variance.targeting = FALSE),
                       mean.model = list(armaOrder = c(1, 1), include.mean = TRUE, archm = FALSE, archpow = 1, arfima = FALSE, external.regressors = as.matrix(m3), archex = FALSE),
                       distribution.model = "std") #
aparch_fit=ugarchfit(aparch_spec,data=R_bk_1520[-(1:3)],solver="hybrid", solver.control = list(),
                 fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all',
                                    trunclag = 1000),
                 numderiv.control = list(grad.eps=1e-4, grad.d=0.0001,
                                         grad.zero.tol=sqrt(.Machine$double.eps/7e-7), hess.eps=1e-4, hess.d=0.1,
                                         hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
aparch_fit
plot(aparch_fit,which="all")
VaR_aparch = fitted(aparch_fit) + sigma(aparch_fit)*qdist("std", p=0.05, mu = 0, sigma = 1, skew  = coef(aparch_fit)["skew"], shape=coef(aparch_fit)["shape"])
VaRTest(0.05, as.numeric(R_bk_1520[-(1:3)]), as.numeric(VaR_aparch))

#gjrgarch fitting
gjrgarch_spec= ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,0), external.regressors = m1, variance.targeting = FALSE),
                          mean.model = list(armaOrder = c(1,0), include.mean = TRUE, archm = FALSE, archpow = 1, arfima = FALSE, external.regressors = m3, archex = FALSE),
                          distribution.model = "std", start.pars = list(), fixed.pars = list())
gjrgarch_fit=ugarchfit(gjrgarch_spec,data=as.data.frame(R_bk_1520[-(1:3)]),solver="hybrid")
gjrgarch_fit
plot(gjrgarch_fit,which="all")
VaR_gjrgarch = fitted(gjrgarch_fit) + sigma(gjrgarch_fit)*qdist("std", p=0.05, mu = 0, sigma = 1, skew  = coef(gjrgarch_fit)["skew"], shape=coef(gjrgarch_fit)["shape"])
VaRTest(0.05, as.numeric(R_bk_1520[-(1:3)]), as.numeric(VaR_gjrgarch))

#tgarch fitting
tgarch_spec= ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(1, 0),
                                       submodel = "TGARCH", external.regressors = m1, variance.targeting = FALSE),
                 mean.model = list(armaOrder = c(1, 1), include.mean = TRUE, archm = FALSE,
                                   archpow = 1, arfima = FALSE, external.regressors = m3, archex = FALSE),
                 distribution.model = "std", start.pars = list(), fixed.pars = list())
tgarch_fit=ugarchfit(tgarch_spec,data=R_bk_1520[-(1:3)],solver="hybrid", solver.control = list())
tgarch_fit
plot(tgarch_fit,which="all")
VaR_tgarch = fitted(tgarch_fit) + sigma(tgarch_fit)*qdist("std", p=0.05, mu = 0, sigma = 1, skew  = coef(tgarch_fit)["skew"], shape=coef(tgarch_fit)["shape"])
VaRTest(0.05, as.numeric(R_bk_1520[-(1:3)]), as.numeric(VaR_tgarch))

tab<-as.data.frame(cbind(VaR_garch,VaR_egarch,VaR_aparch,VaR_tgarch,VaR_gjrgarch))
##############################################################################################################################################
#VaR Forecast (20180516--20200630)
#egarch(1,1)ARMA(1,1) 
egarch_roll = ugarchroll(egarch_spec, as.data.frame(R_bk_1520[-(1:3)]), n.ahead = 1, forecast.length = 500,
           n.start = NULL, refit.every = 30, refit.window = c("recursive", "moving"),
           window.size = NULL, solver = "hybrid", fit.control = list(),
           solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01,0.05),
           cluster = NULL, keep.coef = TRUE)
report(egarch_roll, type = "fpm")
report(egarch_roll, type="VaR", VaR.alpha = 0.01, conf.level = 0.95)
#report(egarch_roll, type="VaR", VaR.alpha = 0.05, conf.level = 0.9)
plot(egarch_roll,which="all",main = "ARMA(1,1)-EGARCH(1,1)")
plot(egarch_roll,which=3,VaR.alpha=0.01)
plot(egarch_roll,which=4,VaR.alpha=0.05)

#garch(1,0) ARMA(1,0)
garch_roll = ugarchroll(garch_spec, as.data.frame(R_bk_1520[-(1:3)]), n.ahead = 1, forecast.length = 500,
                         n.start = NULL, refit.every = 30, refit.window = c("recursive", "moving"),
                         window.size = NULL, solver = "hybrid", fit.control = list(),
                         solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01,0.05),
                         cluster = NULL, keep.coef = TRUE)
report(garch_roll, type = "fpm")
report(garch_roll, type="VaR", VaR.alpha = 0.05, conf.level = 0.95)
plot(garch_roll,which="all",main = "ARMA(1,0)-GARCH(1,0)")
plot(garch_roll,which=3,VaR.alpha=0.01)
plot(garch_roll,which=4,VaR.alpha=0.01)

#aparch(1,1) ARMA(1,1)
aparch_roll = ugarchroll(aparch_spec, as.data.frame(R_bk_1520[-(1:3)]), n.ahead = 1, forecast.length = 500,
                        n.start = NULL, refit.every = 30, refit.window = c("recursive", "moving"),
                        window.size = NULL, solver = "hybrid", fit.control = list(),
                        solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01,0.05),
                        cluster = NULL, keep.coef = TRUE) #hybrid, solnp,lbfgs,gosolnp
report(aparch_roll, type = "fpm")
report(aparch_roll, type="VaR", VaR.alpha = 0.05, conf.level = 0.95)
#report(aparch_roll, type="VaR", VaR.alpha = 0.01, conf.level = 0.95) #下标出界
plot(aparch_roll,which="all")
plot(aparch_roll,which=3,VaR.alpha=0.01)
plot(aparch_roll,which=4,VaR.alpha=0.05)

#arma(1,0)-gjrgarch(1,0) std
gjrgarch_roll = ugarchroll(gjrgarch_spec, as.data.frame(R_bk_1520[-(1:3)]), n.ahead = 1, forecast.length = 500,
                         n.start = NULL, refit.every = 30, refit.window = c("recursive", "moving"),
                         window.size = NULL, solver = "hybrid", fit.control = list(),
                         solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01,0.05,0.1),
                         cluster = NULL, keep.coef = TRUE)#hybrid
report(gjrgarch_roll, type = "fpm")
#report(gjrgarch_roll, type="VaR", VaR.alpha = 0.01, conf.level = 0.95) #出界
report(gjrgarch_roll, type="VaR", VaR.alpha = 0.05, conf.level = 0.95)
plot(gjrgarch_roll,which=3,VaR.alpha=0.01)
plot(gjrgarch_roll,which=4,VaR.alpha=0.05)

#arma(1,1)-tgarch(1,0)
tgarch_roll = ugarchroll(tgarch_spec, as.data.frame(R_bk_1520[-(1:3)]), n.ahead = 1, forecast.length = 500,
                        n.start = NULL, refit.every = 30, refit.window = c("recursive", "moving"),
                        window.size = NULL, solver = "hybrid", fit.control = list(),
                        solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01,0.05),
                        cluster = NULL, keep.coef = TRUE)
report(tgarch_roll, type = "fpm")
report(tgarch_roll, type="VaR", VaR.alpha = 0.05, conf.level = 0.95)
plot(tgarch_roll,which="all")
plot(tgarch_roll,which=3,VaR.alpha=0.01)
plot(tgarch_roll,which=4,VaR.alpha=0.01)