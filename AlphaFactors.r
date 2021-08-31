s = Sys.time()
setwd("D:/rcode/qf5206/")
#VFactor
VFactor = function(X,Y){
  N = ncol(X)
  fct = X/Y
  V = fct
  for(i in 1:nrow(X)){
    V[i,] = rank(fct[i,],ties.method="min")
    V[i,] = (2*V[i,]-N-1)/N
  }
  return(V)
}
#MFactor
MFactor = function(X){
  M = nrow(X)
  N = ncol(X)
  Mfct = X[22:M,]
  for(i in 1:(M-21)){
    X_sub = X[i:(i+20),]
    for(j in 1:N){
      Mfct[i,j] = max(abs(X_sub[,j]))
    }
    Mfct[i,] = rank(Mfct[i,],ties.method="min")
    Mfct[i,] = (2*Mfct[i,]-N-1)/N
  }
  return(Mfct)
}
############################################################################################
in_universe = read.csv("in_univ.csv",header = FALSE)
tickers = read.csv("tickers.csv",header = FALSE)
Date = read.csv("dateline.csv",header = FALSE)
adjp = read.csv("adjusted.csv",header = FALSE)
colnames(tickers) = c("Stock_Name","GICS","NAICS")
dimnames(in_universe) = list(as.character(c(2002:2020)),tickers$Stock_Name)
dimnames(adjp) = list(as.character(Date$V1),tickers$Stock_Name)
index_nagics = which(is.na(tickers$GICS))
tickers_gics = tickers[-index_nagics,]
adjp_gics = adjp[,-index_nagics]
uni_gics = in_universe[,-index_nagics]
GICS = substr(tickers_gics$GICS,1,6)
R_daily = log(adjp_gics[2:nrow(adjp_gics),]/adjp_gics[1:(nrow(adjp_gics)-1),])
R_10days = log(adjp_gics[11:nrow(adjp_gics),]/adjp_gics[1:(nrow(adjp_gics)-10),])
#############################################################################
#Calculate the volatility of prior 10days returns over the past 21 days
sigma_10days = array(dim = c(nrow(R_10days)-21,ncol(R_10days)),
                     dimnames = list(rownames(R_10days[22:nrow(R_10days),]),colnames(R_10days)))
for(j in 1:(nrow(R_10days)-21)){
  for(k in 1:ncol(R_10days)){
    sigma_10days[j,k] = sd(R_10days[j:(j+20),k])
  }
}
# write.csv(sigma_10days,"assign3_sigma_10days.csv")
# sigma_10days = read.csv("assign3_sigma_10days.csv",header = T)
# rownames(sigma_10days)<-sigma_10days[,1]
# sigma_10days<-sigma_10days[,-1]
##############################################################################
dateline = c("20060103","20070103","20080102","20090102","20100104",
             "20110103","20120103","20130102","20140102","20150102",
             "20160104","20170103","20180102","20190102","20200102","20200630")
beta_V_10days<-vector()
beta_M_10days<-vector()
sg_V_10days<-vector()
sg_M_10days<-vector()
t_V_10days<-vector()
t_M_10days<-vector()
TAR_10d<-vector()
TAR_10d_c<-vector()
sigma_annual_10d<-vector()
sigma_annual_10d_c<-vector()
Info_10d = list(short = '',long = '')
#To store the information of the last trading day last year (Enddate of 2007---Enddate of 2020(though no use for 2020))
#used for long/short on the first trading day of the next year
Info_endyear = c(list(NULL),14) 
cost = 5*10^(-4)

dim1 = c("Intercept","Beta_V","Beta_M")
ll = 1
for(yy in 2006:2020){
  index_uni_y = which(uni_gics[rownames(uni_gics)==yy,]==1)
  GICS_y = GICS[index_uni_y]
  freq_y = table(GICS_y)
  dat = Date$V1[which(Date$V1==dateline[ll]):(which(Date$V1==dateline[ll+1])-1)]
  dat2 = Date$V1[(which(Date$V1==dateline[ll])-21):(which(Date$V1==dateline[ll+1])-1)]
  length1 = length(dat)
  length2 = length(freq_y)
  sigma_10days_y = sigma_10days[which(rownames(sigma_10days)==dateline[ll])
                                :(which(rownames(sigma_10days)==dateline[ll+1])-1),index_uni_y]
  sigma_10days_y[sigma_10days_y==0]<-10^(-6)
  R_10days_y = R_10days[which(rownames(R_10days)==dateline[ll]):(which(rownames(R_10days)==dateline[ll+1])-1),index_uni_y]
  R_daily_y = R_daily[(which(rownames(R_daily)==dateline[ll])-21):(which(rownames(R_daily)==dateline[ll+1])-1),index_uni_y]
  R_10days_y[is.na(R_10days_y)] = 0
  R_daily_y[is.na(R_daily_y)] = 0
  R_10days_ind = array(dim = c(length1,length2),dimnames = list(dat,names(freq_y)))
  R_daily_ind = array(dim = c(length1+21,length2),dimnames = list(dat2,names(freq_y)))
  R_10days_sub = R_10days_y
  R_daily_sub = R_daily_y
  #Compute subtracted returns
  for(j in 1:length2){
    I = which(GICS_y==names(freq_y[j]))
    for(k in 1:length1){
      R_10days_ind[k,j] = sum(R_10days_y[k,I])/length(I)
      R_10days_sub[k,I] = R_10days_y[k,I]-R_10days_ind[k,j]
    }
    for(h in 1:(length1+21)){
      R_daily_ind[h,j] = sum(R_daily_y[h,I])/length(I)
      R_daily_sub[h,I] = R_daily_y[h,I]-R_daily_ind[h,j]
    }
  }
  #Compute V, M factors
  V_10days = VFactor(R_10days_sub,sigma_10days_y)
  MFct = MFactor(R_daily_sub)
  R_daily_reg = R_daily_sub[22:nrow(R_daily_sub),]
  coeff_10days = array(dim = c(length1-1,3),dimnames = list(dat[1:(length(dat)-1)],dim1))
  for(i in 1:(length1-1)){
    fit1 <- lm(t(R_daily_reg[i+1,])~t(V_10days[i,])+t(MFct[i,]))
    coeff_10days[i,] = fit1$coefficients
  }
  avg_V_10days = mean(coeff_10days[,2])
  avg_M_10days = mean(coeff_10days[,3])
  sd_V_10days = sd(coeff_10days[,2])
  sd_M_10days = sd(coeff_10days[,3])
  beta_V_10days<-c(beta_V_10days,avg_V_10days)
  beta_M_10days<-c(beta_M_10days,avg_M_10days)
  sg_V_10days<-c(sg_V_10days,sd_V_10days)
  sg_M_10days<-c(sg_M_10days,sd_M_10days)
  t_V_10days<-c(t_V_10days,sqrt(251)*avg_V_10days/sd_V_10days)
  t_M_10days<-c(t_M_10days,sqrt(251)*avg_M_10days/sd_M_10days)
  
###############################Part B######################################################
  if(ll>1){ 
    port_10d = c(rep(list(NULL),length1))
    E_10d = V_10days*beta_V_10days[ll-1]+MFct*beta_M_10days[ll-1] #Calculate the expected returns
    BottomNum = floor(length(index_uni_y)*0.2)
    R_port_10d = matrix(0,length1,1)
    R_port_10d_c = matrix(0,length1,1)
    R_a_10d = 1
    R_a_10d_c = 1
    for(i in 1:(length1)){
      Info_10d$short = which(rank(E_10d[i,],ties.method="min")<=BottomNum)
      Info_10d$long = which(rank(E_10d[i,],ties.method="min")>4*BottomNum)# & rank(E_10d[i,],ties.method="min")<=5*BottomNum)
      LongNum = length(Info_10d$long) #Number of stocks in long position
      ShortNum = length(Info_10d$short)
      #Daily Portfolio Returns with no transaction costs
      R_port_10d[i,] = (sum(E_10d[i,Info_10d$long])-sum(E_10d[i,Info_10d$short]))/LongNum
      port_10d[[i]] = Info_10d
      
      #Daily Portfolio Returns with 5bps transaction costs
      if(i>1){
        # num_10d = 2*(length(Info_10d$long)-length(intersect(port_10d[[i-1]][["long"]],port_10d[[i]][["long"]]))
        #         + length(Info_10d$short)-length(intersect(port_10d[[i-1]][["short"]],port_10d[[i]][["short"]])))
        num_10d = 2*(length(setdiff(port_10d[[i-1]][["long"]],port_10d[[i]][["long"]]))
                     +length(setdiff(port_10d[[i-1]][["short"]],port_10d[[i]][["short"]])))
        R_port_10d_c[i,] = (sum(E_10d[i,Info_10d$long])-sum(E_10d[i,Info_10d$short])-cost*num_10d)/LongNum 
      }
      else{
        #Calculate the daily portfolio return (5bps costs) on the first trading day for each year
        if(ll==2){ #yy==2007
          R_port_10d_c[i,] = (sum(E_10d[i,Info_10d$long])-sum(E_10d[i,Info_10d$short])-cost*(LongNum+ShortNum))/LongNum
        }
        else{
          #num_10d = 2*(length(setdiff(Info_endyear[[ll-2]][["long"]],port_10d[[i]][["long"]]))+length(setdiff(Info_endyear[[ll-2]][["short"]],port_10d[[i]][["short"]])))
          num_10d = LongNum+length(Info_endyear[[ll-2]][["long"]])
                    -2*length(intersect(Info_endyear[[ll-2]][["long"]],port_10d[[i]][["long"]]))
                    + ShortNum+length(Info_endyear[[ll-2]][["short"]])
                    -2*length(intersect(Info_endyear[[ll-2]][["short"]],port_10d[[i]][["short"]]))
          R_port_10d_c[i,] = (sum(E_10d[i,Info_10d$long])-sum(E_10d[i,Info_10d$short])-cost*num_10d)/LongNum 
        }
      }
      R_a_10d = R_a_10d*(1+R_port_10d[i,])
      R_a_10d_c = R_a_10d_c*(1+R_port_10d_c[i,])
    }
    TAR_10d<-c(TAR_10d,R_a_10d)
    TAR_10d_c<-c(TAR_10d_c,R_a_10d_c)
    sigma_annual_10d<-c(sigma_annual_10d,sd(R_port_10d))
    sigma_annual_10d_c<-c(sigma_annual_10d_c,sd(R_port_10d_c))
    Info_endyear[[ll-1]] = port_10d[[length1]]
  }
  ll = ll+1
}
table_10d = cbind(beta_V_10days,sg_V_10days,t_V_10days,
                  beta_M_10days,sg_M_10days,t_M_10days,
                  TAR_10d=c(NA,TAR_10d-1),sigma_annual_10d=c(NA,sigma_annual_10d),
                  TAR_10d_c=c(NA,TAR_10d_c-1),sigma_annual_10d_c=c(NA,sigma_annual_10d_c))
rownames(table_10d) = c(2006:2020)
write.csv(table_10d,"Assign3_10days_final.csv")
e = Sys.time()
print(e-s)