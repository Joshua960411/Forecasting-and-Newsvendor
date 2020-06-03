rm(list=ls())
library('forecast')
library('smooth')
library('quantreg')
library('doMC')
registerDoMC(16)
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}
iter<-20000
#dgp

##ar
arima_series<-sim.ssarima(frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),obs = 5001,nsim=iter,constant = 500,mean=0,sd=200)$data

#function

##linear
profit_function_linear<-function(order,demand){
  profit<-(order>=demand)*(17*demand-5*order)+(order<demand)*(5*order+7*demand)
  return(profit)
}

quant_linear<-1/2

##estimator
test_length=1

mini_linear1<-function(data,par){
  order<-par[1]+par[2]*data[,5]
  total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function_linear(x[1],x[2])))
  return(-total_profit)
}

mini_linear2<-function(data,par){
  order<-par[1]+par[2]*data[,2]+par[3]*data[,3]+par[4]*data[,4]+par[5]*data[,5]+par[6]*data[,6]+par[7]*data[,7]+par[8]*data[,8]
  total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function_linear(x[1],x[2])))
  return(-total_profit)
}

mini_linear3<-function(data,par){
  order<-par*data[,5]+par*(1-par)*data[,6]+par*(1-par)^2*data[,7]+par*(1-par)^3*data[,8]
  total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function_linear(x[1],x[2])))
  return(-total_profit)
}


##generate function
ppl_l<-function(series,forecast){
  ((profit_function_linear(series,series)-profit_function_linear(forecast,series))/profit_function_linear(series,series))
}

##40
series<-arima_series[201:240,]

a_40<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  s1<-rep(c(1,0,0,0),rep)
  s2<-rep(c(0,1,0,0),rep)
  s3<-rep(c(0,0,1,0),rep)
  L1_data=lag(data,-1)
  L2_data=lag(data,-2)
  L3_data=lag(data,-3)
  L4_data=lag(data,-4)
  set_data=cbind(data,s1,s2,s3)
  set_data=cbind(set_data,L1_data,L2_data,L3_data,L4_data)
  colnames(set_data)<-c('data','s1','s2','s3','L1','L2','L3','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length,interval="parametric",level=quant_linear*2-1)$upper
  arima_u1<-arima(test,order=c(1,0,0))
  arima_u1_l<-forecast(arima_u1,test_length)$mean+qnorm(quant_linear,mean=0,sd=sd(arima_u1$residuals))
  arima_u2<-ssarima(test,frequency=4,orders=list(ar=c(3,1)),lags = c(1,4),constant = TRUE)
  arima_u2_l<-forecast(arima_u2,test_length,interval="parametric",level=quant_linear*2-1)$upper
  es_u_l<-es(test,'ANN')
  es_u<-forecast(es_u_l,test_length)$mean+qnorm(quant_linear,mean=0,sd=sd(es_u_l$residuals))
  #cf
  coe1<-lm(data ~ L1, data=set_data)$coefficients
  cf_l1_par<-optim(par = coe1, fn = mini_linear1, method = 'L-BFGS-B', data = set_data)$par
  cf_l1<-cf_l1_par[1]+cf_l1_par[2]*data[40]
  coe2<-lm(data ~., data=set_data)$coefficients
  cf_l2_par<-optim(par = coe2, fn = mini_linear2, method = 'L-BFGS-B', data = set_data)$par
  cf_l2<-cf_l2_par[1]+cf_l2_par[2]*1+cf_l2_par[3]*0+cf_l2_par[4]*0+cf_l2_par[5]*data[40]+cf_l2_par[6]*data[39]+cf_l2_par[7]*data[38]+cf_l2_par[8]*data[37]
  cf_l3_par<-optim(par = 0.5, fn = mini_linear3, method = 'L-BFGS-B', data = set_data)$par
  cf_l3<-cf_l3_par*data[40]+cf_l3_par*(1-cf_l3_par)*data[39]+cf_l3_par*(1-cf_l3_par)^2*data[38]+cf_l3_par*(1-cf_l3_par)^3*data[37]
  ##list
  list(arima_k_l,arima_u1_l,arima_u2_l,es_u,cf_l1,cf_l2,cf_l3)
}

##120
series<-arima_series[201:320,]

a_120<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  s1<-rep(c(1,0,0,0),rep)
  s2<-rep(c(0,1,0,0),rep)
  s3<-rep(c(0,0,1,0),rep)
  L1_data=lag(data,-1)
  L2_data=lag(data,-2)
  L3_data=lag(data,-3)
  L4_data=lag(data,-4)
  set_data=cbind(data,s1,s2,s3)
  set_data=cbind(set_data,L1_data,L2_data,L3_data,L4_data)
  colnames(set_data)<-c('data','s1','s2','s3','L1','L2','L3','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length,interval="parametric",level=quant_linear*2-1)$upper
  arima_u1<-arima(test,order=c(1,0,0))
  arima_u1_l<-forecast(arima_u1,test_length)$mean+qnorm(quant_linear,mean=0,sd=sd(arima_u1$residuals))
  arima_u2<-ssarima(test,frequency=4,orders=list(ar=c(3,1)),lags = c(1,4),constant = TRUE)
  arima_u2_l<-forecast(arima_u2,test_length,interval="parametric",level=quant_linear*2-1)$upper
  es_u_l<-es(test,'ANN')
  es_u<-forecast(es_u_l,test_length)$mean+qnorm(quant_linear,mean=0,sd=sd(es_u_l$residuals))
  #cf
  coe1<-lm(data ~ L1, data=set_data)$coefficients
  cf_l1_par<-optim(par = coe1, fn = mini_linear1, method = 'L-BFGS-B', data = set_data)$par
  cf_l1<-cf_l1_par[1]+cf_l1_par[2]*data[120]
  coe2<-lm(data ~., data=set_data)$coefficients
  cf_l2_par<-optim(par = coe2, fn = mini_linear2, method = 'L-BFGS-B', data = set_data)$par
  cf_l2<-cf_l2_par[1]+cf_l2_par[2]*1+cf_l2_par[3]*0+cf_l2_par[4]*0+cf_l2_par[5]*data[120]+cf_l2_par[6]*data[119]+cf_l2_par[7]*data[118]+cf_l2_par[8]*data[117]
  cf_l3_par<-optim(par = 0.5, fn = mini_linear3, method = 'L-BFGS-B', data = set_data)$par
  cf_l3<-cf_l3_par*data[120]+cf_l3_par*(1-cf_l3_par)*data[119]+cf_l3_par*(1-cf_l3_par)^2*data[118]+cf_l3_par*(1-cf_l3_par)^3*data[117]
  ##list
  list(arima_k_l,arima_u1_l,arima_u2_l,es_u,cf_l1,cf_l2,cf_l3)
}

##480
series<-arima_series[201:680,]

a_480<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  s1<-rep(c(1,0,0,0),rep)
  s2<-rep(c(0,1,0,0),rep)
  s3<-rep(c(0,0,1,0),rep)
  L1_data=lag(data,-1)
  L2_data=lag(data,-2)
  L3_data=lag(data,-3)
  L4_data=lag(data,-4)
  set_data=cbind(data,s1,s2,s3)
  set_data=cbind(set_data,L1_data,L2_data,L3_data,L4_data)
  colnames(set_data)<-c('data','s1','s2','s3','L1','L2','L3','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length,interval="parametric",level=quant_linear*2-1)$upper
  arima_u1<-arima(test,order=c(1,0,0))
  arima_u1_l<-forecast(arima_u1,test_length)$mean+qnorm(quant_linear,mean=0,sd=sd(arima_u1$residuals))
  arima_u2<-ssarima(test,frequency=4,orders=list(ar=c(3,1)),lags = c(1,4),constant = TRUE)
  arima_u2_l<-forecast(arima_u2,test_length,interval="parametric",level=quant_linear*2-1)$upper
  es_u_l<-es(test,'ANN')
  es_u<-forecast(es_u_l,test_length)$mean+qnorm(quant_linear,mean=0,sd=sd(es_u_l$residuals))
  #cf
  coe1<-lm(data ~ L1, data=set_data)$coefficients
  cf_l1_par<-optim(par = coe1, fn = mini_linear1, method = 'L-BFGS-B', data = set_data)$par
  cf_l1<-cf_l1_par[1]+cf_l1_par[2]*data[480]
  coe2<-lm(data ~., data=set_data)$coefficients
  cf_l2_par<-optim(par = coe2, fn = mini_linear2, method = 'L-BFGS-B', data = set_data)$par
  cf_l2<-cf_l2_par[1]+cf_l2_par[2]*1+cf_l2_par[3]*0+cf_l2_par[4]*0+cf_l2_par[5]*data[480]+cf_l2_par[6]*data[479]+cf_l2_par[7]*data[478]+cf_l2_par[8]*data[477]
  cf_l3_par<-optim(par = 0.5, fn = mini_linear3, method = 'L-BFGS-B', data = set_data)$par
  cf_l3<-cf_l3_par*data[480]+cf_l3_par*(1-cf_l3_par)*data[479]+cf_l3_par*(1-cf_l3_par)^2*data[478]+cf_l3_par*(1-cf_l3_par)^3*data[477]
  ##list
  list(arima_k_l,arima_u1_l,arima_u2_l,es_u,cf_l1,cf_l2,cf_l3)
}

##1200
series<-arima_series[201:1400,]

a_1200<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  s1<-rep(c(1,0,0,0),rep)
  s2<-rep(c(0,1,0,0),rep)
  s3<-rep(c(0,0,1,0),rep)
  L1_data=lag(data,-1)
  L2_data=lag(data,-2)
  L3_data=lag(data,-3)
  L4_data=lag(data,-4)
  set_data=cbind(data,s1,s2,s3)
  set_data=cbind(set_data,L1_data,L2_data,L3_data,L4_data)
  colnames(set_data)<-c('data','s1','s2','s3','L1','L2','L3','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length,interval="parametric",level=quant_linear*2-1)$upper
  arima_u1<-arima(test,order=c(1,0,0))
  arima_u1_l<-forecast(arima_u1,test_length)$mean+qnorm(quant_linear,mean=0,sd=sd(arima_u1$residuals))
  arima_u2<-ssarima(test,frequency=4,orders=list(ar=c(3,1)),lags = c(1,4),constant = TRUE)
  arima_u2_l<-forecast(arima_u2,test_length,interval="parametric",level=quant_linear*2-1)$upper
  es_u_l<-es(test,'ANN')
  es_u<-forecast(es_u_l,test_length)$mean+qnorm(quant_linear,mean=0,sd=sd(es_u_l$residuals))
  #cf
  coe1<-lm(data ~ L1, data=set_data)$coefficients
  cf_l1_par<-optim(par = coe1, fn = mini_linear1, method = 'L-BFGS-B', data = set_data)$par
  cf_l1<-cf_l1_par[1]+cf_l1_par[2]*data[1200]
  coe2<-lm(data ~., data=set_data)$coefficients
  cf_l2_par<-optim(par = coe2, fn = mini_linear2, method = 'L-BFGS-B', data = set_data)$par
  cf_l2<-cf_l2_par[1]+cf_l2_par[2]*1+cf_l2_par[3]*0+cf_l2_par[4]*0+cf_l2_par[5]*data[1200]+cf_l2_par[6]*data[1199]+cf_l2_par[7]*data[1198]+cf_l2_par[8]*data[1197]
  cf_l3_par<-optim(par = 0.5, fn = mini_linear3, method = 'L-BFGS-B', data = set_data)$par
  cf_l3<-cf_l3_par*data[1200]+cf_l3_par*(1-cf_l3_par)*data[1199]+cf_l3_par*(1-cf_l3_par)^2*data[1198]+cf_l3_par*(1-cf_l3_par)^3*data[1197]
  ##list
  list(arima_k_l,arima_u1_l,arima_u2_l,es_u,cf_l1,cf_l2,cf_l3)
}

##4800
series<-arima_series[201:5000,]

a_4800<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  s1<-rep(c(1,0,0,0),rep)
  s2<-rep(c(0,1,0,0),rep)
  s3<-rep(c(0,0,1,0),rep)
  L1_data=lag(data,-1)
  L2_data=lag(data,-2)
  L3_data=lag(data,-3)
  L4_data=lag(data,-4)
  set_data=cbind(data,s1,s2,s3)
  set_data=cbind(set_data,L1_data,L2_data,L3_data,L4_data)
  colnames(set_data)<-c('data','s1','s2','s3','L1','L2','L3','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length,interval="parametric",level=quant_linear*2-1)$upper
  arima_u1<-arima(test,order=c(1,0,0))
  arima_u1_l<-forecast(arima_u1,test_length)$mean+qnorm(quant_linear,mean=0,sd=sd(arima_u1$residuals))
  arima_u2<-ssarima(test,frequency=4,orders=list(ar=c(3,1)),lags = c(1,4),constant = TRUE)
  arima_u2_l<-forecast(arima_u2,test_length,interval="parametric",level=quant_linear*2-1)$upper
  es_u_l<-es(test,'ANN')
  es_u<-forecast(es_u_l,test_length)$mean+qnorm(quant_linear,mean=0,sd=sd(es_u_l$residuals))
  #cf
  coe1<-lm(data ~ L1, data=set_data)$coefficients
  cf_l1_par<-optim(par = coe1, fn = mini_linear1, method = 'L-BFGS-B', data = set_data)$par
  cf_l1<-cf_l1_par[1]+cf_l1_par[2]*data[4800]
  coe2<-lm(data ~., data=set_data)$coefficients
  cf_l2_par<-optim(par = coe2, fn = mini_linear2, method = 'L-BFGS-B', data = set_data)$par
  cf_l2<-cf_l2_par[1]+cf_l2_par[2]*1+cf_l2_par[3]*0+cf_l2_par[4]*0+cf_l2_par[5]*data[4800]+cf_l2_par[6]*data[4799]+cf_l2_par[7]*data[4798]+cf_l2_par[8]*data[4797]
  cf_l3_par<-optim(par = 0.5, fn = mini_linear3, method = 'L-BFGS-B', data = set_data)$par
  cf_l3<-cf_l3_par*data[4800]+cf_l3_par*(1-cf_l3_par)*data[4799]+cf_l3_par*(1-cf_l3_par)^2*data[4798]+cf_l3_par*(1-cf_l3_par)^3*data[4797]
  ##list
  list(arima_k_l,arima_u1_l,arima_u2_l,es_u,cf_l1,cf_l2,cf_l3)
}

#generate
##40
series<-arima_series[241,]

re_40<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list())) %dopar%{
  arima_k_l<-ppl_l(series[i],a_40[[1]][[i]])
  arima_u1_l<-ppl_l(series[i],a_40[[2]][[i]])
  arima_u2_l<-ppl_l(series[i],a_40[[3]][[i]])
  es_u<-ppl_l(series[i],a_40[[4]][[i]])
  cf_l1<-ppl_l(series[i],a_40[[5]][[i]])
  cf_l2<-ppl_l(series[i],a_40[[6]][[i]])
  cf_l3<-ppl_l(series[i],a_40[[7]][[i]])
  list(arima_k_l,arima_u1_l,arima_u2_l,es_u,cf_l1,cf_l2,cf_l3)
}

##120
series<-arima_series[321,]

re_120<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list())) %dopar%{
  arima_k_l<-ppl_l(series[i],a_120[[1]][[i]])
  arima_u1_l<-ppl_l(series[i],a_120[[2]][[i]])
  arima_u2_l<-ppl_l(series[i],a_120[[3]][[i]])
  es_u<-ppl_l(series[i],a_120[[4]][[i]])
  cf_l1<-ppl_l(series[i],a_120[[5]][[i]])
  cf_l2<-ppl_l(series[i],a_120[[6]][[i]])
  cf_l3<-ppl_l(series[i],a_120[[7]][[i]])
  list(arima_k_l,arima_u1_l,arima_u2_l,es_u,cf_l1,cf_l2,cf_l3)
}

##480
series<-arima_series[681,]

re_480<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list())) %dopar%{
  arima_k_l<-ppl_l(series[i],a_480[[1]][[i]])
  arima_u1_l<-ppl_l(series[i],a_480[[2]][[i]])
  arima_u2_l<-ppl_l(series[i],a_480[[3]][[i]])
  es_u<-ppl_l(series[i],a_480[[4]][[i]])
  cf_l1<-ppl_l(series[i],a_480[[5]][[i]])
  cf_l2<-ppl_l(series[i],a_480[[6]][[i]])
  cf_l3<-ppl_l(series[i],a_480[[7]][[i]])
  list(arima_k_l,arima_u1_l,arima_u2_l,es_u,cf_l1,cf_l2,cf_l3)
}

##1200
series<-arima_series[1401,]

re_1200<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list())) %dopar%{
  arima_k_l<-ppl_l(series[i],a_1200[[1]][[i]])
  arima_u1_l<-ppl_l(series[i],a_1200[[2]][[i]])
  arima_u2_l<-ppl_l(series[i],a_1200[[3]][[i]])
  es_u<-ppl_l(series[i],a_1200[[4]][[i]])
  cf_l1<-ppl_l(series[i],a_1200[[5]][[i]])
  cf_l2<-ppl_l(series[i],a_1200[[6]][[i]])
  cf_l3<-ppl_l(series[i],a_1200[[7]][[i]])
  list(arima_k_l,arima_u1_l,arima_u2_l,es_u,cf_l1,cf_l2,cf_l3)
}

##4800
series<-arima_series[5001,]

re_4800<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list())) %dopar%{
  arima_k_l<-ppl_l(series[i],a_4800[[1]][[i]])
  arima_u1_l<-ppl_l(series[i],a_4800[[2]][[i]])
  arima_u2_l<-ppl_l(series[i],a_4800[[3]][[i]])
  es_u<-ppl_l(series[i],a_4800[[4]][[i]])
  cf_l1<-ppl_l(series[i],a_4800[[5]][[i]])
  cf_l2<-ppl_l(series[i],a_4800[[6]][[i]])
  cf_l3<-ppl_l(series[i],a_4800[[7]][[i]])
  list(arima_k_l,arima_u1_l,arima_u2_l,es_u,cf_l1,cf_l2,cf_l3)
}


save(re_40,re_120,re_480,re_1200,re_4800,file='information.Rdata')
