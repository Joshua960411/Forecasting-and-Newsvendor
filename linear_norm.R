rm(list=ls())
gc()
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
arima_series1<-sim.ssarima(frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),obs = 5001,nsim=iter,constant = 500,mean=0,sd=200)$data
arima_series2<-sim.ssarima(frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),obs = 5001,nsim=iter,constant = 500,mean=0,sd=200)$data
arima_series3<-sim.ssarima(frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),obs = 5001,nsim=iter,constant = 500,mean=0,sd=200)$data

#function

##linear
profit_function_linear1<-function(order,demand){
  profit<-(order>=demand)*(17*demand-5*order)+(order<demand)*(5*order+7*demand)
  return(profit)
}

quant_linear1<-1/2

profit_function_linear2<-function(order,demand){
  profit<-(order>=demand)*(23*demand-11*order)+(order<demand)*(19*order-7*demand)
  return(profit)
}

quant_linear2<-19/30

profit_function_linear3<-function(order,demand){
  profit<-(order>=demand)*(17*demand-7*order)+(order<demand)*(3*order+7*demand)
  return(profit)
}

quant_linear3<-3/10

##estimator
test_length=1

mini_linear1<-function(data,par){
  order<-par[1]+par[2]*data[,2]+par[3]*data[,3]
  total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function_linear1(x[1],x[2])))
  return(-total_profit)
}

mini_linear2<-function(data,par){
  order<-par[1]+par[2]*data[,2]+par[3]*data[,3]
  total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function_linear2(x[1],x[2])))
  return(-total_profit)
}

mini_linear3<-function(data,par){
  order<-par[1]+par[2]*data[,2]+par[3]*data[,3]
  total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function_linear3(x[1],x[2])))
  return(-total_profit)
}

##generate function
result_l1<-function(series,forecast){
  ie<-series-forecast
  ppl<-((profit_function_linear1(series,series)-profit_function_linear1(forecast,series))/profit_function_linear1(series,series))
  sl<-forecast>=series
  fr<-min(forecast,series)/series
  list<-c(ie,ppl,sl,fr)
  return(list)
}

result_l2<-function(series,forecast){
  ie<-series-forecast
  ppl<-((profit_function_linear2(series,series)-profit_function_linear2(forecast,series))/profit_function_linear2(series,series))
  sl<-forecast>=series
  fr<-min(forecast,series)/series
  list<-c(ie,ppl,sl,fr)
  return(list)
}

result_l3<-function(series,forecast){
  ie<-series-forecast
  ppl<-((profit_function_linear3(series,series)-profit_function_linear3(forecast,series))/profit_function_linear3(series,series))
  sl<-forecast>=series
  fr<-min(forecast,series)/series
  list<-c(ie,ppl,sl,fr)
  return(list)
}


##40
series1<-arima_series1[201:240,]

a_40_1<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series1[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear1,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear1*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear1)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[40]+qr_l_par[3]*data[37]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear1, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[40]+cf_l_par[3]*data[37]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

series2<-arima_series2[201:240,]

a_40_2<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series2[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear2,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear2*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear2)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[40]+qr_l_par[3]*data[37]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear2, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[40]+cf_l_par[3]*data[37]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

series3<-arima_series3[201:240,]

a_40_3<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series3[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear3,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear3*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear3)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[40]+qr_l_par[3]*data[37]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear3, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[40]+cf_l_par[3]*data[37]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

##120
series1<-arima_series1[201:320,]

a_120_1<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series1[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear1,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear1*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear1)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[120]+qr_l_par[3]*data[117]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear1, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[120]+cf_l_par[3]*data[117]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

series2<-arima_series2[201:320,]

a_120_2<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series2[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear2,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear2*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear2)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[120]+qr_l_par[3]*data[117]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear2, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[120]+cf_l_par[3]*data[117]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

series3<-arima_series3[201:320,]

a_120_3<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series3[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear3,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear3*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear3)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[120]+qr_l_par[3]*data[117]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear3, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[120]+cf_l_par[3]*data[117]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}


##480
series1<-arima_series1[201:680,]

a_480_1<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series1[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear1,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear1*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear1)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[480]+qr_l_par[3]*data[477]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear1, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[480]+cf_l_par[3]*data[477]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

series2<-arima_series2[201:680,]

a_480_2<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series2[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear2,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear2*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear2)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[480]+qr_l_par[3]*data[477]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear2, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[480]+cf_l_par[3]*data[477]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

series3<-arima_series3[201:680,]

a_480_3<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series3[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear3,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear3*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear3)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[480]+qr_l_par[3]*data[477]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear3, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[480]+cf_l_par[3]*data[477]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

##1200
series1<-arima_series1[201:1400,]

a_1200_1<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series1[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear1,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear1*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear1)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[1200]+qr_l_par[3]*data[1197]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear1, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[1200]+cf_l_par[3]*data[1197]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

series2<-arima_series2[201:1400,]

a_1200_2<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series2[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear2,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear2*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear2)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[1200]+qr_l_par[3]*data[1197]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear2, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[1200]+cf_l_par[3]*data[1197]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

series3<-arima_series3[201:1400,]

a_1200_3<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series3[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear3,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear3*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear3)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[1200]+qr_l_par[3]*data[1197]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear3, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[1200]+cf_l_par[3]*data[1197]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

##4800
series1<-arima_series1[201:5000,]

a_4800_1<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series1[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear1,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear1*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear1)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[4800]+qr_l_par[3]*data[4797]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear1, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[4800]+cf_l_par[3]*data[4797]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

series2<-arima_series2[201:5000,]

a_4800_2<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series2[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear2,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear2*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear2)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[4800]+qr_l_par[3]*data[4797]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear2, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[4800]+cf_l_par[3]*data[4797]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}

series3<-arima_series3[201:5000,]

a_4800_3<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series3[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,L1_data,L4_data)
  colnames(set_data)<-c('data','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant = 500,mean=0,sd=200)
  arima_k_l<-forecast(arima_k,test_length)$mean+qnorm(quant_linear3,0,200)
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_l<-forecast(arima_p,test_length,interval="parametric",level=quant_linear3*2-1)$upper
  #qr
  qr_l_par<-rq(data ~ .,data=set_data,tau=quant_linear3)$coefficients
  qr_l<-qr_l_par[1]+qr_l_par[2]*data[4800]+qr_l_par[3]*data[4797]
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_l_par<-optim(par = coe, fn = mini_linear3, method = 'L-BFGS-B', data = set_data)$par
  cf_l<-cf_l_par[1]+cf_l_par[2]*data[4800]+cf_l_par[3]*data[4797]
  ##list
  list(arima_k_l,arima_p_l,qr_l,cf_l)
}


#generate
##40
series1<-arima_series1[241,]
series2<-arima_series2[241,]
series3<-arima_series3[241,]

re_40_p_l<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())) %dopar%{
  arima_k1<-result_l1(series1[i],a_40_1[[1]][[i]])
  arima_p1<-result_l1(series1[i],a_40_1[[2]][[i]])
  qr1<-result_l1(series1[i],a_40_1[[3]][[i]])
  cf1<-result_l1(series1[i],a_40_1[[4]][[i]])
  arima_k2<-result_l2(series2[i],a_40_2[[1]][[i]])
  arima_p2<-result_l2(series2[i],a_40_2[[2]][[i]])
  qr2<-result_l2(series2[i],a_40_2[[3]][[i]])
  cf2<-result_l2(series2[i],a_40_2[[4]][[i]])
  arima_k3<-result_l3(series3[i],a_40_3[[1]][[i]])
  arima_p3<-result_l3(series3[i],a_40_3[[2]][[i]])
  qr3<-result_l3(series3[i],a_40_3[[3]][[i]])
  cf3<-result_l3(series3[i],a_40_3[[4]][[i]])
  list(arima_k1,arima_p1,qr1,cf1,arima_k2,arima_p2,qr2,cf2,arima_k3,arima_p3,qr3,cf3)
}


##120
series1<-arima_series1[321,]
series2<-arima_series2[321,]
series3<-arima_series3[321,]

re_120_p_l<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())) %dopar%{
  arima_k1<-result_l1(series1[i],a_120_1[[1]][[i]])
  arima_p1<-result_l1(series1[i],a_120_1[[2]][[i]])
  qr1<-result_l1(series1[i],a_120_1[[3]][[i]])
  cf1<-result_l1(series1[i],a_120_1[[4]][[i]])
  arima_k2<-result_l2(series2[i],a_120_2[[1]][[i]])
  arima_p2<-result_l2(series2[i],a_120_2[[2]][[i]])
  qr2<-result_l2(series2[i],a_120_2[[3]][[i]])
  cf2<-result_l2(series2[i],a_120_2[[4]][[i]])
  arima_k3<-result_l3(series3[i],a_120_3[[1]][[i]])
  arima_p3<-result_l3(series3[i],a_120_3[[2]][[i]])
  qr3<-result_l3(series3[i],a_120_3[[3]][[i]])
  cf3<-result_l3(series3[i],a_120_3[[4]][[i]])
  list(arima_k1,arima_p1,qr1,cf1,arima_k2,arima_p2,qr2,cf2,arima_k3,arima_p3,qr3,cf3)
}


##480
series1<-arima_series1[681,]
series2<-arima_series2[681,]
series3<-arima_series3[681,]

re_480_p_l<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())) %dopar%{
  arima_k1<-result_l1(series1[i],a_480_1[[1]][[i]])
  arima_p1<-result_l1(series1[i],a_480_1[[2]][[i]])
  qr1<-result_l1(series1[i],a_480_1[[3]][[i]])
  cf1<-result_l1(series1[i],a_480_1[[4]][[i]])
  arima_k2<-result_l2(series2[i],a_480_2[[1]][[i]])
  arima_p2<-result_l2(series2[i],a_480_2[[2]][[i]])
  qr2<-result_l2(series2[i],a_480_2[[3]][[i]])
  cf2<-result_l2(series2[i],a_480_2[[4]][[i]])
  arima_k3<-result_l3(series3[i],a_480_3[[1]][[i]])
  arima_p3<-result_l3(series3[i],a_480_3[[2]][[i]])
  qr3<-result_l3(series3[i],a_480_3[[3]][[i]])
  cf3<-result_l3(series3[i],a_480_3[[4]][[i]])
  list(arima_k1,arima_p1,qr1,cf1,arima_k2,arima_p2,qr2,cf2,arima_k3,arima_p3,qr3,cf3)
}


##1200
series1<-arima_series1[1401,]
series2<-arima_series2[1401,]
series3<-arima_series3[1401,]

re_1200_p_l<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())) %dopar%{
  arima_k1<-result_l1(series1[i],a_1200_1[[1]][[i]])
  arima_p1<-result_l1(series1[i],a_1200_1[[2]][[i]])
  qr1<-result_l1(series1[i],a_1200_1[[3]][[i]])
  cf1<-result_l1(series1[i],a_1200_1[[4]][[i]])
  arima_k2<-result_l2(series2[i],a_1200_2[[1]][[i]])
  arima_p2<-result_l2(series2[i],a_1200_2[[2]][[i]])
  qr2<-result_l2(series2[i],a_1200_2[[3]][[i]])
  cf2<-result_l2(series2[i],a_1200_2[[4]][[i]])
  arima_k3<-result_l3(series3[i],a_1200_3[[1]][[i]])
  arima_p3<-result_l3(series3[i],a_1200_3[[2]][[i]])
  qr3<-result_l3(series3[i],a_1200_3[[3]][[i]])
  cf3<-result_l3(series3[i],a_1200_3[[4]][[i]])
  list(arima_k1,arima_p1,qr1,cf1,arima_k2,arima_p2,qr2,cf2,arima_k3,arima_p3,qr3,cf3)
}


##4800
series1<-arima_series1[5001,]
series2<-arima_series2[5001,]
series3<-arima_series3[5001,]

re_4800_p_l<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list(),list())) %dopar%{
  arima_k1<-result_l1(series1[i],a_4800_1[[1]][[i]])
  arima_p1<-result_l1(series1[i],a_4800_1[[2]][[i]])
  qr1<-result_l1(series1[i],a_4800_1[[3]][[i]])
  cf1<-result_l1(series1[i],a_4800_1[[4]][[i]])
  arima_k2<-result_l2(series2[i],a_4800_2[[1]][[i]])
  arima_p2<-result_l2(series2[i],a_4800_2[[2]][[i]])
  qr2<-result_l2(series2[i],a_4800_2[[3]][[i]])
  cf2<-result_l2(series2[i],a_4800_2[[4]][[i]])
  arima_k3<-result_l3(series3[i],a_4800_3[[1]][[i]])
  arima_p3<-result_l3(series3[i],a_4800_3[[2]][[i]])
  qr3<-result_l3(series3[i],a_4800_3[[3]][[i]])
  cf3<-result_l3(series3[i],a_4800_3[[4]][[i]])
  list(arima_k1,arima_p1,qr1,cf1,arima_k2,arima_p2,qr2,cf2,arima_k3,arima_p3,qr3,cf3)
}


save(re_40_p_l,re_120_p_l,re_480_p_l,re_1200_p_l,re_4800_p_l,file='linear_norm.Rdata')

