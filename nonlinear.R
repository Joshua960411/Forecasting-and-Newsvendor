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
arima_series<-sim.ssarima(frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),obs = 5001,nsim=iter,constant=500,mean=0,sd=200)$data

#function

##nonlinear
profit_function_nonlinear<-function(order,demand){
  profit<-(order>=demand)*(24*demand-12*order+5*min((order-demand),rnorm(1,30,5)))+(order<demand)*(12*order-0.01*(demand-order)^2)
  return(profit)
}

profit_function_nonlinear_t<-function(order,demand,sed){
  profit<-(order>=demand)*(24*demand-12*order+5*min((order-demand),sed))+(order<demand)*(12*order-0.01*(demand-order)^2)
  return(profit)
}


requan<-0.56


##estimator
test_length=1

mini_nonlinear<-function(data,par){
  order<-par[1]+par[2]*data[,2]+par[3]*data[,3]+par[4]*data[,4]+par[5]*data[,5]+par[6]*data[,6]
  total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function_nonlinear(x[1],x[2])))
  return(-total_profit)
}


##generate function
result_n<-function(series,forecast,sed){
  ppl<-((profit_function_nonlinear_t(series,series,sed)-profit_function_nonlinear_t(forecast,series,sed))/profit_function_nonlinear_t(series,series,sed))
  sl<-forecast>=series
  list<-c(ppl,sl)
  return(list)
}


##40
series<-arima_series[201:240,]

a_40<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  s1<-rep(c(1,0,0,0),rep)
  s2<-rep(c(0,1,0,0),rep)
  s3<-rep(c(0,0,1,0),rep)
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,s1,s2,s3)
  set_data=cbind(set_data,L1_data,L4_data)
  colnames(set_data)<-c('data','s1','s2','s3','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant=500,mean=0,sd=200)
  arima_k_n<-forecast(arima_k,test_length,interval="parametric",level=requan*2-1)$upper
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_n<-forecast(arima_p,test_length,interval="parametric",level=requan*2-1)$upper
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_n_par<-optim(par = coe, fn = mini_nonlinear,method = 'L-BFGS-B', data = set_data)$par
  cf_n<-cf_n_par[1]+cf_n_par[2]*1+cf_n_par[3]*0+cf_n_par[4]*0+cf_n_par[5]*data[40]+cf_n_par[6]*data[37]
  ##list
  list(arima_k_n,arima_p_n,cf_n)
}

##120
series<-arima_series[201:320,]

a_120<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  s1<-rep(c(1,0,0,0),rep)
  s2<-rep(c(0,1,0,0),rep)
  s3<-rep(c(0,0,1,0),rep)
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,s1,s2,s3)
  set_data=cbind(set_data,L1_data,L4_data)
  colnames(set_data)<-c('data','s1','s2','s3','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant=500,mean=0,sd=200)
  arima_k_n<-forecast(arima_k,test_length,interval="parametric",level=requan*2-1)$upper
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_n<-forecast(arima_p,test_length,interval="parametric",level=requan*2-1)$upper
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_n_par<-optim(par = coe, fn = mini_nonlinear,method = 'L-BFGS-B', data = set_data)$par
  cf_n<-cf_n_par[1]+cf_n_par[2]*1+cf_n_par[3]*0+cf_n_par[4]*0+cf_n_par[5]*data[120]+cf_n_par[6]*data[117]
  ##list
  list(arima_k_n,arima_p_n,cf_n)
}


##480
series<-arima_series[201:680,]

a_480<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  s1<-rep(c(1,0,0,0),rep)
  s2<-rep(c(0,1,0,0),rep)
  s3<-rep(c(0,0,1,0),rep)
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,s1,s2,s3)
  set_data=cbind(set_data,L1_data,L4_data)
  colnames(set_data)<-c('data','s1','s2','s3','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant=500,mean=0,sd=200)
  arima_k_n<-forecast(arima_k,test_length,interval="parametric",level=requan*2-1)$upper
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_n<-forecast(arima_p,test_length,interval="parametric",level=requan*2-1)$upper
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_n_par<-optim(par = coe, fn = mini_nonlinear,method = 'L-BFGS-B', data = set_data)$par
  cf_n<-cf_n_par[1]+cf_n_par[2]*1+cf_n_par[3]*0+cf_n_par[4]*0+cf_n_par[5]*data[480]+cf_n_par[6]*data[477]
  ##list
  list(arima_k_n,arima_p_n,cf_n)
}


##1200
series<-arima_series[201:1400,]

a_1200<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  s1<-rep(c(1,0,0,0),rep)
  s2<-rep(c(0,1,0,0),rep)
  s3<-rep(c(0,0,1,0),rep)
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,s1,s2,s3)
  set_data=cbind(set_data,L1_data,L4_data)
  colnames(set_data)<-c('data','s1','s2','s3','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant=500,mean=0,sd=200)
  arima_k_n<-forecast(arima_k,test_length,interval="parametric",level=requan*2-1)$upper
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_n<-forecast(arima_p,test_length,interval="parametric",level=requan*2-1)$upper
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_n_par<-optim(par = coe, fn = mini_nonlinear,method = 'L-BFGS-B', data = set_data)$par
  cf_n<-cf_n_par[1]+cf_n_par[2]*1+cf_n_par[3]*0+cf_n_par[4]*0+cf_n_par[5]*data[1200]+cf_n_par[6]*data[1197]
  ##list
  list(arima_k_n,arima_p_n,cf_n)
}


##4800
series<-arima_series[201:5000,]

a_4800<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list()),.packages = c("forecast",'quantreg')) %dopar%{
  data<-ts(series[,i],frequency=4)
  length<-length(data)
  test<-ts(data[5:length],frequency=4)
  rep<-length/4
  s1<-rep(c(1,0,0,0),rep)
  s2<-rep(c(0,1,0,0),rep)
  s3<-rep(c(0,0,1,0),rep)
  L1_data=lag(data,-1)
  L4_data=lag(data,-4)
  set_data=cbind(data,s1,s2,s3)
  set_data=cbind(set_data,L1_data,L4_data)
  colnames(set_data)<-c('data','s1','s2','s3','L1','L4')
  set_data=ts(set_data[5:length,])
  #arima
  arima_k<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.3,0.5),constant=500,mean=0,sd=200)
  arima_k_n<-forecast(arima_k,test_length,interval="parametric",level=requan*2-1)$upper
  arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),constant=TRUE)
  arima_p_n<-forecast(arima_p,test_length,interval="parametric",level=requan*2-1)$upper
  #cf
  coe<-lm(data ~., data=set_data)$coefficients
  cf_n_par<-optim(par = coe, fn = mini_nonlinear,method = 'L-BFGS-B', data = set_data)$par
  cf_n<-cf_n_par[1]+cf_n_par[2]*1+cf_n_par[3]*0+cf_n_par[4]*0+cf_n_par[5]*data[4800]+cf_n_par[6]*data[4797]
  ##list
  list(arima_k_n,arima_p_n,cf_n)
}


#generate
##40
series<-arima_series[241,]

re_40<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list())) %dopar%{
  sed<-rnorm(1,30,5)
  arima_k_n<-result_n(series[i],a_40[[1]][[i]],sed)
  arima_p_n<-result_n(series[i],a_40[[2]][[i]],sed)
  cf_n<-result_n(series[i],a_40[[3]][[i]],sed)
  list(arima_k_n,arima_p_n,cf_n)
}


##120
series<-arima_series[321,]

re_120<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list())) %dopar%{
  sed<-rnorm(1,30,5)
  arima_k_n<-result_n(series[i],a_120[[1]][[i]],sed)
  arima_p_n<-result_n(series[i],a_120[[2]][[i]],sed)
  cf_n<-result_n(series[i],a_120[[3]][[i]],sed)
  list(arima_k_n,arima_p_n,cf_n)
}


##480
series<-arima_series[681,]

re_480<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list())) %dopar%{
  sed<-rnorm(1,30,5)
  arima_k_n<-result_n(series[i],a_480[[1]][[i]],sed)
  arima_p_n<-result_n(series[i],a_480[[2]][[i]],sed)
  cf_n<-result_n(series[i],a_480[[3]][[i]],sed)
  list(arima_k_n,arima_p_n,cf_n)
}


##1200
series<-arima_series[1401,]

re_1200<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list())) %dopar%{
  sed<-rnorm(1,30,5)
  arima_k_n<-result_n(series[i],a_1200[[1]][[i]],sed)
  arima_p_n<-result_n(series[i],a_1200[[2]][[i]],sed)
  cf_n<-result_n(series[i],a_1200[[3]][[i]],sed)
  list(arima_k_n,arima_p_n,cf_n)
}


##4800
series<-arima_series[5001,]

re_4800<-foreach(i =1:iter,.combine='comb',.multicombine=TRUE,.init=list(list(),list(),list())) %dopar%{
  sed<-rnorm(1,30,5)
  arima_k_n<-result_n(series[i],a_4800[[1]][[i]],sed)
  arima_p_n<-result_n(series[i],a_4800[[2]][[i]],sed)
  cf_n<-result_n(series[i],a_4800[[3]][[i]],sed)
  list(arima_k_n,arima_p_n,cf_n)
}

save(re_40,re_120,re_480,re_1200,re_4800,file='nonlinear.Rdata')
