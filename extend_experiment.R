rm(list=ls())
gc()
library('rlist')
library('forecast')
library('smooth')
library('quantreg')

iter<-10000
train_list<-c(40,120,480,1200,4800)
test_length<-1

#set_parameters
series1<-matrix(nrow = 5001, ncol = iter)
for (i in 1:iter) {
  series1[,i]<-arima.sim(list(ar = 0.7), n = 5001,sd=200)+3000
}

series2<-matrix(nrow = 5001, ncol = iter)
for (i in 1:iter) {
  series2[,i]<-arima.sim(list(ar = c(0.7,0.2)), n = 5001,sd=200)+3000
}

series3<-sim.ssarima(frequency=4,orders=list(ar=c(1,1)),lags = c(1,4),AR=c(0.4,0.5),obs = 5001,nsim=iter,constant = 1000,mean=0,sd=200)$data
series4<-sim.ssarima(frequency=7,orders=list(ar=c(2,1)),lags = c(1,7),AR=c(0.4,0.3,0.5),obs = 5001,nsim=iter,constant = 600,mean=0,sd=200)$data

series_list<-list(series1,series2,series3,series4)

cost_par<-matrix(nrow = 4, ncol = 20)
cost_par[,1]<-c(20,10,-3,-7)
cost_par[,2]<-c(20,8,-3,-7)
cost_par[,3]<-c(20,8,3,7)
cost_par[,4]<-c(20,10,3,2)
cost_par[,5]<-c(20,10,1,2)
cost_par[,6]<-c(20,16,1,2)
cost_par[,7]<-c(20,6,-2,2)
cost_par[,8]<-c(20,8,-6,2)
cost_par[,9]<-c(20,7,-2,5)
cost_par[,10]<-c(20,5,-2,5)
cost_par[,11]<-c(20,5,-2,-5)
cost_par[,12]<-c(20,6,-1,5)
cost_par[,13]<-c(20,5,4,5)
cost_par[,14]<-c(20,4,4,6)
cost_par[,15]<-c(20,4,-4,6)
cost_par[,16]<-c(20,4,4,-6)
cost_par[,17]<-c(20,11,4,6)
cost_par[,18]<-c(20,11,2,6)
cost_par[,19]<-c(20,11,4,2)
cost_par[,20]<-c(20,11,1,-8)

quan_par<-c()
for (i in 1:20) {
  quan_par<-c(quan_par,(cost_par[1,i]-cost_par[2,i]+cost_par[4,i])/(cost_par[1,i]+cost_par[3,i]+cost_par[4,i]))
}

#compute
##series_1
ts<-series_list[[1]]
para_combin_1<-list()
for (j in 1:20) {
  cost<-cost_par[,j]
  quan<-quan_par[j]
  profit_function<-function(order,demand){
    profit<-(order>=demand)*((cost[1]+cost[3])*demand-(cost[2]+cost[3])*order)+(order<demand)*((cost[1]-cost[2]+cost[4])*order-cost[4]*demand)
    return(profit)
  }
  mini_function<-function(data,par){
    order<-par[1]+par[2]*data[,2]
    total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function(x[1],x[2])))
    return(-total_profit)
  }
  result_function<-function(series,forecast){
    ppl<-((profit_function(series,series)-profit_function(forecast,series))/profit_function(series,series))
    sl<-forecast>=series
    fr<-min(forecast,series)/series
    list<-c(ppl,sl,fr)
    return(list)
  }
  leng_combin<-list()
  for (i in 1:5){
    leng<-train_list[i]
    series<-ts[201:(200+leng),]
    consult<-ts[(201+leng),]
    arima_for<-list()
    qr_for<-list()
    cf_for<-list()
    for (t in 1:iter) {
      data<-ts(series[,t])
      real_con<-consult[t]
      test<-ts(data[2:leng])
      L1_data=lag(data,-1)
      set_data=cbind(data,L1_data)
      colnames(set_data)<-c('data','L1')
      set_data=ts(set_data[2:leng,])
      #arima
      arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,0)),lags =c(1,4),constant=T)
      arima_p_l<-as.numeric(forecast(arima_p,test_length,interval="parametric",level=quan*2-1)$upper)
      #qr
      qr_l_par<-rq(data ~ .,data=set_data,tau=quan)$coefficients
      qr_l<-as.numeric(qr_l_par[1]+qr_l_par[2]*data[leng])
      #cf
      coe<-lm(data ~., data=set_data)$coefficients
      cf_l_par<-optim(par = coe, fn = mini_function, method = 'L-BFGS-B', data = set_data)$par
      cf_l<-as.numeric(cf_l_par[1]+cf_l_par[2]*data[leng])
      ##list
      arima_for<-list.append(arima_for,result_function(real_con,arima_p_l))
      qr_for<-list.append(qr_for,result_function(real_con,qr_l))
      cf_for<-list.append(cf_for,result_function(real_con,cf_l))
    }
    it_combin<-list(arima_for,qr_for,cf_for)
    leng_combin<-list.append(leng_combin,it_combin)
  }
  para_combin_1<-list.append(para_combin_1,leng_combin)
}

##series_2
ts<-series_list[[2]]
para_combin_2<-list()
for (j in 1:20) {
  cost<-cost_par[,j]
  quan<-quan_par[j]
  profit_function<-function(order,demand){
    profit<-(order>=demand)*((cost[1]+cost[3])*demand-(cost[2]+cost[3])*order)+(order<demand)*((cost[1]-cost[2]+cost[4])*order-cost[4]*demand)
    return(profit)
  }
  mini_function<-function(data,par){
    order<-par[1]+par[2]*data[,2]+par[3]*data[,3]
    total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function(x[1],x[2])))
    return(-total_profit)
  }
  result_function<-function(series,forecast){
    ppl<-((profit_function(series,series)-profit_function(forecast,series))/profit_function(series,series))
    sl<-forecast>=series
    fr<-min(forecast,series)/series
    list<-c(ppl,sl,fr)
    return(list)
  }
  leng_combin<-list()
  for (i in 1:5){
    leng<-train_list[i]
    series<-ts[201:(200+leng),]
    consult<-ts[(201+leng),]
    arima_for<-list()
    qr_for<-list()
    cf_for<-list()
    for (t in 1:iter) {
      data<-ts(series[,t])
      real_con<-consult[t]
      test<-ts(data[3:leng])
      L1_data=lag(data,-1)
      L2_data=lag(data,-2)
      set_data=cbind(data,L1_data,L2_data)
      colnames(set_data)<-c('data','L1','L2')
      set_data=ts(set_data[3:leng,])
      #arima
      arima_p<-ssarima(test,frequency=4,orders=list(ar=c(2,0)),lags =c(1,4),constant=T)
      arima_p_l<-as.numeric(forecast(arima_p,test_length,interval="parametric",level=quan*2-1)$upper)
      #qr
      qr_l_par<-rq(data ~ .,data=set_data,tau=quan)$coefficients
      qr_l<-as.numeric(qr_l_par[1]+qr_l_par[2]*data[leng]+qr_l_par[3]*data[(leng-1)])
      #cf
      coe<-lm(data ~., data=set_data)$coefficients
      cf_l_par<-optim(par = coe, fn = mini_function, method = 'L-BFGS-B', data = set_data)$par
      cf_l<-as.numeric(cf_l_par[1]+cf_l_par[2]*data[leng]+cf_l_par[3]*data[(leng-1)])
      ##list
      arima_for<-list.append(arima_for,result_function(real_con,arima_p_l))
      qr_for<-list.append(qr_for,result_function(real_con,qr_l))
      cf_for<-list.append(cf_for,result_function(real_con,cf_l))
    }
    it_combin<-list(arima_for,qr_for,cf_for)
    leng_combin<-list.append(leng_combin,it_combin)
  }
  para_combin_2<-list.append(para_combin_2,leng_combin)
}

##series_3
ts<-series_list[[3]]
para_combin_3<-list()
for (j in 1:20) {
  cost<-cost_par[,j]
  quan<-quan_par[j]
  profit_function<-function(order,demand){
    profit<-(order>=demand)*((cost[1]+cost[3])*demand-(cost[2]+cost[3])*order)+(order<demand)*((cost[1]-cost[2]+cost[4])*order-cost[4]*demand)
    return(profit)
  }
  mini_function<-function(data,par){
    order<-par[1]+par[2]*data[,2]+par[3]*data[,3]
    total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function(x[1],x[2])))
    return(-total_profit)
  }
  result_function<-function(series,forecast){
    ppl<-((profit_function(series,series)-profit_function(forecast,series))/profit_function(series,series))
    sl<-forecast>=series
    fr<-min(forecast,series)/series
    list<-c(ppl,sl,fr)
    return(list)
  }
  leng_combin<-list()
  for (i in 1:5){
    leng<-train_list[i]
    series<-ts[201:(200+leng),]
    consult<-ts[(201+leng),]
    arima_for<-list()
    qr_for<-list()
    cf_for<-list()
    for (t in 1:iter) {
      data<-ts(series[,t])
      real_con<-consult[t]
      test<-ts(data[5:leng])
      L1_data=lag(data,-1)
      L4_data=lag(data,-4)
      set_data=cbind(data,L1_data,L4_data)
      colnames(set_data)<-c('data','L1','L4')
      set_data=ts(set_data[5:leng,])
      #arima
      arima_p<-ssarima(test,frequency=4,orders=list(ar=c(1,1)),lags =c(1,4),constant=T)
      arima_p_l<-as.numeric(forecast(arima_p,test_length,interval="parametric",level=quan*2-1)$upper)
      #qr
      qr_l_par<-rq(data ~ .,data=set_data,tau=quan)$coefficients
      qr_l<-as.numeric(qr_l_par[1]+qr_l_par[2]*data[leng]+qr_l_par[3]*data[(leng-3)])
      #cf
      coe<-lm(data ~., data=set_data)$coefficients
      cf_l_par<-optim(par = coe, fn = mini_function, method = 'L-BFGS-B', data = set_data)$par
      cf_l<-as.numeric(cf_l_par[1]+cf_l_par[2]*data[leng]+cf_l_par[3]*data[(leng-3)])
      ##list
      arima_for<-list.append(arima_for,result_function(real_con,arima_p_l))
      qr_for<-list.append(qr_for,result_function(real_con,qr_l))
      cf_for<-list.append(cf_for,result_function(real_con,cf_l))
    }
    it_combin<-list(arima_for,qr_for,cf_for)
    leng_combin<-list.append(leng_combin,it_combin)
  }
  para_combin_3<-list.append(para_combin_3,leng_combin)
}

##series_4
ts<-series_list[[4]]
para_combin_4<-list()
for (j in 1:20) {
  cost<-cost_par[,j]
  quan<-quan_par[j]
  profit_function<-function(order,demand){
    profit<-(order>=demand)*((cost[1]+cost[3])*demand-(cost[2]+cost[3])*order)+(order<demand)*((cost[1]-cost[2]+cost[4])*order-cost[4]*demand)
    return(profit)
  }
  mini_function<-function(data,par){
    order<-par[1]+par[2]*data[,2]+par[3]*data[,3]+par[4]*data[,4]
    total_profit<-sum(apply(cbind(order,data[,1]),1,function(x) profit_function(x[1],x[2])))
    return(-total_profit)
  }
  result_function<-function(series,forecast){
    ppl<-((profit_function(series,series)-profit_function(forecast,series))/profit_function(series,series))
    sl<-forecast>=series
    fr<-min(forecast,series)/series
    list<-c(ppl,sl,fr)
    return(list)
  }
  leng_combin<-list()
  for (i in 1:5){
    leng<-train_list[i]
    series<-ts[201:(200+leng),]
    consult<-ts[(201+leng),]
    arima_for<-list()
    qr_for<-list()
    cf_for<-list()
    for (t in 1:iter) {
      data<-ts(series[,t])
      real_con<-consult[t]
      test<-ts(data[8:leng])
      L1_data=lag(data,-1)
      L2_data=lag(data,-2)
      L7_data=lag(data,-7)
      set_data=cbind(data,L1_data,L2_data,L7_data)
      colnames(set_data)<-c('data','L1','L2','L7')
      set_data=ts(set_data[8:leng,])
      #arima
      arima_p<-ssarima(test,frequency=7,orders=list(ar=c(2,1)),lags =c(1,7),constant=T)
      arima_p_l<-as.numeric(forecast(arima_p,test_length,interval="parametric",level=quan*2-1)$upper)
      #qr
      qr_l_par<-rq(data ~ .,data=set_data,tau=quan)$coefficients
      qr_l<-as.numeric(qr_l_par[1]+qr_l_par[2]*data[leng]+qr_l_par[3]*data[(leng-1)]+qr_l_par[4]*data[(leng-6)])
      #cf
      coe<-lm(data ~., data=set_data)$coefficients
      cf_l_par<-optim(par = coe, fn = mini_function, method = 'L-BFGS-B', data = set_data)$par
      cf_l<-as.numeric(cf_l_par[1]+cf_l_par[2]*data[leng]+cf_l_par[3]*data[(leng-1)]+cf_l_par[4]*data[(leng-6)])
      ##list
      arima_for<-list.append(arima_for,result_function(real_con,arima_p_l))
      qr_for<-list.append(qr_for,result_function(real_con,qr_l))
      cf_for<-list.append(cf_for,result_function(real_con,cf_l))
    }
    it_combin<-list(arima_for,qr_for,cf_for)
    leng_combin<-list.append(leng_combin,it_combin)
  }
  para_combin_4<-list.append(para_combin_4,leng_combin)
}

save(para_combin_1,para_combin_2,para_combin_3,para_combin_4,file='extend_result.Rdata')
