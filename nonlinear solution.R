profit_function_nonlinear<-function(order,demand){
  profit<-(order>=demand)*(24*demand-12*order+5*min((order-demand),rnorm(1,30,5)))+(order<demand)*(12*order-0.01*(demand-order)^2)
  return(profit)
}

order<-seq(-600,600,1)
quan_list<-c()
for (k in 1:50){
  list<-c()
  quan<-0
  for (i in order){
    precord<-0
    meanprofit<-0
    for (j in 1:1000){
      de<-rnorm(1,0,200)
      profit<-profit_function_nonlinear(i,de)
      precord<-precord+profit
    }
    meanprofit<-precord/1000
    list<-c(list,meanprofit)
  }
  quan<-which.max(list)/length(order)
  quan_list<-c(quan_list,quan)
}
requan<-mean(quan_list)
