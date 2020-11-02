x = scan("F:/DataSUC/code101520/WindDA1month.txt")
Day = 30

#scenarios for m = 1
# set.seed(1234)
# m = 1 * 3 * 24
# write.table(matrix(pmax(rnorm(3*24*50*Day,4*x[(2161-m):(4320-m)],4*0.05*x[(2161-m):(4320-m)]),0),3,50*24*Day),sep=",",col.names=FALSE,row.names=FALSE,
#             file="F:/DataSUC/code11.15/code11.15/118buswind5e.csv")

#scenarios for m = 1
set.seed(1234)
m = 10
f = rep(0,3*24)
f[1] = 1/m
fm = rep(f,m)
x_pred = filter(x[(2161-m*3*24):(4320-1)],rev(fm),sides = 1)
x_pred = na.omit(x_pred)
length(x_pred)
write.table(matrix(pmax(rnorm(3*24*50*Day,4*x_pred,4*0.05*x_pred),0),3,50*24*Day),sep=",",col.names=FALSE,row.names=FALSE,
            file="F:/DataSUC/code101520/118buswind5em5.csv")

#ture value
set.seed(12345)
write.table(matrix(pmax(rnorm(3*24*1000*Day,4*x_pred,4*0.05*x_pred),0),3,1000*24*Day),sep=",",row.names = F,
            col.names = F,file="F:/DataSUC/code101520/118buswindT5e.csv")

#deterministic
set.seed(123)
write.table(matrix(x,24,1),sep=",",row.names = F,
            col.names = F,file="C:/Users/wang.keq/Desktop/code (2)/code11.15/sixbuswindD.csv")

x_pred[length(x_pred)]
a <- c(1,2,3,4,5,6)
filter(a, c(0,0,1),sides = 1)
