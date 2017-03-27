args <- commandArgs(trailingOnly = TRUE)

x<-unlist(strsplit(args[2],split=","))
x<-as.numeric(x)
args[3]
paril<-1-pnorm(as.numeric(args[1]),mean(x),sd(x))
poutlier<-pbinom(0,length(x),paril)
1-poutlier

x<-unlist(strsplit(args[5],split=","))
x<-as.numeric(x)
args[6]
#1-exp(log(pnorm(as.numeric(args[4]),mean(x),sd(x)))*length(x))
paril<-1-pnorm(as.numeric(args[4]),mean(x),sd(x))
poutlier<-pbinom(0,length(x),paril)
1-poutlier

x<-unlist(strsplit(args[8],split=","))
x<-as.numeric(x)
args[9]
#1-exp(log(pnorm(as.numeric(args[7]),mean(x),sd(x)))*length(x))
paril<-1-pnorm(as.numeric(args[7]),mean(x),sd(x))
poutlier<-pbinom(0,length(x),paril)
1-poutlier
