hp<-read.csv('~/Dropbox (Princeton)/2021 Rose EvolHistAaeg/processed_data/humid_period_trace.csv',header=F)
xs<-seq(0,25000,100)
xs<-sort(xs)
GC27<-approxfun(hp$V4[hp$V3=='GC27'],hp$V5[hp$V3=='GC27'],rule=2)(xs)
GC37<-approxfun(hp$V4[hp$V3=='GC37'],hp$V5[hp$V3=='GC37'],rule=2)(xs)
GC49<-approxfun(hp$V4[hp$V3=='GC49'],hp$V5[hp$V3=='GC49'],rule=2)(xs)
GC68<-approxfun(hp$V4[hp$V3=='GC68'],hp$V5[hp$V3=='GC68'],rule=2)(xs)
humid<-data.frame(BP=xs,precip_raw=(GC27+GC37+GC49+GC68)/4)
library(zoo)
humid$precip<-rollmean(humid$precip_raw,k=10,fill='extend')
plot(humid$BP,humid$precip,type='l')
write.csv(humid,file='~/Dropbox (Princeton)/2021 Rose EvolHistAaeg/processed_data/humid_period.csv',quote=F,row.names=F)