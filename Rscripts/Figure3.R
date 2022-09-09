setwd('~/Dropbox (Princeton)/2021 Rose EvolHistAaeg')

plotVit<-function(str,h=1,newPlot=T,doxlab=T,main=''){
	post<-read.delim(paste('processed_data/',str,sep=''),header=F)
	if(newPlot){
		currx=''
		plot(1,type='n',xlim=c(0,max(post$V3)),xlab=currx,ylab='',axes=F,main=main,font.main=1)
		if(doxlab){
			axis(1,at=c(0,1e8,2e8,3e8),labels=c(0,100,200,300))
			currx='Chr. 1 position (Mb)'
		}
	}
	segments(0,h,max(post$V3),h,lwd=5,col='blue',lend=1)
	# segments(245e6,h,308e6,h,lwd=5,col='grey',lend=1)
	# segments(145e6,h,155e6,h,lwd=5,col='grey',lend=1)
	if(length(which(post[,6]=='1,1'))>0) segments(post[post[,6]=='1,1',2],h,post[post[,6]=='1,1',3],h,col='gold',lwd=5,lend=1)
	if(length(which(post[,6]=='0,2'))>0) segments(post[post[,6]=='0,2',2],h,post[post[,6]=='0,2',3],h,col='orange',lwd=5,lend=1)
}

pdf(file='Figures/subplots/f3.pdf',width=1,height=3.25,useDingbats=F)
par(mfrow=c(2,1),mar=c(2,2,0,0),mgp=c(1.1,0.2,0),tck=-0.01,bty='l')

f3s<-read.delim('processed_data/threepop_verily_norepeat_maf05.txt',header=F,sep=' ')
f3s<-f3s[match(unique(f3s[,1]),f3s[,1]),]
rownames(f3s)<-f3s[,1]

trios<-c('KUM;SAN,OHI','KUM;OHI,BKK','KUM;NGO,OHI','KUM;KIN,SAN','KUM;KIN,BKK','KUM;NGO,KIN')
plot(f3s[trios,2],pch=19,ylim=c(-0.005,0.001),bty='l',axes=F,xlab='',ylab='f3 statistic',cex=0.5,cex.lab=0.5)
arrows(1:6,f3s[trios,2]-1.96*f3s[trios,3],1:6,f3s[trios,2]+1.96*f3s[trios,3],code=3,angle=90,length=0.01)
axis(1,at=1:6,labels=c('SAN','BKK','NGO','SAN','BKK','NGO'),cex.axis=0.7,las=2)
mtext('KUM;OHI,X',side=1,line=1,at=2,cex=0.4)
mtext('KUM;KIN,X',side=1,line=1,at=5,cex=0.4)
axis(2,at=c(-0.005,0),cex.axis=0.5)
abline(h=0,lty=2)
box(bty='l')

trios<-c('OGD;SAN,OHI','OGD;OHI,BKK','OGD;NGO,OHI','OGD;KIN,SAN','OGD;KIN,BKK','OGD;NGO,KIN')
plot(f3s[trios,2],pch=19,ylim=c(-0.005,0.001),bty='l',axes=F,xlab='',ylab='f3 statistic',cex=0.5,cex.lab=0.5)
arrows(1:6,f3s[trios,2]-1.96*f3s[trios,3],1:6,f3s[trios,2]+1.96*f3s[trios,3],code=3,angle=90,length=0.01)
axis(1,at=1:6,labels=c('SAN','BKK','NGO','SAN','BKK','NGO'),cex.axis=0.7,las=2)
mtext('OGD;OHI,X',side=1,line=1,at=2,cex=0.4)
mtext('OGD;KIN,X',side=1,line=1,at=5,cex=0.4)
axis(2,at=c(-0.005,0),cex.axis=0.5)
abline(h=0,lty=2)

box(bty='l')
dev.off()


pdf(file='Figures/subplots/vits.pdf',width=6,height=6)
par(mfrow=c(2,3),mar=c(3,3,1,1),mgp=c(1.2,0.2,0),tck=-0.01,bty='l')


samps<-list.files('processed_data/kum_posterior/chr1/')

plot(1,xlim=c(0,475e6),ylim=c(0,20),type='n',bty='n',axes=F,xlab='',ylab='')
for(i in 1:length(samps)){
	plotVit(paste('kum_posterior/chr1/',samps[i],sep=''),newPlot=F,h=i)
}

plot(1,xlim=c(0,475e6),ylim=c(0,20),type='n',bty='n',axes=F,xlab='',ylab='')
for(i in 1:length(samps)){
	plotVit(paste('kum_posterior/chr2/',samps[i],sep=''),newPlot=F,h=i)
}

plot(1,xlim=c(0,475e6),ylim=c(0,20),type='n',bty='n',axes=F,xlab='',ylab='')
for(i in 1:length(samps)){
	plotVit(paste('kum_posterior/chr3/',samps[i],sep=''),newPlot=F,h=i)
}

samps<-list.files('processed_data/ogd_posterior/chr1/')

plot(1,xlim=c(0,475e6),ylim=c(0,20),type='n',bty='n',axes=F,xlab='',ylab='')
for(i in 1:length(samps)){
	plotVit(paste('ogd_posterior/chr1/',samps[i],sep=''),newPlot=F,h=i)
}

plot(1,xlim=c(0,475e6),ylim=c(0,20),type='n',bty='n',axes=F,xlab='',ylab='')
for(i in 1:length(samps)){
	plotVit(paste('ogd_posterior/chr2/',samps[i],sep=''),newPlot=F,h=i)
}

plot(1,xlim=c(0,475e6),ylim=c(0,20),type='n',bty='n',axes=F,xlab='',ylab='')
for(i in 1:length(samps)){
	plotVit(paste('ogd_posterior/chr3/',samps[i],sep=''),newPlot=F,h=i)
}
dev.off()

pdf(file='Figures/subplots/urban.pdf',width=3.25,height=3.25)
par(mar=c(2,2,0,0),mgp=c(1.2,0.2,0),tck=-0.01,bty='l')

popsize<-read.csv('processed_data/urban_pops.csv')

plot(popsize$Year, popsize$Kumasi,type='l',lwd=3,col='blue',xlab='Year',ylab='Human population')
lines(popsize$Year, popsize$Ouagadougou,type='l',lwd=3,col='darkblue')
abline(v=2020,lty=2,lwd=2)

c1<-read.delim('processed_data/kum_posterior/chr1_estimates.txt',header=F)
conf<-quantile(c1[2:80,1],c(0.05,0.95))/15
arrows(2018-conf[1],3e6,2018-conf[2],3e6,lwd=2,col='blue',code=3,angle=90,length= 0.05)
text(2018-conf[2]-10,3e6,labels='Chr. 1',col='blue')

c2<-read.delim('processed_data/kum_posterior/chr2_estimates.txt',header=F)
conf<-quantile(c2[2:80,1],c(0.05,0.95))/15
arrows(2018-conf[1],2.7e6,2018-conf[2],2.7e6,lwd=2,col='blue',code=3,angle=90,length= 0.05)
text(2018-conf[2]-10,2.7e6,labels='Chr. 2',col='blue')

c3<-read.delim('processed_data/kum_posterior/chr3_estimates.txt',header=F)
conf<-quantile(c3[2:80,1],c(0.05,0.95))/15
arrows(2018-conf[1],2.4e6,2018-conf[2],2.4e6,lwd=2,col='blue',code=3,angle=90,length= 0.05)
text(2018-conf[2]-10,2.4e6,labels='Chr. 3',col='blue')

kumest<-mean(c(c1[1,1],c2[1,1],c3[1,1]))/15
arrows(2018-kumest,3e6,2018-kumest,1.5e6,lwd=2,col='blue',length=0.1)
abline(v=2020,lty=2,lwd=2)


c1<-read.delim('processed_data/ogd_posterior/chr1_estimates.txt',header=F)
conf<-quantile(c1[2:80,1],c(0.05,0.95))/15
arrows(2018-conf[1],4e6,2018-conf[2],4e6,lwd=2,col='darkblue',code=3,angle=90,length= 0.05)
text(2018-conf[2]-10,4e6,labels='Chr. 1',col='darkblue')

c2<-read.delim('processed_data/ogd_posterior/chr2_estimates.txt',header=F)
conf<-quantile(c2[2:80,1],c(0.05,0.95))/15
arrows(2018-conf[1],3.7e6,2018-conf[2],3.7e6,lwd=2,col='darkblue',code=3,angle=90,length= 0.05)
text(2018-conf[2]-10,3.7e6,labels='Chr. 2',col='darkblue')

c3<-read.delim('processed_data/ogd_posterior/chr3_estimates.txt',header=F)
conf<-quantile(c3[2:80,1],c(0.05,0.95))/15
arrows(2018-conf[1],3.4e6,2018-conf[2],3.4e6,lwd=2,col='darkblue',code=3,angle=90,length= 0.05)
text(2018-conf[2]-10,3.4e6,labels='Chr. 3',col='darkblue')

ogdest<-mean(c(c1[1,1],c2[1,1],c3[1,1]))/15
arrows(2018-ogdest,4e6,2018-ogdest,1.1e6,lwd=2,col='darkblue',length=0.1)

text(1995,4.7e6,labels='Human-specialist\nadmixture pulse',font=2,col='orange',cex=0.8)
text(1955,5e6,labels='KUM',font=2,col='blue')
text(1955,4.7e6,labels='OGD',font=2,col='darkblue')

dev.off()


getVit<-function(str,vars){
	post<-read.delim(paste('processed_data/',str,sep=''),header=F)
	post$gt<-as.numeric(factor(post[,6],levels=c('2,0','1,1','0,2')))-1
	post[nrow(post),3]<-post[nrow(post),3]-1
	genos<-rep(0,length(vars))
	for(i in 1:nrow(post)){
		genos[which(vars==post[i,2]):which(vars==(post[i,3]+1))]<-post$gt[i]
	}
	return(genos)
}

getHet<-function(str){
	post<-read.delim(paste('processed_data/',str,sep=''),header=F)
	return(post[post[,6]=='1,1',3]-post[post[,6]=='1,1',2])
}


samps<-list.files('processed_data/kum_posterior/chr1/')
k1<-read.delim('processed_data/kum_posterior/kum_chr1.gmap.inputfile',header=F)
k2<-read.delim('processed_data/kum_posterior/kum_chr2.gmap.inputfile',header=F)
k3<-read.delim('processed_data/kum_posterior/kum_chr3.gmap.inputfile',header=F)

v1<-c(0,k1[,2])
v2<-c(0,k2[,2])
v3<-c(0,k3[,2])

kc1<-sapply(paste0('kum_posterior/chr1/',samps),getVit,vars=v1)
kc2<-sapply(paste0('kum_posterior/chr2/',samps),getVit,vars=v2)
kc3<-sapply(paste0('kum_posterior/chr3/',samps),getVit,vars=v3)
kc1aaa<-rowMeans(kc1)/2
kc2aaa<-rowMeans(kc2)/2
kc3aaa<-rowMeans(kc3)/2
kumHet<-c(unlist(sapply(paste0('kum_posterior/chr1/',samps),getHet)),unlist(sapply(paste0('kum_posterior/chr2/',samps),getHet)),unlist(sapply(paste0('kum_posterior/chr3/',samps),getHet)))

samps<-list.files('processed_data/ogd_posterior/chr1/')
o1<-read.delim('processed_data/ogd_posterior/ogd_chr1.gmap.inputfile',header=F)
o2<-read.delim('processed_data/ogd_posterior/ogd_chr2.gmap.inputfile',header=F)
o3<-read.delim('processed_data/ogd_posterior/ogd_chr3.gmap.inputfile',header=F)
oc1<-sapply(paste0('ogd_posterior/chr1/',samps),getVit,vars=v1)
oc2<-sapply(paste0('ogd_posterior/chr2/',samps),getVit,vars=v2)
oc3<-sapply(paste0('ogd_posterior/chr3/',samps),getVit,vars=v3)
oc1aaa<-rowMeans(oc1)/2
oc2aaa<-rowMeans(oc2)/2
oc3aaa<-rowMeans(oc3)/2
ogdHet<-c(unlist(sapply(paste0('ogd_posterior/chr1/',samps),getHet)),unlist(sapply(paste0('ogd_posterior/chr2/',samps),getHet)),unlist(sapply(paste0('ogd_posterior/chr3/',samps),getHet)))

samps1<-list.files('processed_data/ancestry_sims/',pattern='sim.chr1.2000000.*')
samps2<-list.files('processed_data/ancestry_sims/',pattern='sim.chr2.2000000.*')
samps3<-list.files('processed_data/ancestry_sims/',pattern='sim.chr3.2000000.*')
sc1<-sapply(paste0('ancestry_sims/',samps1),getVit,vars=v1)
sc2<-sapply(paste0('ancestry_sims/',samps2),getVit,vars=v2)
sc3<-sapply(paste0('ancestry_sims/',samps3),getVit,vars=v3)
sc1aaa<-rowMeans(sc1)/2
sc2aaa<-rowMeans(sc2)/2
sc3aaa<-rowMeans(sc3)/2

samps1<-list.files('processed_data/ancestry_sims/',pattern='sim2.chr1.2000000.*')
samps2<-list.files('processed_data/ancestry_sims/',pattern='sim2.chr2.2000000.*')
samps3<-list.files('processed_data/ancestry_sims/',pattern='sim2.chr3.2000000.*')
s2c1<-sapply(paste0('ancestry_sims/',samps1),getVit,vars=v1)
s2c2<-sapply(paste0('ancestry_sims/',samps2),getVit,vars=v2)
s2c3<-sapply(paste0('ancestry_sims/',samps3),getVit,vars=v3)
s2c1aaa<-rowMeans(s2c1)/2
s2c2aaa<-rowMeans(s2c2)/2
s2c3aaa<-rowMeans(s2c3)/2



# s2mb<-read.delim('processed_data/ancestry_sims/sim.tracts.2000000.160.1mb.bed',header=F)
# s2mb2<-read.delim('processed_data/ancestry_sims/sim2.tracts.2000000.160.1mb.bed',header=F)
outl<-read.delim('processed_data/outlier_regions.bed',header=F)



#############################
#shared signal between cities
#############################


truev<-cor(c(kc1aaa,kc2aaa,kc3aaa),c(oc1aaa,oc2aaa,oc3aaa))
spin<-function(x) x[1+(sample(1:length(x),1)+1:length(x))%%length(x)]
perms<-replicate(1000,cor(c(spin(kc1aaa),spin(kc2aaa),spin(kc3aaa)),c(oc1aaa,oc2aaa,oc3aaa)))
hist(perms,breaks=100,col='grey',xlim=c(-1,1))
abline(v=truev)
truev
mean(perms>truev)
#cor 0.50
#one tailed p<0.001 



pdf(file='Figures/subplots/tract_distribution.pdf',width=6,height=6,useDingbats=F)

par(mfrow=c(5,3),mar=c(2.5,2.5,0,1),mgp=c(1.2,0.2,0),tck=-0.02,bty='l',cex.axis=0.7,cex.lab=0.7)

plot(c(0,k1[,2]),kc1aaa,lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='KUM prop. specialist',type='n')
axis(2)
rect(outl[outl$V1==1,2],0,outl[outl$V1==1,3],1,col='lightgrey',border=NA)
lines(c(0,k1[,2]),kc1aaa,lwd=1)

plot(c(0,k2[,2]),kc2aaa,type='n',lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='')
rect(outl[outl$V1==2,2],0,outl[outl$V1==2,3],1,col='lightgrey',border=NA)
lines(c(0,k2[,2]),kc2aaa,lwd=1)

plot(c(0,k3[,2]),kc3aaa,type='n',lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='')
rect(outl[outl$V1==3,2],0,outl[outl$V1==3,3],1,col='lightgrey',border=NA)
lines(c(0,k3[,2]),kc3aaa,lwd=1)


plot(c(0,o1[,2]),oc1aaa,lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='OGD prop. specialist',type='n')
axis(2)
rect(outl[outl$V1==1,2],0,outl[outl$V1==1,3],1,col='lightgrey',border=NA)
lines(c(0,o1[,2]),oc1aaa,lwd=1)

plot(c(0,o2[,2]),oc2aaa,type='n',lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='')
rect(outl[outl$V1==2,2],0,outl[outl$V1==2,3],1,col='lightgrey',border=NA)
lines(c(0,o2[,2]),oc2aaa,lwd=1)

plot(c(0,o3[,2]),oc3aaa,type='n',lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='')
rect(outl[outl$V1==3,2],0,outl[outl$V1==3,3],1,col='lightgrey',border=NA)
lines(c(0,o3[,2]),oc3aaa,lwd=1)

plot(c(0,o1[,2]),oc1aaa,lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='KUM + OGD prop. specialist',type='n')
axis(2)
rect(outl[outl$V1==1,2],0,outl[outl$V1==1,3],1,col='lightgrey',border=NA)
lines(c(0,o1[,2]),(oc1aaa+kc1aaa)/2,lwd=1)

plot(c(0,o2[,2]),oc2aaa,type='n',lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='')
rect(outl[outl$V1==2,2],0,outl[outl$V1==2,3],1,col='lightgrey',border=NA)
lines(c(0,o2[,2]),(oc2aaa+kc2aaa)/2,lwd=1)

plot(c(0,o3[,2]),oc3aaa,type='n',lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='')
rect(outl[outl$V1==3,2],0,outl[outl$V1==3,3],1,col='lightgrey',border=NA)
lines(c(0,o3[,2]),(oc3aaa+kc3aaa)/2,lwd=1)

plot(c(0,k1[,2]),sc1aaa,lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='Sim1 prop. specialist',type='n')
axis(2)
rect(outl[outl$V1==1,2],0,outl[outl$V1==1,3],1,col='lightgrey',border=NA)
lines(c(0,k1[,2]),sc1aaa,lwd=1)

plot(c(0,k2[,2]),sc2aaa,type='n',lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='')
rect(outl[outl$V1==2,2],0,outl[outl$V1==2,3],1,col='lightgrey',border=NA)
lines(c(0,k2[,2]),sc2aaa,lwd=1)

plot(c(0,k3[,2]),sc3aaa,type='n',lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='')
rect(outl[outl$V1==3,2],0,outl[outl$V1==3,3],1,col='lightgrey',border=NA)
lines(c(0,k3[,2]),sc3aaa,lwd=1)


plot(c(0,k1[,2]),s2c1aaa,lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='Sim1 prop. specialist',type='n')
axis(2)
rect(outl[outl$V1==1,2],0,outl[outl$V1==1,3],1,col='lightgrey',border=NA)
lines(c(0,k1[,2]),s2c1aaa,lwd=1)

plot(c(0,k2[,2]),s2c2aaa,type='n',lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='')
rect(outl[outl$V1==2,2],0,outl[outl$V1==2,3],1,col='lightgrey',border=NA)
lines(c(0,k2[,2]),s2c2aaa,lwd=1)

plot(c(0,k3[,2]),s2c3aaa,type='n',lwd=1,ylim=c(0,0.7),axes=F,xlab='',ylab='')
rect(outl[outl$V1==3,2],0,outl[outl$V1==3,3],1,col='lightgrey',border=NA)
lines(c(0,k3[,2]),s2c3aaa,lwd=1)


# # plot(s2mb[s2mb $V1=='NC_035107.1','V2'], s2mb[s2mb $V1=='NC_035107.1','V4'],type='n',lwd=2,axes=F,xlab='',ylab='Sim. 1 tract coverage',ylim=c(0,100))
# axis(2,at=c(0,50,100))
# rect(outl[outl$V1==1,2],0,outl[outl$V1==1,3],100,col='lightgrey',border=NA)
# lines(s2mb[s2mb $V1=='NC_035107.1','V2'], s2mb[s2mb $V1=='NC_035107.1','V4'])

# plot(s2mb[s2mb $V1=='NC_035108.1','V2'], s2mb[s2mb $V1=='NC_035108.1','V4'],type='n',lwd=2,axes=F,xlab='',ylab='',ylim=c(0,100))
# rect(outl[outl$V1==2,2],0,outl[outl$V1==2,3],100,col='lightgrey',border=NA)
# lines(s2mb[s2mb $V1=='NC_035108.1','V2'], s2mb[s2mb $V1=='NC_035108.1','V4'])

# plot(s2mb[s2mb $V1=='NC_035109.1','V2'], s2mb[s2mb $V1=='NC_035109.1','V4'],type='n',lwd=2,axes=F,xlab='',ylab='',ylim=c(0,100))
# rect(outl[outl$V1==3,2],0,outl[outl$V1==3,3],100,col='lightgrey',border=NA)
# lines(s2mb[s2mb $V1=='NC_035109.1','V2'], s2mb[s2mb $V1=='NC_035109.1','V4'])


# plot(s2mb2[s2mb2 $V1=='NC_035107.1','V2'], s2mb2[s2mb2 $V1=='NC_035107.1','V4'],type='n',lwd=2,axes=F,xlab='Chr. 1 position (Mb)',ylab='Sim. 2 tract coverage',ylim=c(0,100))
# axis(1,at=c(0,1e8,2e8,3e8),labels=c(0,100,200,300))
# axis(2,at=c(0,50,100))
# rect(outl[outl$V1==1,2],0,outl[outl$V1==1,3],100,col='lightgrey',border=NA)
# lines(s2mb2[s2mb2 $V1=='NC_035107.1','V2'], s2mb2[s2mb2 $V1=='NC_035107.1','V4'])

# plot(s2mb2[s2mb2 $V1=='NC_035108.1','V2'], s2mb2[s2mb2 $V1=='NC_035108.1','V4'],type='n',lwd=2,axes=F,xlab='Chr. 2 position (Mb)',ylab='',ylim=c(0,100))
# axis(1,at=c(0,1e8,2e8,3e8,4e8),labels=c(0,100,200,300,400))
# rect(outl[outl$V1==2,2],0,outl[outl$V1==2,3],100,col='lightgrey',border=NA)
# lines(s2mb2[s2mb2 $V1=='NC_035108.1','V2'], s2mb2[s2mb2 $V1=='NC_035108.1','V4'])

# plot(s2mb2[s2mb2 $V1=='NC_035109.1','V2'], s2mb2[s2mb2 $V1=='NC_035109.1','V4'],type='n',lwd=2,axes=F,xlab='Chr. 3 position (Mb)',ylab='',ylim=c(0,100))
# axis(1,at=c(0,1e8,2e8,3e8,4e8),labels=c(0,100,200,300,400))
# rect(outl[outl$V1==3,2],0,outl[outl$V1==3,3],100,col='lightgrey',border=NA)
# lines(s2mb2[s2mb2 $V1=='NC_035109.1','V2'], s2mb2[s2mb2 $V1=='NC_035109.1','V4'])


dev.off()







#############
#Simulations#
#############



plotSim<-function(str,h=2,newPlot=T,doxlab=T,main=''){
	currTrue<-read.delim(paste('processed_data/ancestry_sims/shuf.',str,'.bed',sep=''),header=F)

	post<-read.delim(paste('processed_data/ancestry_sims/sim.chr1.',str,'.viterbi',sep=''),header=F)
	plot(1,type='n',xlim=c(0,311e6),ylim=c(1,4),xlab='Chr. 1 position (Mb)',ylab='',axes=F,main='',font.main=1)
	axis(1,at=c(0,1e8,2e8,3e8),labels=c(0,100,200,300))
	axis(2,at=c(2,3),labels=c('Sim','AHMM'),las=1)
	segments(0,h,311e6,h,lwd=5,col='blue',lend=1)
	segments(0,h+1,311e6,h+1,lwd=5,col='blue',lend=1)
	segments(post[post[,6]=='1,1',2],h+1,post[post[,6]=='1,1',3],h+1,col='gold',lwd=5,lend=1)
	segments(currTrue[currTrue$V1=='NC_035107.1',2],h,currTrue[currTrue$V1=='NC_035107.1',3],h,col='gold',lwd=5,lend=1)

	post<-read.delim(paste('processed_data/ancestry_sims/sim.chr2.',str,'.viterbi',sep=''),header=F)
	plot(1,type='n',xlim=c(0,474e6),ylim=c(1,4),xlab='Chr. 2 position (Mb)',ylab='',axes=F,main=main,font.main=1)
	axis(1,at=c(0,1e8,2e8,3e8,4e8),labels=c(0,100,200,300,400))
	segments(0,h,474e6,h,lwd=5,col='blue',lend=1)
	segments(0,h+1,474e6,h+1,lwd=5,col='blue',lend=1)
	segments(post[post[,6]=='1,1',2],h+1,post[post[,6]=='1,1',3],h+1,col='gold',lwd=5,lend=1)
	segments(currTrue[currTrue$V1=='NC_035108.1',2],h,currTrue[currTrue$V1=='NC_035108.1',3],h,col='gold',lwd=5,lend=1)

	post<-read.delim(paste('processed_data/ancestry_sims/sim.chr3.',str,'.viterbi',sep=''),header=F)
	plot(1,type='n',xlim=c(0,410e6),ylim=c(1,4),xlab='Chr. 3 position (Mb)',ylab='',axes=F,main='',font.main=1)
	axis(1,at=c(0,1e8,2e8,3e8,4e8),labels=c(0,100,200,300,400))
	segments(0,h,410e6,h,lwd=5,col='blue',lend=1)
	segments(0,h+1,410e6,h+1,lwd=5,col='blue',lend=1)
	segments(post[post[,6]=='1,1',2],h+1,post[post[,6]=='1,1',3],h+1,col='gold',lwd=5,lend=1)
	segments(currTrue[currTrue$V1=='NC_035109.1',2],h,currTrue[currTrue$V1=='NC_035109.1',3],h,col='gold',lwd=5,lend=1)
}





jacc<-read.delim('processed_data/ancestry_sims/all_jaccards.txt',header=F,strings=F,sep=' ')
js<-read.delim('processed_data/ancestry_sims/all_jacshuff.txt',header=F,strings=F,sep=' ')

getTrues<-function(i,str){
	currTrue<-read.delim(paste('processed_data/ancestry_sims/shuf.',str,'.',i,'.bed',sep=''),header=F)
	trues<-currTrue[,3]-currTrue[,2]
}

getPosts<-function(i,str){
	post<-read.delim(paste('processed_data/ancestry_sims/sim.chr1.',str,'.',i,'.viterbi',sep=''),header=F)
	posts<-(post[post[,6]=='1,1',3]-post[post[,6]=='1,1',2])
	post<-read.delim(paste('processed_data/ancestry_sims/sim.chr2.',str,'.',i,'.viterbi',sep=''),header=F)
	posts<-c(posts,(post[post[,6]=='1,1',3]-post[post[,6]=='1,1',2]))
	post<-read.delim(paste('processed_data/ancestry_sims/sim.chr3.',str,'.',i,'.viterbi',sep=''),header=F)
	posts<-c(posts,(post[post[,6]=='1,1',3]-post[post[,6]=='1,1',2]))
}


t1<-unlist(sapply(1:100,getTrues,'10000000.32'))
p1<-unlist(sapply(1:100,getPosts,'10000000.32'))
t4<-unlist(sapply(1:100,getTrues,'2000000.160'))
p4<-unlist(sapply(1:100,getPosts,'2000000.160'))
t2<-unlist(sapply(1:100,getTrues,'1000000.320'))
p2<-unlist(sapply(1:100,getPosts,'1000000.320'))
t3<-unlist(sapply(1:100,getTrues,'500000.640'))
p3<-unlist(sapply(1:100,getPosts,'500000.640'))

t1d<-cbind(rep(log10(10e6),length(t1)),log10(t1))
t2d<-cbind(rep(log10(1e6),length(t2)),log10(t2))
t3d<-cbind(rep(log10(5e5),length(t3)),log10(t3))
t4d<-cbind(rep(log10(2e6),length(t4)),log10(t4))
p1d<-cbind(rep(log10(10e6),length(p1)),log10(p1))
p2d<-cbind(rep(log10(1e6),length(p2)),log10(p2))
p3d<-cbind(rep(log10(5e5),length(p3)),log10(p3))
p4d<-cbind(rep(log10(2e6),length(p4)),log10(p4))

ts<-rbind(t1d,t4d,t2d,t3d)
ps<-rbind(p1d,p4d,p2d,p3d)

agg1<-aggregate(ts[,2]~ts[,1],FUN=function(x) c(median(x),sd(x)/sqrt(length(x))))
agg2<-aggregate(ps[,2]~ps[,1],FUN=function(x) c(median(x),sd(x)/sqrt(length(x))))

library(vioplot)

pdf(file='Figures/subplots/ancestry_sim.pdf',width=6,height=6,useDingbats=F)
par(mar=c(3,3,1,1),mgp=c(1.2,0.2,0),tck=-0.02,bty='l',cex.lab=0.7,cex.axis=0.7)
layout(mat=cbind(rbind(1:3,1:3,1:3,4:6,4:6,4:6,7:9,7:9,7:9,10:12,10:12,10:12),c(rep(13,4),rep(14,4),rep(15,4)),c(rep(13,4),rep(14,4),rep(15,4))))
plotSim('10000000.32.1',doxlab=F,main='10 Mb')
plotSim('2000000.160.1',doxlab=F,main='2 Mb')
plotSim('1000000.320.1',doxlab=F,main='1 Mb')
plotSim('500000.640.1',doxlab=F,main='500 kb')

plot(jitter(log10(c(rep(10e6,100),rep(1e6,100),rep(2e6,100),rep(5e5,100))),fact=1),jacc$V3,xlab='Tract length (log10 scale)',ylab='Jaccard statistic',cex=0.5,col=rgb(1,0,0,0.5),ylim=c(0,1),axes=F)
points(jitter(log10(c(rep(10e6,100),rep(1e6,100),rep(2e6,100),rep(5e5,100))),fact=1),js$V3,xlab='Tract length (log10 scale)',ylab='Jaccard statistic',cex=0.5,col=grey(0,0.5))
axis(2)
axis(1,at=log10(c(5e5,1e6,2e6,10e6)),labels=c('500kb','1Mb','2Mb','10Mb'),las=1)
box(bty='l')
text(log10(7e5),0.9,labels='permuted')
text(log10(7e5),1,labels='AHMM',col='red')

vioplot(list(t3d[,2],p3d[,2],t2d[,2],p2d[,2],t4d[,2],p4d[,2],t1d[,2],p1d[,2]),ylim=c(4,8),col=c('grey','red','grey','red','grey','red','grey','red'),colMed='black',range=0,ylab='Log10 tract length',bty='l')
text(1.5,7.5,labels='simulated')
text(1.5,8,labels='AHMM',col='red')


plot(agg1[,2][,1],agg2[,2][,1],xlim=c(5.5,7.1),ylim=c(5.5,7.1),type='b',xlab='Median simulated tract length',ylab='Median AHMM tract length')
abline(0,1,lty=2)
text(5.75,7,labels=expression(R^2==0.99))
dev.off()

summary(lm(agg1[,2][,1]~agg2[,2][,1]))



