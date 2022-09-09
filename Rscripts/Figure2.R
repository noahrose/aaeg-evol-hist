setwd('~/Dropbox (Princeton)/2021 Rose EvolHistAaeg')

mu=4.85e-9
gen=1/15

st<-read.csv('processed_data/slavetrade.csv')

precip<-read.csv('processed_data/humid_period.csv')[,c('BP','precip')]
precip<-rbind(c(0,0),precip,c(25000,0))

plotPair<-function(pair,boots,overplot=F,ylimhi=1.2){
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.combined.final.txt',sep=''))
	res$time<-gen*res[1:nrow(res),2]/mu
	res$rel<-(2*res$lambda_01)/(res$lambda_00+res$lambda_11)
	plot(log10(res$time+3),res$rel,type='s',ylim=c(0,ylimhi),col=NULL,main=pair,ylab='Relative cross-coalescence',xlab='log10 years before present',lwd=0,axes=F,xlim=c(1.5,6))
	axis(1)
	axis(2,at=c(0,0.5,1))
	# rect(log10(160),-1,log10(350),1,col='lightgrey',border='NA')
	# rect(log10(5000),-1,log10(15000),1,col='lightgrey',border='NA')
	
	for(i in 1:boots){
		res<-read.delim(paste('processed_data/',pair,'.boots/boot_',i,'.combined.final.txt',sep=''))
		res$time<-gen*res[1:nrow(res),2]/mu
		res$rel<-(2*res$lambda_01)/(res$lambda_00+res$lambda_11)
		lines(log10(res$time+3),res$rel,type='s',col='grey',lwd=1)
	}
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.combined.final.txt',sep=''))
	res$time<-gen*res[1:nrow(res),2]/mu
	res$rel<-(2*res$lambda_01)/(res$lambda_00+res$lambda_11)
	lines(log10(res$time+3),res$rel,type='s',col='black',lwd=1)
	abline(h=1,lty=2)
	if(overplot) par(new=T)
	return(res)
}

plotIM<-function(pair,boots,overplot=F,ylimhi=1.05){
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.IM.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt',sep=''))
	res$time<-1e-8*gen*res[1:nrow(res),1]/mu
	res$rel<-res$M
	plot(log10(res$time+3),res$rel,type='s',ylim=c(0,ylimhi),col=NULL,main=pair,ylab='Cumulative migration',xlab='log10 years before present',lwd=0,axes=F,xlim=c(1.5,6))
	axis(1)
	axis(2,at=c(0,0.5,1))
	# rect(log10(160),-1,log10(350),1,col='lightgrey',border='NA')
	# rect(log10(5000),-1,log10(15000),1,col='lightgrey',border='NA')
	
	for(i in 1:boots){
		res<-read.delim(paste('processed_data/',pair,'.boots/boot_',i,'.IM.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt',sep=''))
		res$time<-1e-8*gen*res[1:nrow(res),1]/mu
		res$rel<-res$M
		lines(log10(res$time+3),res$rel,type='s',col='grey',lwd=1)
	}
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.IM.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt',sep=''))
	res$time<-1e-8*gen*res[1:nrow(res),1]/mu
	res$rel<-res$M
		lines(log10(res$time+3),res$rel,type='s',col='black',lwd=1)
	abline(h=1,lty=2)
	if(overplot) par(new=T)
	return(res)
}

plotIMunscaled<-function(pair,boots,overplot=F,ylimhi=1.05){
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.IM.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt',sep=''))
	res$time<-1e-8*res[,1]
	res$rel<-res$M
	plot(log10(res$time+1e-6),res$rel,type='s',ylim=c(0,ylimhi),col=NULL,main=pair,ylab='Cumulative migration',xlab='log10 unscaled coalescent time',lwd=0,axes=F,xlim=c(-6,0))
	axis(1)
	axis(2,at=c(0,0.5,1))
	# rect(log10(160),-1,log10(350),1,col='lightgrey',border='NA')
	# rect(log10(5000),-1,log10(15000),1,col='lightgrey',border='NA')
	
	for(i in 1:boots){
		res<-read.delim(paste('processed_data/',pair,'.boots/boot_',i,'.IM.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt',sep=''))
		res$time<-1e-8*res[,1]
		res$rel<-res$M
		lines(log10(res$time+1e-6),res$rel,type='s',col='grey',lwd=1)
	}
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.IM.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt',sep=''))
	res$time<-1e-8*res[,1]
	res$rel<-res$M
	lines(log10(res$time+1e-6),res$rel,type='s',col='black',lwd=1)
	abline(h=1,lty=2)
	if(overplot) par(new=T)
	return(res)
}



plotMig<-function(pair,boots,overplot=F,ylimhi=1e-4){
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.IM.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt',sep=''))
	res$time<-1e-8*gen*res[1:nrow(res),1]/mu
	res$rel<-res$m
	plot(log10(res$time+3),res$rel,type='s',ylim=c(0,ylimhi),col=NULL,main=pair,ylab='Migration rate',xlab='log10 years before present',lwd=0,axes=F,xlim=c(1.5,6))
	axis(1)
	axis(2,at=c(0,0.0001,0.0002))
	# rect(log10(160),-1,log10(350),1,col='lightgrey',border='NA')
	# rect(log10(5000),-1,log10(15000),1,col='lightgrey',border='NA')
	
	for(i in 1:boots){
		res<-read.delim(paste('processed_data/',pair,'.boots/boot_',i,'.IM.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt',sep=''))
		res$time<-1e-8*gen*res[1:nrow(res),1]/mu
		res$rel<-res$m
		lines(log10(res$time+3),res$rel,type='s',col='grey',lwd=1)
	}
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.IM.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt',sep=''))
	res$time<-1e-8*gen*res[1:nrow(res),1]/mu
	res$rel<-res$m
		lines(log10(res$time+3),res$rel,type='s',col='black',lwd=1)
	abline(h=1,lty=2)
	if(overplot) par(new=T)
	return(res)
}


plotST<-function(overplot=T,axis=T){
	plot(log10(st$bp+1),st$number,col=1,lwd=1,xlim=c(1.5,6),ylim=c(0,2.7e6),xlab='',ylab='',xaxt=NULL,yaxt=NULL,axes=F,main='',type='n')
	x=log10(st$bp)
	y=st$number
	y2 <- rep(y, each=2)
	y2 <- y2[-length(y2)]
	x2 <- rep(x, each=2)[-1]
	polygon(x2, y2, border=NA, col=rgb(0.9,0.6,0.9))

	if(axis){
		axis(4,at=c(0,1e6,2e6),labels=c(0,1,2),tck=-0.01, col=rgb(0.9,0.6,0.9),col.ticks=rgb(0.9,0.6,0.9),col.axis=rgb(0.9,0.6,0.9),line=0.5)
		mtext('Embarkations (millions of enslaved people)',side=4,line=1.5,cex=0.5, col=rgb(0.9,0.6,0.9),at=1e6)
	}
	if(overplot) par(new=T)
}

plotPrecip<-function(overplot=T,offset=0,axis=T){
	plot(log10(precip$BP+1),precip$precip,type='n',xlim=c(1.5,6),ylim=c(0,1300),axes=F,xlab='',ylab='')
	polygon(log10(precip$BP+1),precip$precip,col=rgb(157/256,223/256,0),border=NA)
	if(axis){
		axis(4,at=c(100,1000),tck=-0.01, col= rgb(157/256,223/256,0),col.ticks= rgb(157/256,223/256,0),col.axis= rgb(157/256,223/256,0),line=offset)
		mtext('Precipitation in the Sahara (mm/year)',side=4,line=offset+1,cex=0.5, col= rgb(157/256,223/256,0))
	}
	if(overplot) par(new=T)
}

plotNe<-function(pair,boots,overplot=F,ylimhi=7,col1='red',col2='blue'){
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.combined.final.txt',sep=''))
	res$time<-gen*res[1:nrow(res),2]/mu
	res$NE1=(1/res$lambda_00)/(2*mu)
	res$NE2=(1/res$lambda_11)/(2*mu)
	plot(log10(res$time+3),log10(res$NE1),type='s',ylim=c(2,ylimhi),col=NULL,main=pair,ylab='Effective population size',xlab='log10 years before present',lwd=0,axes=F,xlim=c(1.5,6))
	axis(1)
	axis(2)
	# rect(log10(160),-1,log10(350),1,col='lightgrey',border='NA')
	# rect(log10(5000),-1,log10(15000),1,col='lightgrey',border='NA')
	
	for(i in 1:boots){
		res<-read.delim(paste('processed_data/',pair,'.boots/boot_',i,'.combined.final.txt',sep=''))
		res$time<-gen*res[1:nrow(res),2]/mu
		res$NE1=(1/res$lambda_00)/(2*mu)
		res$NE2=(1/res$lambda_11)/(2*mu)
		lines(log10(res$time+3),log10(res$NE1),type='s',col='grey',lwd=1)
		lines(log10(res$time+3),log10(res$NE2),type='s',col='darkgrey',lwd=1)
	}
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.combined.final.txt',sep=''))
	res$time<-gen*res[1:nrow(res),2]/mu
	res$NE1=(1/res$lambda_00)/(2*mu)
	res$NE2=(1/res$lambda_11)/(2*mu)
	lines(log10(res$time+3),log10(res$NE1),type='s',col=col1,lwd=1)
	lines(log10(res$time+3),log10(res$NE2),type='s',col=col2,lwd=1)
	if(overplot) par(new=T)
	return(res)
}


plotNeunscaled<-function(pair,boots,overplot=F,ylimhi=7,col1='red',col2='blue'){
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.combined.final.txt',sep=''))
	res$time<-res[,2]
	res$NE1=(1/res$lambda_00)/(2*mu)
	res$NE2=(1/res$lambda_11)/(2*mu)
	plot(log10(res$time+1e-6),log10(res$NE1),type='s',ylim=c(2,ylimhi),col=NULL,main=pair,ylab='Effective population size',xlab='log10 unscaled coalescent time',lwd=0,axes=F,xlim=c(-6,0))
	axis(1)
	axis(2)
	# rect(log10(160),-1,log10(350),1,col='lightgrey',border='NA')
	# rect(log10(5000),-1,log10(15000),1,col='lightgrey',border='NA')
	
	for(i in 1:boots){
		res<-read.delim(paste('processed_data/',pair,'.boots/boot_',i,'.combined.final.txt',sep=''))
		res$time<-res[,2]
		res$NE1=(1/res$lambda_00)/(2*mu)
		res$NE2=(1/res$lambda_11)/(2*mu)
		lines(log10(res$time+1e-6),log10(res$NE1),type='s',col='grey',lwd=1)
		lines(log10(res$time+1e-6),log10(res$NE2),type='s',col='darkgrey',lwd=1)
	}
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.combined.final.txt',sep=''))
	res$time<-res[,2]
	res$NE1=(1/res$lambda_00)/(2*mu)
	res$NE2=(1/res$lambda_11)/(2*mu)
	lines(log10(res$time+1e-6),log10(res$NE1),type='s',col=col1,lwd=1)
	lines(log10(res$time+1e-6),log10(res$NE2),type='s',col=col2,lwd=1)
	if(overplot) par(new=T)
	return(res)
}



st$rate<-0
st$rate[1:17]<-st$number[1:17]/(st$left[1:17]-st$left[2:18])
strate<-approxfun(st$bp,st$rate,rule=2,method='constant')

getOverlap<-function(pair,currgen=1/10,currmu=7e-9){
	res<-read.delim(paste('processed_data/',pair,'.boots/',pair,'.IM.b1_1e-08.b2_1e-06.MSMC_IM.estimates.txt',sep=''))
	res$time<-1e-8*currgen*res[1:nrow(res),1]/currmu
	currm<-approxfun(res$time,res$m,rule=2,method='constant')
	currmig<-currm(1:5000)/sum(currm(1:5000))
	currst<-strate(1:5000)/sum(strate(1:5000))
	return(sum((currmig*currst)^0.5))
}

overmat<-matrix(nrow=24,ncol=20)
museq=seq(1e-9,1.25e-8,length.out=24)
gseq<-1/(25:5)
ratios<-c()
vals<-c()
for(i in 1:24){
	for(j in 1:20){
		overmat[i,j]<-getOverlap('NGO.SAN',currgen=gseq[j],currmu=museq[i])
		ratios<-c(ratios,gseq[j]/museq[i])
		vals<-c(vals,overmat[i,j])
	}
}
1/(15*ratios[which.max(vals)])
ratios[which.max(vals)]
museq[which.max(apply(overmat,1,max))]
gseq[which.max(apply(overmat,2,max))]

library(plotrix)

NGO.PKT<-plotPair('NGO.PKT',10)





pdf(file='Figures/subplots/Fig2.pdf',width=8,height=6,useDingbats=F)
par(mfrow=c(2,2),mar=c(3,6,3,3),mgp=c(1.5,0.2,0),tck=-0.01,bty='l',xaxs='i',yaxs='i')

par(mar=c(3,5,3,5),mgp=c(1.5,0.2,0),tck=-0.01,bty='l',xaxs='i',yaxs='i')


image(overmat,col=heat.colors(256),x=museq,y=gseq,xlab='Mutation rate (bp/generation)',ylab='Generation time (years)',axes=F)
abline(b=13720930,a=0,lwd=2,lty=2)
gradient.rect(7e-9,0.21,9e-9,0.22,col=heat.colors(256))
par(xpd=T)
text(6.8e-9,0.215,labels='0')
text(9.2e-9,0.215,labels='1')
text(8e-9,0.225,labels='Bhattacharyya coefficient',cex=0.8)
points(2.8e-9,0.036,pch=17,cex=1)
points(1.25e-8,0.036,pch=17,cex=1)
points(4.85e-9,1/15,pch=19,cex=1)
axis(1,at=c(2e-9,6e-9,1e-8),labels=c(2,6,10))
axis(2)
par(xpd=F)
box(bty='l')


plotST(axis=T)
plotMig('NGO.SAN',10,ylimhi=2.5e-4)


plotPrecip(axis=F)
plotNe('NGO.PKT',10,ylimhi=7)
box(bty='l')

plotPrecip(axis=T)
plotIM('NGO.PKT',10)
box(bty='l')

dev.off()

pdf(file='Figures/subplots/NGO.PKT.unscaled.pdf',width=8,height=3,useDingbats=F)
par(mfrow=c(1,2),tck=-0.01,bty='l',xaxs='i',yaxs='i')
par(mar=c(3,5,3,5),mgp=c(1.5,0.2,0),tck=-0.01,bty='l',xaxs='i',yaxs='i')

plotNeunscaled('NGO.PKT',10,ylimhi=7)
box(bty='l')

plotIMunscaled('NGO.PKT',10)
box(bty='l')

dev.off()


pdf(file='Figures/subplots/MSMCgrid.pdf',width=8,height=8,useDingbats=F)
par(mfrow=c(4,4),mar=c(3,3,1,1),mgp=c(1.5,0.2,0),tck=-0.01,bty='l',xaxs='i',yaxs='i',cex.axis=0.7,cex.lab=0.7,oma=c(0,0,0,5))

plotPrecip(axis=F)
plotST(axis=F)
plotNe('NGO.PKT',10,ylimhi=7)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotNe('NGO.OHI',10,ylimhi=7)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotNe('NGO.SAN',10,ylimhi=7)
box(bty='l')

plotPrecip(axis=T,offset=3)
plotST(axis=T)
plotNe('NGO.BKK',10,ylimhi=7)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotPair('NGO.PKT',10,ylimhi= 1.5)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotPair('NGO.OHI',10,ylimhi= 1.5)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotPair('NGO.SAN',10,ylimhi= 1.5)
box(bty='l')

plotPrecip(axis=T,offset=3)
plotST(axis=T)
plotPair('NGO.BKK',10,ylimhi=1.5)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotIM('NGO.PKT',10)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotIM('NGO.OHI',10)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotIM('NGO.SAN',10)
box(bty='l')

plotPrecip(axis=T,offset=3)
plotST(axis=T)
plotIM('NGO.BKK',10)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotMig('NGO.PKT',10,ylimhi=2.5e-4)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotMig('NGO.OHI',10,ylimhi=2.5e-4)
box(bty='l')

plotPrecip(axis=F)
plotST(axis=F)
plotMig('NGO.SAN',10,ylimhi=2.5e-4)
box(bty='l')

plotPrecip(axis=T,offset=3)
plotST(axis=T)
plotMig('NGO.BKK',10,ylimhi=2.5e-4)
box(bty='l')

dev.off()

plotIM('BKK.SAN',10)
