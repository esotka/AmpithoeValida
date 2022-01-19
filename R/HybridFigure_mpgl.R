### plot of hybrid zone similar to Fig 3 of Teeter etal 2010.Evolution 64:472–485.
rm(list=ls())
library(matrixStats);library(lattice)
pdf("output/HybridFigure_mpgl.pdf")
dat <- read.delim('data/genotypeLikelihoods.mpgl',sep=" ",header=F)
#dat <- t(dat[,-1]) # remove indiv IDS
meta <- read.table('data/meta_inds.txt',header=T,sep=",")
meta$reg <- factor(meta$reg)
meta$reg <- factor(meta$reg,levels = levels(meta$reg)[c(2,1,3,5,4,6)])
num.loci <- dim(dat)[1]
eNA <- as.matrix(dat[1:num.loci,meta$reg=="eUS"])
eNA.af <- rowSums(round(eNA))/(2*dim(eNA)[2])
pac <- as.matrix(dat[1:num.loci,meta$reg%in%c("Jap")])
pac.af <- rowSums(round(pac))/(2*dim(pac)[2])

#par(mfrow=c(2,1))
plot(abs(eNA.af-pac.af))

### take only those loci whose allele freqs differ by 0.95 (abs. percentage)
dat2 <- dat[1:num.loci,]
dat2 <- dat2[,order(meta$reg)]
dat2 <- as.matrix(dat2[abs(eNA.af-pac.af)>.6,])
dat2 <- round(dat2) ### rounds genotype frequencies to 0, 1 or 2
### switch alleles to same color
dat3 <- c()
diffloci <- (eNA.af-pac.af)[abs(eNA.af-pac.af)>.6]
for (i in 1:dim(dat2)[1])
{
  if(diffloci[i]<0) 
  {tmp <- 2-dat2[i,]
    dat3 <- rbind(dat3,tmp)}
  else 
    {dat3 <- rbind(dat3,dat2[i,])}
}

par(fig=c(.1,.70,.13,.87),mar=c(0,0,0,0))
image(dat3,xaxt="n",yaxt="n",main="")
segments(-1,82/349,2,82/349,lwd=2) ## eNA
segments(-1,103/349,2,103/349,lwd=2) ## ARG
segments(-1,126/349,2,126/349,lwd=2) ## HUM
segments(-1,177/349,2,177/349,lwd=2) ## SFB
segments(-1,296/349,2,296/349,lwd=2) ## JAP
mtext(at=c(40,90,114,150,236,325)/349,side=2,c("eastNA","Arg","Hum","SFB","Jap","westNA"),las=2)

##mtDNA pacific vs atlantic
par(fig=c(.71,.83,.1,.9),new=T)
mt <- ifelse(meta$mitoClade=="A",.9,ifelse(meta$mitoClade=="B",.5,.1))
mt <- mt[order(meta$reg)]
plot(x=mt,y=1:length(mt),type="l",xaxt="n",yaxt="n",xlab = '', ylab = '',
     xlim=c(0,1))
mtext(side=1,"COI allele",cex=.7)
mtext(side=1,at=c(.1,.5,.9),c("C","B","A"),cex=.7,line=-1)

# % atlantic alleles
par(fig=c(.84,.97,.1,.9),new=T)
s <- colSums(dat3)/(2*dim(dat3)[1])
plot(x=s,y=1:length(s),type="l",xaxt="n",yaxt="n",xlab = '', ylab = '',
     xlim=c(-.1,1.1))
mtext(side=1,"% Atlantic alleles",cex=.7)
mtext(side=1,at=c(.05,.95),c(0,1),cex=.7,line=-1)
#print(histogram(~s | meta$reg))

dev.off()


