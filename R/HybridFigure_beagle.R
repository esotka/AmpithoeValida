### plot of hybrid zone similar to Fig 3 of Teeter etal 2010.Evolution 64:472–485.
rm(list=ls())
library(matrixStats);library(lattice); library(parallel)
pdf("output/HybridFigure_beagle.pdf")
# convert calls - it takes a long while; do it once
#source('R/beagle_calls.R')
#Acalls <- beagle_gl(fn="data/261_indiv.beagle.gz",bamlist= "data/inds261",rowstart=0,nrows=-1,support=log(3))
#write.table(Acalls,"data/261_indiv.beagleCalls_log3.txt",quote=F,row.names=F,sep="\t")

dat <- read.delim("data/261_indiv.beagleCalls_log3.txt",sep="\t",header=T)
rownames(dat) <- dat[,1]; dat <- dat[,-1]
n <- dim(dat)[2]
### keep loci with calls at 90% or more individuals
tmp <- rowSums(is.na(dat))
dat2 <- dat[tmp<=(.2*n),] # no missing data

## convert to numbers
dat2.tr <- data.frame(t(dat2))
dat2.num <- c()
for (i in 1:dim(dat2.tr)[2])
{
  tmp <- as.numeric(as.factor(dat2.tr[,i]))-1
  dat2.num <- cbind(dat2.num,tmp)
}

dat2.num <- t(dat2.num)
colnames(dat2.num) <- colnames(dat)

meta <- read.table('data/meta_inds.txt',header=T,sep=",")
meta <- meta[meta$ind%in%colnames(dat),]
inds <- meta$ind
meta$reg <- factor(meta$reg)
meta$reg <- factor(meta$reg,levels = levels(meta$reg)[c(2,1,3,5,4,6)])
num.loci <- dim(dat2.num)[1]


eNA <- as.matrix(dat2.num[1:num.loci,meta$reg=="eUS"])
eNA.n <- (dim(eNA))[2]-rowSums(is.na(eNA))
eNA.af <- rowSums(eNA,na.rm=T)/(2*eNA.n)
pac <- as.matrix(dat2.num[1:num.loci,meta$reg%in%c("Jap","wUS")])
pac.n <- (dim(pac))[2]-rowSums(is.na(pac))
#pac.af <- rowSums(round(pac))/(2*dim(pac)[2])
pac.af <- rowSums(pac,na.rm=T)/(2*pac.n)

#par(mfrow=c(2,1))
plot(abs(eNA.af-pac.af))

### take only those loci whose allele freqs differ by 0.8 (abs. percentage)
dat2.num <- dat2.num[abs(eNA.af-pac.af)>.9,]
dat2.num <- dat2.num[,order(meta$reg)]
#dat2.num <- as.matrix(dat2.num[abs(eNA.af-pac.af)>.8,])
image(dat2.num,xaxt="n",yaxt="n",main="")


#dat2 <- round(dat2) ### rounds genotype frequencies to 0, 1 or 2
### switch alleles to same color
dat3 <- c()
diffloci <- (eNA.af-pac.af)[abs(eNA.af-pac.af)>.9]
for (i in 1:dim(dat2.num)[1])
{
  if(diffloci[i]<0) 
  {tmp <- 2-dat2.num[i,]
    dat3 <- rbind(dat3,tmp)}
  else 
    {dat3 <- rbind(dat3,dat2.num[i,])}
}

par(fig=c(.1,.70,.13,.87),mar=c(0,0,0,0))
image(dat3,xaxt="n",yaxt="n",main="")


regnames <- meta$reg[order(meta$reg)]
nums <- 1:length(meta$reg)
for (n in 1:length(levels(meta$reg)))
{
  tmp <- nums[regnames==levels(meta$reg)[n]]
  segments(-1,max(tmp)/length(meta$reg),2,max(tmp)/length(meta$reg),lwd=2)
  print(mean(tmp))
  mtext(at=mean(tmp)/length(meta$reg),side=2,levels(meta$reg)[n],las=2)
}
#segments(-1,82/349,2,82/349,lwd=2) ## eNA
#segments(-1,103/349,2,103/349,lwd=2) ## ARG
#segments(-1,126/349,2,126/349,lwd=2) ## HUM
#segments(-1,177/349,2,177/349,lwd=2) ## SFB
#segments(-1,296/349,2,296/349,lwd=2) ## JAP
#mtext(at=c(40,90,114,150,236,325)/349,side=2,c("eastNA","Arg","Hum","SFB","Jap","westNA"),las=2)

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
s <- colSums(dat3,na.rm=T)/(2*dim(dat3)[1])
plot(x=s,y=1:length(s),type="l",xaxt="n",yaxt="n",xlab = '', ylab = '',
     xlim=c(-.1,1.1))
mtext(side=1,"% Atlantic alleles",cex=.7)
mtext(side=1,at=c(.05,.95),c(0,1),cex=.7,line=-1)
#print(histogram(~s | meta$reg))

dev.off()


