library("RColorBrewer");library("adegenet")
### plot PCAs
rm(list=ls())

dat <- read.csv('data/genotypeLikelihoods.txt',header=T)
meta <- read.table('data/meta_inds.txt',header=T,sep=",")
dat2 <- dat[,-1]
pdf('output/PCA.pdf')

pc.cr <- prcomp(dat2)
print(summary(pc.cr))
col.reg <- brewer.pal(6,"Dark2")[meta$reg]
plot(pc.cr$x[,1],pc.cr$x[,2],type="n")
text(pc.cr$x[,1],pc.cr$x[,2],col=col.reg,pch=20,cex=.7,meta$reg)

# within the US - remove recent immigrants (#1 and 54): ALN001 and CBV013
nat <- dat2[meta$reg%in%c("eUS"),]
nat <- nat[!rownames(nat)%in%c(1,54),]
nat.pop <- meta$pop[rownames(meta)%in%rownames(nat)];nat.pop <- factor(nat.pop)
pc.cr <- prcomp(nat)
print(summary(pc.cr))
col.reg <- brewer.pal(6,"Dark2")[nat.pop]
plot(pc.cr$x[,1],pc.cr$x[,2],type="n")
text(pc.cr$x[,1],pc.cr$x[,2],col=col.reg,pch=20,cex=.7,nat.pop)

# within the eUS + argentina
nat <- dat2[meta$reg%in%c("eUS","Arg"),]
nat <- nat[!rownames(nat)%in%c(1,54),]
nat.pop <- meta$pop[rownames(meta)%in%rownames(nat)];nat.pop <- factor(nat.pop)
pc.cr <- prcomp(nat)
summary(pc.cr)
col.reg <- brewer.pal(6,"Dark2")[nat.pop]
plot(pc.cr$x[,1],pc.cr$x[,2],type="n")
text(pc.cr$x[,1],pc.cr$x[,2],col=col.reg,pch=20,cex=.7,nat.pop)

# within the wUS - remove invasive (#313 TMB004)

nat <- dat2[meta$reg%in%c("wUS"),]
nat <- nat[!rownames(nat)%in%c(313),]
nat.pop <- meta$pop[rownames(meta)%in%rownames(nat)];nat.pop <- factor(nat.pop)
pc.cr <- prcomp(nat)
summary(pc.cr)
col.reg <- brewer.pal(6,"Dark2")[nat.pop]
plot(pc.cr$x[,1],pc.cr$x[,2],type="n")
text(pc.cr$x[,1],pc.cr$x[,2],col=col.reg,pch=20,cex=.7,nat.pop)

## Clade C
nat <- dat2[meta$pop%in%c("MNG","MOU","SOU","EHS","TMB","WIL"),]
nat <- nat[!rownames(nat)%in%c(313),]
nat.pop <- meta$pop[rownames(meta)%in%rownames(nat)];nat.pop <- factor(nat.pop)
pc.cr <- prcomp(nat)
summary(pc.cr)
col.reg <- brewer.pal(6,"Dark2")[nat.pop]
plot(pc.cr$x[,1],pc.cr$x[,2],type="n")
text(pc.cr$x[,1],pc.cr$x[,2],col=col.reg,pch=20,cex=.7,nat.pop)

dev.off()