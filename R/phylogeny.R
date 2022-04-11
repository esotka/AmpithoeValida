library("ape")
### draw Bayes tree and add ML and Bayes support

rm(list=ls())
tr <- read.tree("data/bayes.tr")
hap <- readLines("data/ListOfDuplicateSequences_665ind.txt",n=52)
# hap names and count number of haplotypes
tmp <- c()
n <- c()
for(i in 1:length(hap))
{
  tmp <- c(tmp,substr(hap[grep(tr$tip.label[i],hap)],1,6))
  n <- c(n,length(unlist(strsplit(hap[grep(tr$tip.label[i],hap)]," ")))-1)
  }
tmp[grep("KP316287",hap)] <- "Ampithoe_caddi"
tmp[grep("KP316288",hap)] <- "Ampithoe_dalli"
tmp[grep("KP316293",hap)] <- "Ampithoe_longimana"
tmp[grep("KP316294",hap)] <- "Ampithoe_lacertosa"
tmp[grep("KP316301",hap)] <- "Ampithoe_sectimanus" 
tmp2 <- ifelse(n==1,
               tmp,
               paste(tmp," (",n,")",sep=""))

pdf("output/phylogeny.pdf",width=5,height=5)
par(mar=c(1,1,1,1))
plot(tr,show.tip.label = F,x.lim=c(-.01,0.3))

tiplabels(tmp2,frame = "none",adj=-.1,cex=.5)

# number, bayes, ml bs
nodes <- c(70,100,70,
           71,97,70,
           72,99,0,
           54,100,100,
           69,98,100,
           63,93,69,
           64,99,88)
nodes <- as.data.frame(matrix(nodes,ncol=3,byrow=T))          
colnames(nodes) <- c("node","bayes","ml_bs")
nodelabels(node=nodes$node,pch=19,cex=ifelse(nodes$bayes>=95,1,0))
nodelabels(node=nodes$node,pch=20,cex=ifelse(nodes$ml_bs>=70,1,0),col="grey")
text(0.14,42,"Clade A - Atlantic",adj=0,col="orange")
segments(0.13,47,0.13,36,lwd=5,col="orange")
text(0.13,28,"Clade B - Pacific",adj=-.05,col="darkgreen")
segments(0.125,35,0.125,19,lwd=5,col="darkgreen")
text(0.11,9,"Clade C - Pacific",adj=-0,col="darkblue")
segments(0.105,18,0.105,1,lwd=5,col="darkblue")

dev.off()


