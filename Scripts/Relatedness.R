library(dplyr)
library(vegan)
library(lares)
library(ggplot2)
library(ggpubr)
setwd("~/Documents/FLKeys/GitHub/Data") 
# Note: relatedness files include all samples passing sequencing quality filtering, and include potential clone and sibling groups



########--------- Relatedness in Porites astreoides -----------###################
rel=read.table("Porites_rel.res",sep="\t",header=T) # Import genetic distance matrix (relatedness), output of ngsRelate
Cluster=read.csv("Porites_cluster.csv")[,-1] # Import admixture cluster assignment
# creating an empty square matrix
relm=matrix(0,nrow=length(unique(rel$a))+1,ncol=length(unique(rel$a))+1)
# filling up the square matrix with entries from "rab" column
for (a in unique(rel$a)) {
  for (b in unique(rel$b)) {
    if (b<=a) {next}
    relm[a+1,b+1]=relm[b+1,a+1]=rel[rel$a==a & rel$b==b,"rab"]
  }
}
diag(relm)=1
poritesZ=(1-relm)
poritesZ.ord=capscale(poritesZ~1)
poritesZ.ords=scores(poritesZ.ord,display="sites")
axes2plot=c(1,2) # which PCAs to plot
Zscores=data.frame(poritesZ.ord$CA$u[,axes2plot])
cols=env.porites0$k3cluster.admix
colvec <- c("blue4", "goldenrod1","lightblue", "grey60")
plot(Zscores,col=cols+1,pch=16,main="Porites Relatedness")
points(Zscores, pch=16, cex=1, col=colvec[Cluster$cluster])
row.names(poritesZ)=Cluster[,1]
hc=hclust(as.dist(poritesZ),"ave")
hcd <- as.dendrogram(hc, hang = 0.1)
dend1 <- color_branches(hcd, k = 3, col=c("lightblue", "goldenrod1","blue4"))
plot(dend1, ylab = "Relatedness", main="Cluster Dendrogram", leaflab = "none")







########--------- Relatedness in Agaricia agaricites -----------###################
# Must explore each lineage separately due to high differentiation between lineages
yrel=read.table("Ag_yellow.res",sep="\t",header=T) # Import genetic distance matrix (relatedness), output of ngsRelate
brel=read.table("Ag_blue.res",sep="\t",header=T)
irel=read.table("Ag_indigo.res",sep="\t",header=T)

# For Agaricia yellow:
# creating an empty square matrix
yrelm=matrix(0,nrow=length(unique(yrel$a))+1,ncol=length(unique(yrel$a))+1)
# filling up the square matrix with entries from "rab" column
for (a in unique(yrel$a)) {
  for (b in unique(yrel$b)) {
    if (b<=a) {next}
    yrelm[a+1,b+1]=yrelm[b+1,a+1]=yrel[yrel$a==a & yrel$b==b,"rab"]
  }
}
diag(yrelm)=1
AyellowZ=(1-yrelm)
hc=hclust(as.dist(AyellowZ),"ave")
hcd <- as.dendrogram(hc, hang = 0.1)
dend1 <- color_branches(hcd, k = 1, col="goldenrod1")
plot(dend1, ylab = "Relatedness", main="Agaricia yellow Cluster Dendrogram", leaflab = "none")



# For Agaricia blue:
# creating an empty square matrix
brelm=matrix(0,nrow=length(unique(brel$a))+1,ncol=length(unique(brel$a))+1)
# filling up the square matrix with entries from "rab" column
for (a in unique(brel$a)) {
  for (b in unique(brel$b)) {
    if (b<=a) {next}
    brelm[a+1,b+1]=brelm[b+1,a+1]=brel[brel$a==a & brel$b==b,"rab"]
  }
}
diag(brelm)=1
AblueZ=(1-brelm)
hc=hclust(as.dist(AblueZ),"ave")
hcd <- as.dendrogram(hc, hang = 0.1)
dend1 <- color_branches(hcd, k = 1, col="lightblue")
plot(dend1, ylab = "Relatedness", main="Agaricia blue Cluster Dendrogram", leaflab = "none")


# For Agaricia indigo:
# creating an empty square matrix
irelm=matrix(0,nrow=length(unique(irel$a))+1,ncol=length(unique(irel$a))+1)
# filling up the square matrix with entries from "rab" column
for (a in unique(irel$a)) {
  for (b in unique(irel$b)) {
    if (b<=a) {next}
    irelm[a+1,b+1]=irelm[b+1,a+1]=irel[irel$a==a & irel$b==b,"rab"]
  }
}
diag(irelm)=1
AindigoZ=(1-irelm)
hc=hclust(as.dist(AindigoZ),"ave")
hcd <- as.dendrogram(hc, hang = 0.1)
dend1 <- color_branches(hcd, k = 1, col="blue4")
plot(dend1, ylab = "Relatedness", main="Agaricia indigo Cluster Dendrogram", leaflab = "none")
