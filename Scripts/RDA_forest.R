library(ggplot2)
library(vegan)
library(gradientForest)
library(dplyr)
library(maptools)
setwd("~/Documents/Data") 
source("~/Documents/RDA-forest_functions.R")


# Load the shoreline  of the Florida Keys (download from www.ngdc.noaa.gov/mgg/shorelines)
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "gshhg-bin-2.3.6/gshhs_f.b"
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-83.5, -80), ylim = c(24.2, 25.7)) %>%
  fortify()






#####--------- Import species or lineage files ---------#######

###-------------------- For Agaricia -------------########
species <- "Agaricia"
ll=load("Agaricia_env.RData")
IBS=as.matrix(read.table("Agaricia.ibsMat"))
samples=rownames(env)
dimnames(IBS)=list(samples,samples)

inds=as.data.frame(samples)
colnames(inds) = "ind"
row.names(inds)<-inds$ind

###---- For Agaricia yellow lineage, >90% admixture assignment to yellow cluster
lineage<-"yellow"
inds=read.table("Agyellow",sep="\t")
colnames(inds) = "ind"
row.names(inds)<-inds$ind
env=env[row.names(inds),]

IBS=as.matrix(read.table("Agyellow.ibsMat"))
samples=rownames(inds)
dimnames(IBS)=list(samples,samples)


###---- For Agaricia blue lineage, >60% admixture assignment to blue cluster
lineage<-"blue"
inds=read.table("Agblue",sep="\t")
colnames(inds) = "ind"
row.names(inds)<-inds$ind
env=env[row.names(inds),]

IBS=as.matrix(read.table("Agblue.ibsMat"))
samples=rownames(inds)
dimnames(IBS)=list(samples,samples)


###---- For Agaricia indigo lineage, >90% admixture assignment to indigo cluster
lineage<-"indigo"
inds=read.table("Agindigo",sep="\t")
colnames(inds) = "ind"
row.names(inds)<-inds$ind
env=env[row.names(inds),]

IBS=as.matrix(read.table("Agindigo.ibsMat"))
samples=rownames(inds)
dimnames(IBS)=list(samples,samples)



###---- Import admixture proportions between 3 clusters 
admix0=read.csv("Agaricia_admix.csv")[,-1]
rownames(admix0)<-admix0$Sample
admix<-admix0[row.names(inds),]


###---- Import site info 
sites=read.table("Agpop_coords.txt")
sites1 <- merge(inds, sites, by="ind", sort=F)
rownames(sites1)<-sites1$ind
latlon=sites1[,c("Latitude","Longitude")]
names(latlon)=c("lat","lon")






###-------------------- For Porites -------------########
species <- "Porites"
ll=load("Porites_env.RData")
rownames(env)=env$Sample
env$Sample=NULL
IBS=as.matrix(read.table("Porites.ibsMat"))
samples=rownames(env)
dimnames(IBS)=list(samples,samples)


# Porites yellow lineage (all lineages admix > 0.6)
lineage<-"yellow"
inds=read.table("Poyellow",sep="\t")
colnames(inds) = "ind"
row.names(inds)<-inds$ind
env=env[row.names(inds),]

IBS=as.matrix(read.table("Poyellow.ibsMat"))
samples=rownames(inds)
dimnames(IBS)=list(samples,samples)



# Porites blue lineage
lineage<-"blue"
inds=read.table("Poblue",sep="\t")
colnames(inds) = "ind"
row.names(inds)<-inds$ind
env=env[row.names(inds),]

IBS=as.matrix(read.table("Poblue.ibsMat"))
samples=rownames(inds)
dimnames(IBS)=list(samples,samples)


# Porites indigo lineage
lineage<-"indigo"
inds=read.table("Poindigo",sep="\t")
colnames(inds) = "ind"
row.names(inds)<-inds$ind
env=env[row.names(inds),]

IBS=as.matrix(read.table("Poindigo.ibsMat"))
samples=rownames(inds)
dimnames(IBS)=list(samples,samples)




###---- Import admixture proportions between 3 clusters 
admix0=read.csv("Porites_admix.csv")[,-1]
rownames(admix0)<-admix0$Sample
admix<-admix0[row.names(inds),]


###---- Import site info 
sites=read.table("Porites_pops.txt")
sites1 <- merge(inds, sites, by="ind", sort=F)
rownames(sites1)<-sites1$ind
latlon=sites1[,c("Latitude","Longitude")]
names(latlon)=c("lat","lon")







#####--------- Check pop structure of lineages and remove clones if necessary ---------#######
hc=hclust(as.dist(IBS),"ave")
plot(hc,cex=0.5) # clustering of samples by IBS (great to detect clones or closely related individuals)
abline(h=0.15,col="red") # this seems like a reasonable  "low-hang" threshold for calling related groups

cuts=cutree(hc,h=0.15)
goods=c();i=1
for (i in unique(cuts)) {
  goods[i]=names(cuts[cuts==i])[1]
}
length(goods)  # how many samples are left?
# Subsetting all data for only the retained samples
IBS=IBS[goods,goods] 
latlon=latlon[goods,]
admix=admix[goods,]
env=env[goods,]
inds=inds[goods,]
sites1=sites1[goods,]

# Evaluate PCoA for outliers before continuing
ord.all=capscale(as.dist(IBS)~1)
plot(ord.all,scaling=1)
points(ord.all,scaling=1, pch=16)
summary(ord.all)





#####--------- Create spatial predictors  ---------#######
# principal components of space - lat and lon rotated to be uncorrelated, and centered 
xy=scores(capscale(dist(latlon)~1),scaling=1)$sites
colnames(xy)=c("xx","yy")
# Moran eigenvector maps (MEMs) - capturing possible spatial trends
mems=data.frame(pcnm(dist(xy))$vectors)
# adding xy and first 5 MEMs (will be called PCNM1-5) to env
env=cbind(env,xy)
env=cbind(env,mems[,1:5])
colnames(env)
# remembering the names of spatial predictors
space=colnames(env)[grep("PCNM|xx|yy",colnames(env))]






#####--------- Variable selection  ---------#######
# See https://github.com/z0on/RDA-forest for more explanation

admix.cov=admix[1:3]
admix.cov=admix.cov[,!(names(admix.cov) %in% lineage)] # Select admixture columns for correcting population structure
covars=cbind(admix.cov[,1], admix.cov[,2])
mm=mtrySelection(Y=IBS,X=env,nreps=11,covariates=covars, prop.positive.cutoff=0.5,top.pcs=25)

# boxplot of importance differences at different mtry 
ggplot(mm$delta,aes(var,values))+
  geom_boxplot()+
  coord_flip()+
  geom_hline(yintercept=0,col="red")

# bar chart of proportion of positive change in response to higher mtry, good variables would be the ones above the red line
ggplot(mm$prop.positive,aes(var,prop.positive))+
  geom_bar(stat="identity")+
  coord_flip()+
  geom_hline(yintercept=0.5,col="red")





#####--------------- Spatial bootstrap of mtry-passing variables ---------#######

ll=load("rasters_XY.RData") # load table of environmental values for a grid of spatial points where we want to predict how our creatures would adapt.
ll
# "rasters" "XY"
plot(XY,pch=".",asp=1) # view the grid (Florida Keys seascape)

# if there are no new points to predict, just skip the newX option in the call to spatialBootstrap; 
# the predictions will be made for original data points then.
sb=spatialBootstrap(Y=IBS,X=env,newX=rasters,nreps=25,covariates=covars,top.pcs=25)

# plot importance boxplot including space variables
ggplot(sb$all.importances,aes(variable,importance))+geom_boxplot()+coord_flip()+theme_bw()
sum(sb$median.importance) #proportion of variation explained

# plot importance boxplot without space variables
ggplot(sb$all.importances[!(sb$all.importances$variable %in% space),],aes(variable,importance))+geom_boxplot()+coord_flip()+theme_bw()
sum(sb$median.importance[!(names(sb$median.importance) %in% space)]) #proportion of variation explained without spatial variables







#####------------------- Plotting genetic turnover (adaptation) map ---------#######
# see https://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf for more explanation

turnovers=sb$turnovers
write.csv(turnovers, file=paste(species, lineage, "turnovers.csv", sep='_'))

raster.vars=colnames(turnovers)

# principal component analysis of the predicted genetic turnover patterns
pc <- prcomp(turnovers)
plot(pc$sdev)
pcs2show=c(1,2,3)

# color flippage flags - change between -1 and 1 to possibly improve color representation in the final map
flip.pc1=(1)
flip.pc2=(1)
flip.pc3=(1)

flip=""
if(flip.pc1==(-1)) { flip=paste(flip,"1",sep="")}
if(flip.pc2==(-1)) { flip=paste(flip,"2",sep="")}
if(flip.pc3==(-1)) { flip=paste(flip,"3",sep="")}

# magic to color three PCA dimensions on a map
pc1 <- flip.pc1*pc$x[, pcs2show[1]]
pc2 <- flip.pc2*pc$x[, pcs2show[2]]
pc3 <- flip.pc3*pc$x[, pcs2show[3]]
b <- pc1 - pc2
g <- -pc1
r <- pc3 + pc2 - pc1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
nvs <- dim(pc$rotation)[pcs2show[1]]
lv <- length(raster.vars)
vind <- rownames(pc$rotation) %in% raster.vars
scal <- 30
xrng <- range(pc$x[, 1], pc$rotation[, pcs2show[1]]/scal) * 1.9
yrng <- range(pc$x[, 2], pc$rotation[, pcs2show[2]]/scal) * 1.9
man.colors=rgb(r, g, b, max = 255)

# -------- plotting the map of predicted adaptive communities

coltxt="coral" 
important=bests[1:3]
lv=length(important)
par(mfrow=c(1,2))
plot((pc$x[, pcs2show[1:2]]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1,xaxt="none",yaxt="none",bty="none",xlab="",ylab="")
jit <- rnorm(lv,0.005,0.001)
arrows(rep(0, lv), rep(0, lv), pc$rotation[important, pcs2show[1]]/scal, pc$rotation[important, pcs2show[2]]/scal, length = 0.0625,col=coltxt)
text(pc$rotation[important, 1]/scal + jit * sign(pc$rotation[important, pcs2show[1]]), pc$rotation[important, pcs2show[2]]/scal + jit * sign(pc$rotation[important, pcs2show[2]]), labels = important,cex=0.7,col=coltxt)
plot(XY, pch=15,cex = 0.5, asp = 1, col = man.colors)
map(coasts,add=T,col="grey80",fill=T,border="grey80",lwd=1)

# contrasting colors = habitats requiring differential adaptation, likely driven by factors in the legend.


