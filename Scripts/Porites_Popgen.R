######## Agaricia agaricites population genetics
library(dplyr)
library(vegan)
library(ggplot2)
library(tess3r)
setwd("~/Documents/FLKeys/GitHub/Data") 

ll=load("Porites_env.RData") # Import environmental data
rownames(env)=env$Sample
env$Sample=NULL
IBS=as.matrix(read.table("Porites.ibsMat")) # Import genetic distance matrix
samples=rownames(env)
dimnames(IBS)=list(samples,samples)
inds=as.data.frame(samples)
colnames(inds) = "ind"
row.names(inds)<-inds$ind
Admixture=read.csv("Porites_admix.csv") # Import admixture proportions
Admixture=as.matrix(Admixture[,c(1:3)])
Cluster=read.csv("Porites_cluster.csv")[,-1] # Import admixture cluster assignment
sites=read.table("Porites_pops.txt") # Import site info
sites1 <- merge(inds, sites, by="ind", sort=F)
latlon=sites1[,c("Latitude","Longitude")]
names(latlon)=c("lat","lon")






#----------------------- Plot admixture barplots
my.colors <- c("goldenrod1", "blue4", "lightblue","olivedrab", "cyan3","hotpink","gold","orange")
my.palette <- CreatePalette(my.colors, 8)

# For K=2
qm=as.qmatrix(as.matrix(read.table("Porites_k2.qopt")))
bp=barplot(qm,border=NA,space=0,col.palette = my.palette)
axis(1, at = 1:nrow(qm), las = 3, cex.axis = .4) 


# For K=3
qm=as.qmatrix(as.matrix(read.table("Porites_k3.qopt")))
bp=barplot(qm,border=NA,space=0,col.palette = my.palette)
axis(1, at = 1:nrow(qm), las = 3, cex.axis = .4) 







#----------------------- Plot a PCoA
ord=capscale(IBS~1)
ords=scores(ord,display="sites")
axes2plot=c(1,2) # which PCAs to plot
scores=data.frame(ord$CA$u[,axes2plot])
summary(ord) # MDS1 and MDS2 proportion explained

# Plot colored by cluster
colvec <- c("goldenrod1", "lightblue", "blue4", "grey80")
cols=Cluster$cluster
plot(scores,col=cols+1,pch=16,main="Porites IBS")
points(scores, pch=16, cex=1, col=colvec[Cluster$cluster])

# Plot colored by depth
ggplot(scores,aes(scores[,1],scores[,2],colour=env$Depth)) + 
  geom_point(alpha=0.5)+
  scale_colour_gradient(low='skyblue', high='navy')+
  theme_bw()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  coord_equal()




#----------------------   Isolation by distance
geo <- sites1[, c("Latitude", "Longitude")]
geo_dist<- dist(geo)
ibs_dist<- as.dist(IBS)
#Plot IBS vs geo distance
ggplot(NULL, aes(x=geo_dist, y=ibs_dist) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis", trans="log") +
  geom_smooth(method="gam", formula = y ~s(x))+
  labs(title = "Porites- isolation by distance") +
  xlab('Geographic distance') +
  ylab('Genetic distance') +
  theme_bw()



#----------- Violin plot: depth v. cluster
Cluster$cluster <- as.factor(as.character(Cluster$cluster))
Cluster$cluster <- ifelse(Cluster$cluster == 1, 'Indigo',
                                      ifelse(Cluster$cluster == 2, 'Yellow', 
                                             ifelse(Cluster$cluster == 3, 'Blue', 'Hybrid')))
po.temp <- Cluster %>%
  mutate(cluster=factor(cluster,levels=c("Yellow", "Blue", "Indigo", "Hyrbid")) )
po.temp <- cbind(po.temp, env$Depth)
names(po.temp)[3]<-"Depth"
my_comparisons <- list( c("Yellow", "Indigo"), c("Indigo", "Blue"), c("Blue", "Yellow") )

cols=c("Yellow"="goldenrod1", "Blue"="lightblue", "Indigo"="blue4")
ggplot(po.temp, aes(x = cluster, y = Depth, fill = cluster)) + 
  geom_violin(scale = "width") + 
  scale_fill_manual(values = cols) +
  xlab('Cluster') +
  stat_compare_means(comparison = my_comparisons, label="p.signif") + # Add pairwise comparisons p-value
  stat_compare_means(label.y = 25, label.x = 1) +     # Add global p-value
  theme(axis.text = element_text(size = 20))  +
  theme_classic()


