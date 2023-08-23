######## Agaricia agaricites population genetics
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(tess3r)
setwd("~/Documents/FLKeys/Agaricia") 

# Import genetic distance matrix
Y=as.matrix(read.table("Agaricia.ibsMat"))
# Import admixture proportions
Admixture=read.csv("Agaricia_admix.csv")
Admixture=as.matrix(Admixture[,c(1:3)])
# Import admixture cluster assignment
Cluster=read.csv("Agaricia_cluster.csv")[,-1]
# Import environmental data
load("Agaricia_env.RData")



#----------------------- Plot admixture barplots

# For K=2
dir="~/Documents/FLKeys/Agaricia" # path to input files
inName="Agaricia_k2.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
pops="Agaricia_pops.txt" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(inName,sep=""),header=F)
i2p=read.table(paste(pops,sep=""),fill = T, header=T)
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)

qm=as.qmatrix(as.matrix(read.table("Agaricia_k2.qopt")))
bp=barplot(qm,border=NA,space=0)
axis(1, at = 1:nrow(qm), labels = i2p$pop[bp$order], las = 3, cex.axis = .4) 


# For K=3
dir="~/Documents/FLKeys/Agaricia" # path to input files
inName="Agaricia_k3.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
pops="Agaricia_pops.txt" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(inName,sep=""),header=F)
i2p=read.table(paste(pops,sep=""),fill = T, header=T)
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
head(tbl)


my.colors <- c("goldenrod1", "lightblue", "blue4")
my.palette <- CreatePalette(my.colors, 8)
qm=as.qmatrix(as.matrix(read.table("Agaricia_k3.qopt")))
bp=barplot(qm,border=NA,space=0,col.palette = my.palette)
axis(1, at = 1:nrow(qm), labels = i2p$pop[bp$order], las = 3, cex.axis = .4) 






#----------------------- Plot a PCoA
ord=capscale(as.dist(Y)~1)
ords=scores(ord,scaling=1,display="sites")
axes2plot=c(1,2) # which PCAs to plot
scores=data.frame(ord$CA$u[,axes2plot])

summary(ord) # Proportion of genetic variation explained by MDS1 and MDS2 

# Plot colored by cluster
colvec <- c("goldenrod1", "lightblue", "blue4")
cols=Cluster$cluster
axes2plot=c(1,2)
scores=scores(ord,scaling=1,display="sites",choices=axes2plot)
plot(scores,pch=16,col=colvec[Cluster$cluster],asp=1)
abline(h=0,lty=3,col="grey60")
abline(v=0,lty=3,col="grey60")

# Plot colored by depth
ggplot(scores,aes(scores[,1],scores[,2],colour=env$Depth)) + 
  geom_point(alpha=0.5)+
  scale_colour_gradient(low='skyblue', high='navy')+
  theme_bw()+
  xlab(names(scores)[1])+
  ylab(names(scores)[2])+
  coord_equal()



#----------------------   Isolation by distance
meta <- read.csv("Agaricia_meta.csv", row.names=1)
geo <- meta[, c("Latitude", "Longitude")]
geo_dist<- dist(geo)
ibs_dist<- as.dist(Y)
#Plot IBS vs geo distance
ggplot(NULL, aes(x=geo_dist, y=ibs_dist) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis", trans="log") +
  geom_smooth(method="gam", formula = y ~s(x))+
  labs(title = "Agaricia- isolation by distance") +
  xlab('Geographic distance') +
  ylab('Genetic distance') +
  theme_bw()



#----------- Violin plot: depth v. cluster
metad <- cbind(meta, env$Depth)
names(metad)[4] <- "Depth"
metad$cluster.admix <- as.factor(as.character(metad$cluster.admix))
metad$cluster.admix <- ifelse(metad$cluster.admix == 1, 'Yellow',
                                      ifelse(metad$cluster.admix == 2, 'Blue', 'Indigo'))
ag.temp <- metad %>%
  mutate(cluster.admix=factor(cluster.admix,levels=c("Yellow", "Blue", "Indigo")) )
my_comparisons <- list( c("Yellow", "Indigo"), c("Indigo", "Blue"), c("Yellow", "Blue") )

cols=c("Yellow"="goldenrod1", "Blue"="lightblue", "Indigo"="blue4")
ggplot(ag.temp, aes(x = cluster.admix, y = Depth, fill = cluster.admix)) + 
  geom_violin(scale = "width") + 
  scale_fill_manual(values = cols) +
  xlab('Cluster') +
  stat_compare_means(comparison = my_comparisons, label="p.signif") + # Add pairwise comparisons p-value
  stat_compare_means(label.y = 25, label.x = 1) +     # Add global p-value
  theme(axis.text = element_text(size = 20))  +
  theme_classic()


