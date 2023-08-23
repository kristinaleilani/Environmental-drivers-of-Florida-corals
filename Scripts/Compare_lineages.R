setwd("~/Documents/FLKeys/Scratch/Data") 
library(ggpubr)
library(tidyr)
library(maptools)
load("rasters_XY.RData") # Import raster data

# Load the shoreline  of the Florida Keys (download from www.ngdc.noaa.gov/mgg/shorelines)
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "gshhg-bin-2.3.6/gshhs_f.b"
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-83.5, -80), ylim = c(24.2, 25.7)) %>%
  fortify()

# Import turnover grids for each lineage
Ay<-read.csv("Agaricia_yellow_turnovers.csv")[,-1]
Ab<-read.csv("Agaricia_blue_turnovers.csv")[,-1]
Ai<-read.csv("Agaricia_indigo_turnovers.csv")[,-1]
Py<-read.csv("Porites_yellow_turnovers.csv")[,-1]
Pb<-read.csv("Porites_blue_turnovers.csv")[,-1]
Pi<-read.csv("Porites_indigo_turnovers.csv")[,-1]

# Function to measure distance between turnover grids
dist2dataframes=function(X,Y,method="euclidean"){
  if((nrow(X) == nrow(Y)) & (ncol(X) == ncol(Y))) {
    di=c()
    for (i in 1:nrow(X)) {
      #    message(i)
      dd=data.frame(rbind(X[i,],Y[i,]))
      di=c(di,dist(dd,method=method))
    }
    return(di)
  } else { 
    stop("dist2dataframes: dataframes are of different sizes!")
  }
}

# Change grid1 and grid2 depending on wihch lineages you are comparing
grid1=Ai
grid2=Pi

raster.vars=colnames(grid1)
dd=dist(grid1[,raster.vars])
qq=quantile(abs(dd),c(0.5,0.9),na.rm=T)
med=qq[1]
max=qq[2]

colnames(grid2[,raster.vars])

offset=dist2dataframes(grid1[,raster.vars],grid2[,raster.vars])
offset.2max=offset/qq[2]
offset.2max=cbind(XY, offset.2max)

mismatch_plot <- ggplot() +
  geom_raster(data = offset.2max, aes(x = x, y = y, fill = offset.2max)) + # add the raster 
  scale_fill_viridis(option="inferno")+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey90", color='black', lwd = 0.1)+
  guides(size="none") + 
  ggtitle("Ai_Pi") + # Change title depending on which lineages you are comparing
  theme(panel.background = element_rect(fill = "white", colour = "black"),  axis.title = element_blank()) 
mismatch_plot

# Save mismatch figure
save(mismatch_plot, offset.2max, file = "Ai_Pi.Rdata")





#######-------  Import mismatch data for every comparison ------######
load("Ay_Ab.Rdata")
Ay_Ab_plot<- mismatch_plot
Acomp<-offset.2max
names(Acomp)[3]<-"Ay_Ab"
load("Ay_Ai.Rdata")
Ay_Ai_plot<- mismatch_plot
Acomp<-cbind(Acomp, offset.2max[,3])
names(Acomp)[4]<-"Ay_Ai"
load("Ab_Ai.Rdata")
Ab_Ai_plot<- mismatch_plot
Acomp<-cbind(Acomp, offset.2max[,3])
names(Acomp)[5]<-"Ab_Ai"

load("Py_Pb.Rdata")
Py_Pb_plot<- mismatch_plot
Pcomp<-offset.2max
names(Pcomp)[3]<-"Py_Pb"
load("Py_Pi.Rdata")
Py_Pi_plot<- mismatch_plot
Pcomp<-cbind(Pcomp, offset.2max[,3])
names(Pcomp)[4]<-"Py_Pi"
load("Pb_Pi.Rdata")
Pb_Pi_plot<- mismatch_plot
Pcomp<-cbind(Pcomp, offset.2max[,3])
names(Pcomp)[5]<-"Pb_Pi"

load("Ay_Py.Rdata")
Ay_Py_plot<- mismatch_plot
APcomp<-offset.2max
names(APcomp)[3]<-"Ay_Py"
load("Ab_Pb.Rdata")
Ab_Pb_plot<- mismatch_plot
APcomp<-cbind(APcomp, offset.2max[,3])
names(APcomp)[4]<-"Ab_Pb"
load("Ai_Pi.Rdata")
Ai_Pi_plot<- mismatch_plot
APcomp<-cbind(APcomp, offset.2max[,3])
names(APcomp)[5]<-"Ai_Pi"

# Plot all mismatch maps together
ggarrange(Ay_Ab_plot, Ay_Ai_plot, Ab_Ai_plot, 
          Py_Pb_plot, Py_Pi_plot, Pb_Pi_plot, 
          Ay_Py_plot, Ab_Pb_plot, Ai_Pi_plot,
          common.legend = TRUE)






#####---------- Violin plot comparing lineages of Agaricia -----------#########

Acomp<- Acomp %>% pivot_longer(cols=c('Ay_Ab', 'Ay_Ai', 'Ab_Ai'),
                               names_to='Comparisons',
                               values_to='Offset')
Acomp<- Acomp %>%
  mutate( Comparisons=factor(Comparisons,levels=c('Ay_Ab', 'Ay_Ai', 'Ab_Ai')) )
a=ggplot(Acomp, aes(x = Comparisons, y = Offset)) + 
  geom_violin(scale = "width", draw_quantiles = 0.5) + 
  scale_fill_viridis(option="inferno")+
  xlab('Comparisons') +
  ylim(0.5, 2)+
  theme(axis.text = element_text(size = 20))  +
  theme_classic()





#####---------- Violin plot comparing lineages of Porites -----------#########

Pcomp<- Pcomp %>% pivot_longer(cols=c('Py_Pb', 'Py_Pi', 'Pb_Pi'),
                               names_to='Comparisons',
                               values_to='Offset')
Pcomp<- Pcomp %>%
  mutate( Comparisons=factor(Comparisons,levels=c('Py_Pb', 'Py_Pi', 'Pb_Pi')) )
p=ggplot(Pcomp, aes(x = Comparisons, y = Offset)) + 
  geom_violin(scale = "width", draw_quantiles = 0.5) + 
  scale_fill_viridis(option="inferno")+
  xlab('Comparisons') +
  ylim(0.5, 3.25)+
  theme(axis.text = element_text(size = 20))  +
  theme_classic()





#####---------- Violin plot comparing lineages Agaricia with Porites -----------#########

APcomp<- APcomp %>% pivot_longer(cols=c('Ay_Py', 'Ab_Pb', 'Ai_Pi'),
                                 names_to='Comparisons',
                                 values_to='Offset')
APcomp<- APcomp %>%
  mutate( Comparisons=factor(Comparisons,levels=c('Ay_Py', 'Ab_Pb', 'Ai_Pi')) )
ap=ggplot(APcomp, aes(x = Comparisons, y = Offset)) + 
  geom_violin(scale = "width", draw_quantiles = 0.5) + 
  scale_fill_viridis(option="inferno")+
  xlab('Comparisons') +
  ylim(0.5, 2)+
  theme(axis.text = element_text(size = 20))  +
  theme_classic()





#####---------- Plot all violin plots together -----------#########
ggarrange(a,p,ap, align = "hv")

