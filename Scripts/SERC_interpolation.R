library(tidyverse)
library(maptools)
library(raster)
library(dplyr)
library(automap)
setwd("~/Documents/Flkeys/SERC/") 

# Load and format raw data
wq.data <- read.csv('SERC.csv', na.strings=c("","NA","ND")) %>%
  # Reformat dates and station IDs
  mutate(newdate1 = as.Date(as.character(DATE), format="%m/%d/%y"),
         newdate2 = as.Date(as.character(DATE), format="%d-%b-%y"),
         STATION = gsub("i", "", gsub("0\\+", "", gsub("\\+0i", "", STATION)))) %>%
  # Merge dates from the two date formats
  mutate(newdate.merge = as.Date(ifelse(!is.na(newdate1), as.character(newdate1), as.character(newdate2)))) %>%
  # Delineate data by YR-MON and YR-WK
  mutate(YR_WK = as.character(format(newdate.merge, "%Y-%V")), YR_MON = as.character(format(newdate.merge, "%Y-%m"))) %>%
  # Remove rows without dates and that are outside our region of interest
  filter(!is.na(newdate.merge),
         SEGMENT %in% c('OFF', 'MAR', 'LK', 'MK', 'UK', 'BKS', 'BKB', 'WFB', 'SFB', 'EFCB', 'MBS',
                        'CS', 'SCM', 'SCI', 'NCI', 'NCO', 'SNB', 'NNB')) %>%
  # Fix formatting error in TURB.S variable
  mutate(TURB.S = as.numeric(gsub("\\.\\.", "\\.", as.character(TURB.S))))

wq.2<- wq.data %>%
  separate(YR_MON, c("Year", "Month"), "-")
wq.data<-wq.2
## Discrete sampling periods (224 unique sites):
# 1123 days
# 470 weeks
# 217 months
# 98 designated surveys
# 26 years

if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "gshhg-bin-2.3.6/gshhs_f.b"
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-83.5, -80), ylim = c(24.2, 25.7)) %>%
  fortify()



# Plot survey sites for every month in the sampling range
for(mon in unique(wq.data$YR_MON)){
  map.data <- subset(wq.data, YR_MON == mon)
  sampling.map <- ggplot()+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
    ggtitle(zoo::as.yearmon(mon))+
    geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
    geom_point(data=map.data, aes(x=LONDEC, y=LATDEC), size=0.75, color = 'red')+
    scale_color_brewer(palette = 'Set1')+
    coord_fixed(xlim = c(-83.1,-79.9), ylim = c(24.2,26), expand = 0)
  cowplot::save_plot(paste0('~/Documents/Flkeys/SERC/sites_', mon, '.png'),
                     plot = sampling.map, base_width = 5, base_height = 3)
}

# Visualize the survey schedule and frequency
site.time.df <- base::expand.grid(date = seq(min(wq.data$newdate.merge), max(wq.data$newdate.merge), by="days"),
                                  site = unique(wq.data$SITE)) %>%
  mutate(sample = paste0(date, ':', site) %in% paste0(wq.data$newdate.merge, ':', wq.data$SITE),
         y = as.numeric(site)) %>%
  group_by(site) %>%
  dplyr::mutate(x = seq(1:n())) 
sampling.grid <- ggplot(site.time.df)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none')+
  geom_tile(aes(x = x, y = y, fill = sample), col = NA)+
  scale_fill_manual(values = c('white', 'navyblue'))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))
#sampling.grid
cowplot::save_plot(filename = 'sampling_times.pdf',
                   plot = sampling.grid, base_height = 1, base_width = 3)

# Identify and remove sites lacking sufficient data points
total.samples <- group_by(wq.data, SITE) %>%
  dplyr::summarise(n = n())
ggplot(total.samples)+
  aes(x = n) +
  geom_histogram(bins = 30, fill = "#0c4c8a") +
  theme_minimal()
wq.data.filt <- filter(wq.data, SITE %in% total.samples[total.samples$n > 50,]$SITE)


# Summarize data for each site
wq.summ <- group_by(wq.data.filt, SITE) %>%
  dplyr::summarise(across(1:10, first),
                   across(11:54, list(mean = function(x) mean(x, na.rm = T), 
                                      #median = function(x) median(x, na.rm = T),
                                      #sd = function(x) sd(x, na.rm = T),
                                      min = function(x) min(x, na.rm = T),
                                      max = function(x) max(x, na.rm = T))))
                                      #n = function(x) sum(!is.na(x))), .names = "{col}.{fn}"))
wq.summ$LONDEC[wq.summ$SITE == "Middle Ground"] <- -81.8917

# Calculate monthly range for every variable
wq.month <- group_by(wq.data, SITE, Month) %>%
  dplyr::summarise(across(1:10, first),
                   across(11:54, list(monthly = function(x) (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))))
wq.month2 <- group_by(wq.month, SITE) %>%
  dplyr::summarise(across(1:10, first),
                   across(11:54, list(range = function(x) mean(x, na.rm = T))))
wq.month2 <- wq.month2[, -c(2:12)]

# Calculate yearly range for every variable
wq.year <- group_by(wq.data, SITE, Year) %>%
  dplyr::summarise(across(1:10, first),
                   across(11:54, list(yearly = function(x) (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))))
wq.year2 <- group_by(wq.year, SITE) %>%
  dplyr::summarise(across(1:10, first),
                   across(11:54, list(range = function(x) mean(x, na.rm = T))))
wq.year2 <- wq.year2[, -c(2:12)]
wq.summ2 <- left_join(wq.summ, wq.month2, by="SITE")                                   
wq.summ3 <- left_join(wq.summ2, wq.year2, by="SITE")                                   
wq.summ<-wq.summ3

wq.summ[wq.summ == Inf] <- NA
wq.summ[wq.summ == -Inf] <- NA


#Check how many long/lat are NA
row.has.na <- apply(wq.summ[,10:11], 1, function(x){any(is.na(x))})
sum(row.has.na) 
which(is.na(wq.summ[,10:11]), arr.ind=TRUE) 
wq.summ<-wq.summ[-c(47), ] # remove row with NA

# Example plot
ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=wq.summ, aes(x=LONDEC, y=LATDEC, color = NOX.S_mean))+
  scale_color_viridis_c(na.value = NA)+
  coord_fixed(xlim = c(-83.3,-79.9), ylim = c(24.2,26), expand = 0)



# Import polygon of Florida Keys
map.hull.stretch <- read.table('map_hull.txt', sep = " ", header = FALSE)
map.hull.stretch <- as.data.frame(map.hull.stretch)

# Convert hull to spatial object
hull.poly <- SpatialPolygons(list(Polygons(list(Polygon(cbind(map.hull.stretch$V1, map.hull.stretch$V2))), ID=1)), proj4string = CRS("+proj=longlat +datum=NAD83"))
hull.poly.df <- SpatialPolygonsDataFrame(hull.poly, data=data.frame(ID=1))
hull.grid <- raster(hull.poly, res = 1/100)

# Clean up environmental table
wq.summ <- wq.summ %>% dplyr::select(-contains("max")) # Remove max values
wq.summ <- wq.summ %>% dplyr::select(-contains("min")) # Remove min values
wq.summ <- wq.summ %>% dplyr::select(-contains(".S_")) # Remove observations collected at the surface (only retain observations from bottom of the water column)
wq.summ <- subset(wq.summ, select = -c(APA.B_mean, CHLA.B_mean, APA.B_monthly_range, CHLA.B_monthly_range, pH_monthly_range, APA.B_yearly_range, CHLA.B_yearly_range, pH_yearly_range)) # Remove variables with too much missing data
wq.only <- wq.summ[c(9:82)] # Subset water quality variables only

# Convert sampling data to Spatial object
wq.sp2 <- SpatialPoints(wq.only[,2:1], proj4string = CRS("+proj=utm +datum=NAD83"))
wq.spdf2 <- SpatialPointsDataFrame(wq.sp2, as.data.frame(wq.only))




# Make a raster loop
for(c in names(wq.spdf2)) { 
  interp.var <- c
  test <- na.omit(wq.spdf2@data[,c("LATDEC", "LONDEC", c)])
  #make test into a new spatial dataframe
  wq.sp.new <- SpatialPoints(test[,2:1], proj4string = CRS("+proj=utm +datum=NAD83"))
  wq.spdf.new <- SpatialPointsDataFrame(wq.sp.new, as.data.frame(test))
  
  # Auto-fit a model variogram
  variogram<-automap::autofitVariogram(get(interp.var)~1, wq.spdf.new)

  plot(variogram)
  spatial.points <- SpatialPoints(coords = coordinates(hull.grid), 
                                  proj4string = CRS("+proj=utm +datum=NAD83"))
  spatial.grid <- SpatialPixels(points = spatial.points)
  kriging_result <- automap::autoKrige(get(interp.var)~1, wq.spdf.new, spatial.grid)
  # save variogram
  pdf(file = paste0("variogram_", interp.var, ".pdf"))
  plot(variogram)
  dev.off()
  
  # Create Rasterbrick object and mask to FLKeys polygon
  krige.brick <- brick(kriging_result$krige_output)
  krige.mask <- mask(krige.brick, hull.poly)
  names(krige.mask) <- c('prediction', 'variance', 'stdev')
  # save prediction map
  pdf(file = paste0("prediction_", interp.var, ".pdf"))
  dev.off()
  plot(krige.mask)
  dev.off()
  
  krige.predict <- cbind(as.data.frame(krige.mask$prediction), as.data.frame(spatial.grid))
  ggpredict <- ggplot()+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
    geom_tile(data = krige.predict, aes(x = x, y = y, fill = prediction, col = prediction))+
    geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
    scale_fill_viridis_c(na.value = NA)+
    scale_color_viridis_c(na.value = NA)+
    coord_fixed(xlim = c(-83.3,-79.9), ylim = c(24.2,26), expand = 0)
  # # save ggplot prediction map
  pdf(file = paste0("ggpredict_", interp.var, ".pdf"))
  plot(ggpredict)
  dev.off()
  krig.fixed<-data.frame(cbind(krige.predict$x, krige.predict$y, krige.predict$prediction))
  r = rasterFromXYZ(krig.fixed)
  plot(r)
  writeRaster(r, file = paste0(interp.var, "_raster"), format="raster", overwrite=TRUE)
  
}



