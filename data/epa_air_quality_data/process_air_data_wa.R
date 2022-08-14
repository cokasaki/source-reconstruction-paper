library(sf)
library(maps)

# Load the air quality data
wa.2019 <- read.csv("data/epa_air_quality_data/epa_air_quality_data_wa_2019.csv",stringsAsFactors=F)

wa.2019.sf <- st_as_sf(wa.2019,coords=c("SITE_LONGITUDE","SITE_LATITUDE"))

map('state',regions=c('Washington'))

for(i in unique(st_geometry(wa.2019.sf))) {
  plot(i,add=T,pch=19,cex=0.5,col="blue")
}

write.csv(wa.2019[,c("Date","Daily.Mean.PM2.5.Concentration","SITE_LATITUDE","SITE_LONGITUDE","COUNTY","STATE")],"data/processed_air_quality_wa.csv")
write.csv(subset(wa.2019,Date=="01/01/2019")[,c("Date","Daily.Mean.PM2.5.Concentration","SITE_LATITUDE","SITE_LONGITUDE","COUNTY","STATE")],"data/jan1_processed_air_quality_wa.csv")

