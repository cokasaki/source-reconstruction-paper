library(sf)
library(maps)

# Load the air quality data
wa.2019 <- read.csv("data/epa_air_quality_data/epa_air_quality_data_wa_2019.csv",stringsAsFactors=F)
or.2019 <- read.csv("data/epa_air_quality_data/epa_air_quality_data_or_2019.csv",stringsAsFactors=F)
ca.2019 <- read.csv("data/epa_air_quality_data/epa_air_quality_data_ca_2019.csv",stringsAsFactors=F)
idaho.2019 <- read.csv("data/epa_air_quality_data/epa_air_quality_data_idaho_2019.csv",stringsAsFactors=F)
nv.2019 <- read.csv("data/epa_air_quality_data/epa_air_quality_data_nv_2019.csv",stringsAsFactors=F)

air.2019 <- rbind(wa.2019,or.2019,ca.2019,idaho.2019,nv.2019)
air.2019.sf <- st_as_sf(air.2019,coords=c("SITE_LONGITUDE","SITE_LATITUDE"))

map('state',regions=c('Washington','Oregon','California','Idaho','Nevada'))

for(i in unique(st_geometry(air.2019.sf))) {
  plot(i,add=T,pch=19,cex=0.5,col="blue")
}

write.csv(air.2019[,c("Date","Daily.Mean.PM2.5.Concentration","SITE_LATITUDE","SITE_LONGITUDE","STATE","COUNTY")],"data/processed_air_quality_west.csv")
write.csv(subset(air.2019,Date=="01/01/2019")[,c("Date","Daily.Mean.PM2.5.Concentration","SITE_LATITUDE","SITE_LONGITUDE","STATE","COUNTY")],"data/jan1_processed_air_quality_west.csv")

