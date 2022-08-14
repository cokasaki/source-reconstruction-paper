library(sf)
library(maps)
library(here)
library(stringr)


states <- c("Alabama","Arizona","Arkansas","California","Colorado","Connecticut","Delaware","Florida","Georgia","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York","North Carolina","North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont","Virginia","Washington","West Virginia","Wisconsin","Wyoming")
states <- tolower(states)
states <- gsub(" ","_",states,fixed=TRUE)
states <- paste(states,"csv",sep=".")
air.2019 <- NULL
for(state in states){
  st.data <- read.csv(here("data","epa_air_quality_data",state),stringsAsFactors=F)
  air.2019 <- rbind(air.2019,st.data)
}
air.jan01 <- subset(air.2019,date="01/01/2019")

air.2019.sf <- st_as_sf(air.2019,coords=c("SITE_LONGITUDE","SITE_LATITUDE"))
air.jan01.sf <- st_as_sf(air.jan01,coords=c("SITE_LONGITUDE","SITE_LATITUDE"))

map('state')

for(i in unique(st_geometry(air.jan01.sf))) {
  plot(i,add=T,pch=19,cex=0.5,col="blue")
}

write.csv(air.2019[,c("Date","Daily.Mean.PM2.5.Concentration","SITE_LATITUDE","SITE_LONGITUDE","STATE","COUNTY")],"data/epa_air_quality_data/processed_air_quality_us.csv")
write.csv(subset(air.2019,Date=="01/01/2019")[,c("Date","Daily.Mean.PM2.5.Concentration","SITE_LATITUDE","SITE_LONGITUDE","STATE","COUNTY")],"data/epa_air_quality_data/jan1_processed_air_quality_us.csv")
write.csv(subset(air.2019,substr(Date,1,2)=="08")[,c("Date","Daily.Mean.PM2.5.Concentration","SITE_LATITUDE","SITE_LONGITUDE","STATE","COUNTY")],"data/epa_air_quality_data/aug_processed_air_quality_us.csv")

for(i in 1:31){
  write.csv(subset(air.2019,Date==paste("08/",str_pad(i,2,pad="0"),"/2019",sep=""))[,c("Date","Daily.Mean.PM2.5.Concentration","SITE_LATITUDE","SITE_LONGITUDE","STATE","COUNTY")],paste("data/epa_air_quality_data/aug",str_pad(i,2,pad="0"),"_processed_air_quality_us.csv",sep=""))
}

