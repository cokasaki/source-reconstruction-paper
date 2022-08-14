##
#
# NOTE WHEN INSTALLING/UPDATING WGRIB2 ON MAC
# NEED TO USE GCC-# NOT GCC, SINCE GCC POINTS TO CLAING
# THEN ONCE BUILT, MOVE grib2/wgrib2/wgrib2 into /usr/local/bin
#
##


library(rNOMADS)
library(stringr)
library(here)
# we want gfs-anl 004 (0.5deg domain)
# this is produced 4/day at 00 06 12 and 1800 hours
# it also produced output timesteps +00 (+03 and +06)
# we will use the +03 if we need interpolation
# since washington is 8x4 degrees roughly we then have
# 128 spatial nodes * 4/day * 365 days = 186880 data points
# however we want a 2 degree buffer so
# 12x8 degrees = 560 640 data points

extract.data <- function(date,run,offset,vars=NULL,levels=NULL,domain=NULL) {
  model.info <- ArchiveGribGrab("gfsanl",date,run,offset,here("data","nomads"))
  
  # 1.4M nodes for the whole globe
  # only 8k nodes for a grid around WA state
  model.data <- ReadGrib(c(model.info[[1]]$file.name[1]),levels,vars,domain=domain)
  df <- data.frame(date=model.data$model.run.date,variable=model.data$variables,level=model.data$levels,
                   lon=model.data$lon,lat=model.data$lat,value=model.data$value)
  return(df)
}

extract.daily.data <- function(date,hours=c("00","06","12","18"),offsets=c("+00","+03"),...) {
  df <- NULL
  for(offset in offsets){
    for(hour in hours){
      df <- rbind(df,extract.data(data,hour,offset,...))
    }
  }
  return(df)
}

extract.monthly.data <- function(year,month,is.leap.year=FALSE,midpoints=FALSE) {
  num.days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  if(is.leap.year) {
    num.days[2] <- 29
  }
  num.days <- num.days[month]
  month.str <- str_pad(month,2,pad="0")
  df <- extract.daily.wa.data(paste(year,month.str,"01",sep=""),midpoints)
  for(i in 2:num.days) {
    day.str <- str_pad(i,2,pad="0")
    new.day <- extract.daily.wa.data(paste(year,month.str,day.str,sep=""),midpoints)
    df <- rbind(df,new.day)
  }
  
  return(df)
}

# nearly 1M obs for the month of march
# mar.01.data <- extract.daily.wa.data("20190301")
# mar.02.data <- extract.daily.wa.data("20190302")
# mar.03.data <- extract.daily.wa.data("20190303")
# mar.04.data <- extract.daily.wa.data("20190304")
# mar.05.data <- extract.daily.wa.data("20190305")
# mar.06.data <- extract.daily.wa.data("20190306")
# mar.07.data <- extract.daily.wa.data("20190307")
# mar.08.data <- extract.daily.wa.data("20190308")
# mar.09.data <- extract.daily.wa.data("20190309")
# mar.10.data <- extract.daily.wa.data("20190310")
# mar.11.data <- extract.daily.wa.data("20190311")
# mar.12.data <- extract.daily.wa.data("20190312")
# mar.13.data <- extract.daily.wa.data("20190313")
# mar.14.data <- extract.daily.wa.data("20190314")
# mar.15.data <- extract.daily.wa.data("20190315")
# mar.16.data <- extract.daily.wa.data("20190316")
# mar.17.data <- extract.daily.wa.data("20190317")
# mar.18.data <- extract.daily.wa.data("20190318")
# mar.19.data <- extract.daily.wa.data("20190319")
# mar.20.data <- extract.daily.wa.data("20190320")
# mar.21.data <- extract.daily.wa.data("20190321")
# mar.22.data <- extract.daily.wa.data("20190322")
# mar.23.data <- extract.daily.wa.data("20190323")
# mar.24.data <- extract.daily.wa.data("20190324")
# mar.25.data <- extract.daily.wa.data("20190325")
# mar.26.data <- extract.daily.wa.data("20190326")
# mar.27.data <- extract.daily.wa.data("20190327")
# mar.28.data <- extract.daily.wa.data("20190328")
# mar.29.data <- extract.daily.wa.data("20190329")
# mar.30.data <- extract.daily.wa.data("20190330")
# mar.31.data <- extract.daily.wa.data("20190331")




# available levels are given here-ish
# https://www.nco.ncep.noaa.gov/pmb/products/gfs/gfs.t00z.pgrb2.1p00.anl.shtml
levels <- c(1000) 
#c(seq(1000,925,-25),seq(900,100,-50)) use this to get vertical structure
levels <- sapply(levels,function(x) paste(toString(x),"mb",sep=" "))


domain <- c(-129,-112,53,41) 
# use this to get western region
# use this to get US

vars <- c("UGRD","VGRD") 
# use this to get vertical c("VVEL","UGRD","VGRD")


jan.01.ground.midnight.wa <- extract.data("20190101","00","+00",vars=c("UGRD","VGRD"),levels="1000 mb",domain=c(-129,-112,54,41) )
jan.01.ground.midnight.west <- extract.data("20190101","00","+00",vars=c("UGRD","VGRD"),levels="1000 mb",domain=c(-129,-106,54,27) )
jan.01.ground.midnight.us <- extract.data("20190101","00","+00",vars=c("UGRD","VGRD"),levels="1000 mb",domain=c(-129,-61,54,20) )
jan.01.ground.midnight.globe <- extract.data("20190101","00","+00",vars=c("UGRD","VGRD"),levels="1000 mb",domain=c(-180,180,90,-90) )

write.csv(jan.01.ground.midnight.wa,"data/jan01windwa.csv")
write.csv(jan.01.ground.midnight.west,"data/jan01windwest.csv")
write.csv(jan.01.ground.midnight.us,"data/jan01windus.csv")
write.csv(jan.01.ground.midnight.globe,"data/jan01windglobe.csv")


# available levels are given here-ish
# https://www.nco.ncep.noaa.gov/pmb/products/gfs/gfs.t00z.pgrb2.1p00.anl.shtml
levels <- c(seq(1000,925,-25),seq(900,100,-50)) #use this to get vertical structure
levels <- sapply(levels,function(x) paste(toString(x),"mb",sep=" "))
vars <- c("VVEL","UGRD","VGRD") # we might also  conceivably want APCP (total precipitation), LAND (land cover proportion), VEG (% vegetation)
aug.01.midnight.us <- extract.data("20190801","00","+00",vars=vars,levels=levels,domain=c(-129,-61,54,20)  )
write.csv(aug.01.midnight.us,"data/nomads/2019aug01_00+00wind_us.csv")

for(i in 2:31){
  aug.i.midnight.us <- extract.data(paste("201908",str_pad(i,2,pad="0"),sep=""),"00","+00",vars=vars,levels=levels,domain=c(-129,-61,54,20) )
  write.csv(aug.i.midnight.us,paste("data/nomads/2019aug",str_pad(i,2,pad="0"),"_00+00wind_us.csv",sep=""))
}




