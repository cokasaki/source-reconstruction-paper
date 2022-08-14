
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

extract.data <- function(date,run,offset) {
  model.info <- ArchiveGribGrab("gfsanl",date,run,offset,here("data","nomads"))
  
  # available levels are given here-ish
  # https://www.nco.ncep.noaa.gov/pmb/products/gfs/gfs.t00z.pgrb2.1p00.anl.shtml
  levels <- c(seq(1000,925,-25),seq(900,100,-50))
  levels.strs <- sapply(levels,function(x) paste(toString(x),"mb",sep=" "))
  
  # 1.4M nodes for the whole globe
  # only 8k nodes for a grid around WA state
  model.data <- ReadGrib(c(model.info[[1]]$file.name[1]),levels.strs,c("VVEL","UGRD","VGRD"),domain=NULL)
  df <- data.frame(date=model.data$model.run.date,variable=model.data$variables,level=model.data$levels,
                   lon=model.data$lon,lat=model.data$lat,value=model.data$value)
  return(df)
}

extract.daily.data <- function(date,midpoints=FALSE) {
  df.00 <- extract.data(date,"00","+00")
  df.06 <- extract.data(date,"06","+00")
  df.12 <- extract.data(date,"12","+00")
  df.18 <- extract.data(date,"18","+00")
  df <- rbind(df.00,df.06,df.12,df.18)
  
  if(midpoints) {
    df.03 <- extract.data(date,"00","+03")
    df.09 <- extract.data(date,"06","+03")
    df.15 <- extract.data(date,"12","+03")
    df.21 <- extract.data(date,"18","+03")    
    df <- rbind(df,df.03,df.09,df.15,df.21)
  }
  
  return(df)
}

for(i in 1:31){
  df <- extract.daily.data(paste("2019","08",str_pad(i,2,pad="0"),sep=""),TRUE)
  write.csv(df,paste("data/nomads/2019aug",str_pad(i,2,pad="0"),sep=""))
}

