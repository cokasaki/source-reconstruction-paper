
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

extract.wa.data <- function(date,run,offset) {
  model.info <- ArchiveGribGrab("gfsanl",date,run,offset,here("data","nomads"))
  
  # available levels are given here-ish
  # https://www.nco.ncep.noaa.gov/pmb/products/gfs/gfs.t00z.pgrb2.1p00.anl.shtml
  levels <- c(seq(1000,925,-25),seq(900,100,-50))
  levels.strs <- sapply(levels,function(x) paste(toString(x),"mb",sep=" "))
  
  # 1.4M nodes for the whole globe
  # only 8k nodes for a grid around WA state
  model.data <- ReadGrib(c(model.info[[1]]$file.name[1]),levels.strs,c("VVEL","UGRD","VGRD"),domain=c(-129,-112,53,41))
  df <- data.frame(date=model.data$model.run.date,variable=model.data$variables,level=model.data$levels,
                   lon=model.data$lon,lat=model.data$lat,value=model.data$value)
  return(df)
}

extract.daily.wa.data <- function(date,midpoints=FALSE) {
  df.00 <- extract.wa.data(date,"00","+00")
  df.06 <- extract.wa.data(date,"06","+00")
  df.12 <- extract.wa.data(date,"12","+00")
  df.18 <- extract.wa.data(date,"18","+00")
  df <- rbind(df.00,df.06,df.12,df.18)
  
  if(midpoints) {
    df.03 <- extract.wa.data(date,"00","+03")
    df.09 <- extract.wa.data(date,"06","+03")
    df.15 <- extract.wa.data(date,"12","+03")
    df.21 <- extract.wa.data(date,"18","+03")    
    df <- rbind(df,df.03,df.09,df.15,df.21)
  }
  
  return(df)
}

extract.monthly.wa.data <- function(year,month,is.leap.year=FALSE,midpoints=FALSE) {
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
jan.01.data <- extract.daily.wa.data("20190101")
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

jan.01.ground.midnight <- subset(jan.01.data,variable!="VVEL" & level=="1000 mb" & date=="2019-01-01 00:00:00")
write.csv(jan.01.ground.midnight,"data/jan01wind.csv")

