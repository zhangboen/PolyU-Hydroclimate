library(doParallel)
library(foreach)
library(aiRthermo)
library(ncdf4)
library(PCICt)
library(lubridate)
library(RadioSonde)

getTime <- function(nc) {
  origin <- strsplit(nc$dim$time$units,' ')[[1]][3]
  cal <- nc$dim$time$calendar
  origin.pcict <- as.PCICt(origin, cal)
  date_seq <- origin.pcict + nc$dim$time$vals * 86400
  return(date_seq)
}

func <- function(hus_name, ta_name, out_name) {
  # function for calculating CAPE ----
  calc_CMIP5_CAPE <- function(idx) {
    row <- dim0[idx,1]
    col <- dim0[idx,2]
    aPs <- prs
    aTs <- ta[row,col,]
    ahus <- hus[row,col,]
    aws <- q2w(ahus)
    valid_idx <- which(!is.na(aTs))
    capeCin <- CAPE_CIN(precoolType="adiabatic",
                        Ps=aPs[valid_idx], 
                        Ts=aTs[valid_idx], 
                        ws=aws[valid_idx],
                        doLog=0,deltaP=5,
                        upToTop=TRUE)
    return(c(capeCin$cape,capeCin$cin))
  }
  
  # load dataset ----
  nc = nc_open(hus_name)
  hus_mat = ncvar_get(nc, 'hus')
  prs = ncvar_get(nc,'plev')
  nc_close(nc)
  
  nc = nc_open(ta_name)
  ta_mat = ncvar_get(nc, 'ta')

  date_seq <- getTime(nc)
  nc_close(nc)
  
  # cate array to save CAPE and CIN ----
  CAPE_array = array(NA, dim = dim(hus_mat)[c(1,2,4)])
  CIN_array = array(NA, dim = dim(hus_mat)[c(1,2,4)])
  
  # loop calculate for each timestep ----
  for (i in 1:dim(hus_mat)[4]) {
    hus = hus_mat[,,,i]
    ta = ta_mat[,,,i]
    # parallel calculate
    dim0 <- which(!is.na(array(1, c(dim(hus_mat)[1], dim(hus_mat)[2]))), arr.ind = TRUE )
    cl <- makeCluster( 12 )
    registerDoParallel(cl)
    CAPE_CIN_num = foreach(idx=1:nrow(dim0), .combine='rbind', .packages = c("aiRthermo")) %dopar% calc_CMIP5_CAPE(idx)
    stopCluster(cl)
    # transform output to array
    for (idx in 1:nrow(dim0)) {
      row <- dim0[idx,1]
      col <- dim0[idx,2]
      CAPE_array[row,col,i] <- CAPE_CIN_num[idx,1]
      CIN_array[row,col,i] <- CAPE_CIN_num[idx,2]
    }
    print(paste('Finish calculating for month', date_seq[i]))
  }
  
  # define dimension ----
  nc = nc_open(hus_name)
  londim <- nc$dim$lon
  latdim <- nc$dim$lat
  timedim <- nc$dim$time
  nc_close(nc)
  # define variables ----
  fillvalue <- 1e32
  dlname <- "Convective Available Potential Energy"
  cape_def <- ncvar_def("cape", "J/kg", list(londim,latdim,timedim), fillvalue, dlname)
  dlname <- "Convective inhibition"
  cin_def <- ncvar_def('cin', "J/kg", list(londim,latdim,timedim), fillvalue, dlname)
  # create netCDF file to save history ----
  ncout <- nc_create(out_name, list(cape_def, cin_def), force_v4=T)
  ncvar_put(ncout, "cape", CAPE_array)
  ncvar_put(ncout, 'cin', CIN_array)
  nc_close(ncout)
}

setwd('e:/temp')
hus_name = 'hus_Amon_HadGEM2-ES_historical_r1i1p1_198412-200511.nc'
ta_name = 'ta_Amon_HadGEM2-ES_historical_r1i1p1_198412-200511.nc'
out_name = 'cape_cin_Amon_HadGEM2-ES_historical_r1i1p1_198412-200511.nc'
func(hus_name, ta_name, out_name)
