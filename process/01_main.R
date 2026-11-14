library(xts)
library(terra)
library(RCurl)
library(gstat)
library(raster)
library(cowplot)
library(ggplot2)
library(qmap)
source('src/f_data.R')
#
ini_dat <- "1981-01-01"
end_dat <- "2025-12-31"
sqd_dat <- seq(as.Date(ini_dat), as.Date(end_dat), 'day')
sqm_dat <- seq(as.Date(ini_dat), as.Date(end_dat), 'month')
sqy_dat <- seq(as.Date(ini_dat), as.Date(end_dat), 'year')

aoi_ext <- ext(-82, -64, -19, 2)
ncores <- 6
path <- '/PISCO/PISCOp/'
options(timeout = 2000)
#
## 01: data ----------------------
# Download
for (i in 1:length(sqd_dat)) {
  #download_CHIRPm(date = as.Date(sqm_dat[i], format = "date_%Y%m%d"), path = path, aoi = aoi_ext)
  download_CHIRPSde(date = as.Date(sqd_dat[i], format = "date_%Y%m%d"), path = path, aoi = aoi_ext)
}

# CHIRPM
for (i in 1:length(sqm_dat)) {
  mes_sel <- sqm_dat[i]
  day_sel <- sqd_dat[format(sqd_dat, '%Y%m') == format(mes_sel, '%Y%m')]

  chs_day <- rast(list.files('data/raw/CHIRPS/CHIRPSd/', full.names = T, pattern=format(mes_sel, '%Y.%m')))
  chp_mes <- rast(list.files('data/raw/CHIRP-v3.0/CHIRPm/', full.names = T, pattern=format(mes_sel, '%Y.%m')))

  create_chirm_d(chirp_d = chs_day, days = day_sel, path = path)
  create_chirm_m(chirp_m = chp_mes, month = mes_sel, path = path)
}

# Merge
lst_p4s <- rev(list.files('output/06_database/data/', pattern='*.csv', recursive=T, full.names=T))
psd_xyz <- read.csv('output/06_database/xyz/xyz.csv')

for (j in c(14,18:21)) {
  pcp_sel <- list()
  for (i in 1:length(lst_p4s)) {
    nam_cws <- lst_p4s[i]
    pcp_raw <- read.csv(nam_cws)[c(1,j)]
    pcp_sel[[i]] <- setNames(xts(pcp_raw[,2], order.by=as.Date(pcp_raw[,1])), 
                             substr(nam_cws,30,nchar(nam_cws)-4))
  }
  psd_val <- do.call(cbind, pcp_sel)
  saveRDS(list(xyz = psd_xyz, value = psd_val), 
          paste0('data/processed/',names(read.csv(nam_cws))[j],'.RDS'))
}
#
## 02: selection -----------------
obs_qcr <- readRDS('data/processed/qc_obs.RDS')
obs_nmd <- readRDS('data/processed/obs_mod.RDS')
obs_nbc <- readRDS('data/processed/obs_bc.RDS')
obs_hmd <- readRDS('data/processed/hmg_obs_mod.RDS')
obs_hbc <- readRDS('data/processed/hmg_obs_bc.RDS')

met_sel <- read.csv('data/processed/metodo.csv')

met_sel[met_sel$REV %in% c(0,1),6] <- NA

val_sel <- list()
for (i in 1:nrow(met_sel)) {
  imp_all <- cbind(obs_nmd$value[,obs_nmd$xyz$ID[i]],
                   obs_nbc$value[,obs_nbc$xyz$ID[i]],
                   obs_hmd$value[,obs_hmd$xyz$ID[i]],
                   obs_hbc$value[,obs_hbc$xyz$ID[i]])
  val_sel[[i]] <- imp_all[,met_sel$METODO[i]]
}
val_all <- setNames(do.call(cbind, val_sel), met_sel$ID)

met_sel <- data.frame(subset(met_sel, is.na(met_sel$REV)), row.names = NULL)
met_sel <- met_sel[!is.na(met_sel$METODO),]

imp_obs_met_xyz <- data.frame(subset(obs_qcr$xyz, obs_qcr$xyz$ID %in% met_sel$ID), row.names=NULL)
imp_obs_met_val <- val_all[,imp_obs_met_xyz$ID]

saveRDS(list(xyz = imp_obs_met_xyz, value = imp_obs_met_val),
        'data/processed/imp_obs_met.RDS')
#
## 03: grid ----------------------
stt <- readRDS('data/processed/imp_obs_met.RDS')
stm <- readRDS('data/processed/imp_obm_met.RDS')

parallel::mclapply(1:length(sqm_dat),
                   function(i){
                     mes_sel <- sqm_dat[i]
                     day_sel <- sqd_dat[format(sqd_dat, '%Y%m') == format(mes_sel, '%Y%m')]
                     
                     chirpm_d <- rast(list.files('data/raw/CHIRPM/CHIRPMd/', full.names = T, pattern=format(mes_sel, '%Y.%m')))
                     chirpm_m <- rast(list.files('data/raw/CHIRPM/CHIRPMm/', full.names = T, pattern=format(mes_sel, '%Y.%m')))
                     
                     stt_day <- stt$value[format(mes_sel, '%Y-%m')]
                     stt_mes <- stm$value[mes_sel]
                     
                     sts_day <- crop(vect(data.frame(stt$xyz[,c(3,4)], t(as.matrix(stt_day))), geom=c('LON','LAT'),'EPSG:4326'), ext(chirpm_d))
                     sts_mes <- crop(vect(data.frame(stm$xyz[,c(3,4)], t(as.matrix(stt_mes))), geom=c('LON','LAT'),'EPSG:4326'), ext(chirpm_m))
                     
                     rg_dat_day <- mean_double_Station(gauge = sts_day, sat = chirpm_d[[1]], longlat = T)
                     rg_dat_mes <- mean_double_Station(gauge = sts_mes, sat = chirpm_m, longlat = T)
                     
                     # mes
                     CHIRPMm_log <- log1p(chirpm_m)
                     names(chirpm_m) <- "sat"
                     names(CHIRPMm_log) <- "sat"
                     names(rg_dat_mes) <- "gauge"
                     rg_dat_mes$gauge_l1 <- log1p(rg_dat_mes$gauge)
                     kdm <- ROK(gauge = setNames(rg_dat_mes[,2], 'gauge'), cov = CHIRPMm_log)
                     if(sum(as.vector(kdm) > 10000, na.rm = T) > 0) {
                       kdm <- RIDW(gauge = rg_dat_mes[,1], cov = chirpm_m, formula = "gauge~sat")
                     }
                     writeRaster(kdm, sprintf("%s/data/raw/PISCOp/PISCOp_m/PISCOpm.%s.tif",path,format(mes_sel,"%Y.%m.%d")),overwrite=T)
                     
                     # day
                     kdd <- list()
                     for (j in seq_along(day_sel)) {
                       CHIRPMd_p <- chirpm_d[[j]]
                       rg_dat_dys <- rg_dat_day[,j]
                       names(CHIRPMd_p) <- "sat"
                       names(rg_dat_dys) <- "gauge"
                       kd <- RIDW(gauge = rg_dat_dys, cov = CHIRPMd_p, formula = "gauge~sat")
                       kdd[[j]] <- kd
                       writeRaster(kd, sprintf("%s/data/raw/PISCOp/P_PISCOp_d/PISCOpd.%s.tif", path, format(day_sel[j], "%Y.%m.%d")),overwrite=T)
                     }
                     
                     # day cor
                     l_d_total_select <- rast(kdd)
                     names(l_d_total_select) <- day_sel
                     l_d_total_sum <- sum(l_d_total_select)
                     l_m_total_select <- kdm
                     ratio_C <- (l_m_total_select+1)/(l_d_total_sum+1)
                     PISCOpd <- l_d_total_select*ratio_C
                     mapply(function(u) writeRaster(PISCOpd[[u]],
                                                    sprintf("%s/data/raw/PISCOp/PISCOp_d/PISCOpd.%s.tif",path,format(day_sel[[u]], "%Y.%m.%d")),overwrite=T),
                            1:nlyr(PISCOpd))
                   }, mc.cores = ncores)
 
#
## 04: export ---------------
piscopd_v30_lst <- list.files('data/raw/PISCOp/PISCOp_d/', full.names = T)
writeCDF(rast(piscopd_v30_lst), 
         filename = "data/raw/PISCOp/PISCOp_d.nc",
         varname = "precipitation", unit = 'mm', timename = 'day',
         longname = "PISCOp_v3.0: Gridded rainfall dataset for Peru from 1981 to 2025 (Gutierrez L. et al.; 2026)",
         compression=7, overwrite = T)

piscopm_v30_lst <- list.files('data/raw/PISCOp/PISCOp_m/', full.names = T)
writeCDF(rast(piscopm_v30_lst), 
         filename = "data/raw/PISCOp/PISCOp_m.nc",
         varname = "precipitation", unit = 'mm', timename = 'month',
         longname = "PISCOp_v3.0: Gridded rainfall dataset for Peru from 1981 to 2025 (Gutierrez L. et al.; 2026)",
         compression=7, overwrite = T)

