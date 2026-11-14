`%!in%` <- Negate(`%in%`)

download_CHIRPm <- function(date, path, aoi) {
  tmpfile <- tempfile()
  chirp_monthly_path <- sprintf("%sdata/raw/CHIRP-v3.0/CHIRPm/", path)
  server <- "https://data.chc.ucsb.edu/products/CHIRP-v3.0/monthly/latam/tifs/"  #ftp data
  chirp_file <- sprintf("%schirp-v3.0.%s.tif", server, format(date,"%Y.%m"))
  chirp_monthly_file <- paste0(chirp_monthly_path, basename(chirp_file))
  if(!file.exists(chirp_monthly_file)){
    download.file(chirp_file, tmpfile)#, mode = "wget", method = "curl"
    chirp_crop <- crop(rast(tmpfile), aoi)
    writeRaster(chirp_crop, chirp_monthly_file, overwrite = T)
    unlink(tmpfile)
    chirp_crop
  }else{
    chirp_crop <- rast(chirp_monthly_file)
  }
}

download_CHIRPSde <- function(date, path, aoi) {
  chirps_daily_path <- sprintf("%s/data/raw/CHIRPS/CHIRPSd/", path)
  server <- "https://data.chc.ucsb.edu/products/CHIRPS/v3.0/daily/final/rnl/"  #ftp data
  chirp_file <- sprintf("%s/%s/chirps-v3.0.rnl.%s.tif",
                        server,
                        format(date,"%Y"),
                        format(date,"%Y.%m.%d"))
  tmpfile <- paste0(tempfile(),".tif")
  chirp_daily_file <- paste0(chirps_daily_path, gsub('.rnl.','.',basename(chirp_file)))
  if(!file.exists(chirp_daily_file)){
    download.file(chirp_file, tmpfile)#, mode = "wget", method = "curl"
    chirp_crop <- crop(rast(tmpfile), aoi)
    chirp_crop[chirp_crop < 0] = 0
    writeRaster(chirp_crop, chirp_daily_file, overwrite = T)
    unlink(tmpfile)
    chirp_crop
  }else{
    chirp_crop <- rast(chirp_daily_file)
  }
}

create_chirm_d <- function(chirp_d, days, path) {
  # 1. Global params
  e_lumb <- 0.5
  l_month <- as.numeric(format(days, '%m'))[1]
  
  # 2. Read PISCOp_clim and CHIRPMm
  psc_clm <- rast('data/raw/PISCOp_clim2_peb.nc')[[l_month]]
  chp_clm_all <- crop(rast(list.files('data/raw/CHPclim2/', full.names=T)[l_month]), ext(psc_clm))
  chp_clm_all[chp_clm_all<0] <- NA
  chp_clim <- terra::resample(chp_clm_all, psc_clm, method='cubicspline')
  
  CHIRPMds <- list()
  for (index in seq_along(days)) {
    # 3. From 0.05 to 0.1
    chirpx_R <- terra::resample(chirp_d[[index]], psc_clm, method='cubicspline')
    # 4. Create CHIRPMm
    CHIRPMd <- chirpx_R * ((psc_clm + e_lumb)/(chp_clim + e_lumb));  CHIRPMd[CHIRPMd < 0] <- 0
    # 5. Save raster
    writeRaster(CHIRPMd, sprintf("%s/data/raw/CHIRPM/CHIRPMd/CHIRPMd.%s.tif",
                                 path, format(days[index], "%Y.%m.%d")), overwrite=T)
    CHIRPMds[[index]] <- CHIRPMd
  }
  rast(CHIRPMds)
}

create_chirm_m <- function(chirp_m, path, month) {
  # 1. Global params
  e_lumb <- 0.5
  l_month <- as.numeric(format(month, '%m'))
  
  # 2. Read PISCOp_clim and CHIRPMm
  psc_clm <- rast('data/raw/PISCOp_clim2_peb.nc')[[l_month]]
  chp_clm_all <- crop(rast(list.files('data/raw/CHPclim2/', full.names=T)[l_month]), ext(psc_clm))
  chp_clm_all[chp_clm_all<0] <- NA
  chp_clm <- terra::resample(chp_clm_all, psc_clm, method='cubicspline')
  
  # 3. From 0.05 to 0.1
  chirpx_R <- terra::resample(chirp_m, psc_clm, method='cubicspline')
  
  # 4. Create CHIRPMm
  CHIRPMm <- chirpx_R * ((psc_clm + e_lumb)/(chp_clm + e_lumb)); CHIRPMm[CHIRPMm < 0] <- 0
  
  writeRaster(CHIRPMm,sprintf("%s/data/raw/CHIRPM/CHIRPMm/CHIRPMm.%s.tif",path,format(month,"%Y.%m.%d")),overwrite=T)
  CHIRPMm
}

create_chirm_m_tsm <- function(chirp_m, bsn) {
  # 1. Global params
  e_lumb <- 0.5
  
  # 2. Read PISCOp_clim and CHIRPMm
  psc_clm <- rast('data/raw/PISCOp_clim2.nc')
  chp_clm_all <- crop(rast('data/raw/CHPclim2.nc'), ext(psc_clm))
  chp_clm <- terra::resample(chp_clm_all, psc_clm, method='cubicspline')
  
  # 3. series
  psc_clm <- t(terra::extract(psc_clm, bsn, method = 'bilinear', ID = F, fun = mean))
  chp_clm <- t(terra::extract(psd_all, bsn, method = 'bilinear', ID = F, fun = mean))
  chirpx_R <- chirp_m
  
  # 4. Create CHIRPMm
  CHIRPMm <- chirpx_R * ((psc_clm + e_lumb)/(chp_clm + e_lumb)); CHIRPMm[CHIRPMm < 0] <- 0
  
  writeRaster(CHIRPMm,sprintf("%s/data/raw/CHIRPM/CHIRPMm/CHIRPMm.%s.tif",path,format(month,"%Y.%m.%d")),overwrite=T)
  CHIRPMm
}

ROK <- function (gauge, cov, formula, model, crossval = F,nfold, ...) {
  sav_name <- names(gauge)
  names(gauge) <- "gauge"
  sat_cor <- terra::extract(cov,gauge, ID = F) #extraido del sat
  sat_cor[is.na(sat_cor)] = 0
  dTF <- data.frame(obs=gauge$gauge, sat=sat_cor) #creo el df
  gauge$diff <- round(dTF$obs-dTF$sat,1)
  fv <- FitVariogram(formula = diff~1, input_data = as(gauge, 'Spatial'), fix.values = c(0,NA,NA)) #cov_lm <- cov*llm$coefficients[2]+llm$coefficients[1] #corrigo el sat
  cov_lm <- cov
  kd <- krige(diff~1,as(gauge,'Spatial'),as(terra::as.points(cov_lm),'Spatial'),model=fv$var_model)
  rF <- kd[1]
  sp::gridded(rF)=T
  rFF <- rast(rF) + crop(cov_lm, ext(rF))
  rFF <- expm1(rFF)
  rFF[rFF<0] <- 0
  rFF
}

FitVariogram <- function (formula, input_data, model = c("Sph", "Exp", "Gau","Ste"), 
                          kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA,NA, NA), 
                          verbose = FALSE, GLS.model = NA, start_vals = c(NA,NA, NA), 
                          miscFitOptions = list(), boundaries, ...)
{
  if ("alpha" %in% names(list(...)))
    warning("Anisotropic variogram model fitting not supported, see the documentation of autofitVariogram for more details.")
  miscFitOptionsDefaults = list(merge.small.bins = TRUE, min.np.bin = 5)
  miscFitOptions = modifyList(miscFitOptionsDefaults, miscFitOptions)
  longlat = !sp::is.projected(input_data)
  if (is.na(longlat))
    longlat = FALSE
  diagonal = sp::spDists(t(sp::bbox(input_data)), longlat = longlat)[1,
                                                                     2]
  if (!missing(boundaries))
    boundaries = boundaries
  else {
    boundaries = c(2, 4, 6, 9, 12, 15, 25, 35, 50, 65, 80,
                   100) * diagonal * 0.35/100
  }
  if (!is(GLS.model, "variogramModel")) {
    experimental_variogram = variogram(formula, input_data,
                                       boundaries = boundaries, ...)
  }
  else {
    if (verbose)
      cat("Calculating GLS sample variogram\n")
    g = gstat(NULL, "bla", formula, input_data, model = GLS.model,
              set = list(gls = 1))
    experimental_variogram = variogram(g, boundaries = boundaries,
                                       ...)
  }
  if (miscFitOptions[["merge.small.bins"]]) {
    if (verbose)
      cat("Checking if any bins have less than 5 points, merging bins when necessary...\n\n")
    while (TRUE) {
      if (length(experimental_variogram$np[experimental_variogram$np <
                                           miscFitOptions[["min.np.bin"]]]) == 0 | length(boundaries) ==
          1)
        break
      boundaries = boundaries[2:length(boundaries)]
      if (!is(GLS.model, "variogramModel")) {
        experimental_variogram = variogram(formula,
                                           input_data, boundaries = boundaries, ...)
      }
      else {
        experimental_variogram = variogram(g, boundaries = boundaries,
                                           ...)
      }
    }
  }
  if (is.na(start_vals[1])) {
    initial_nugget = min(experimental_variogram$gamma)
  }
  else {
    initial_nugget = start_vals[1]
  }
  if (is.na(start_vals[2])) {
    initial_range = 0.1 * diagonal
  }
  else {
    initial_range = start_vals[2]
  }
  if (is.na(start_vals[3])) {
    initial_sill = mean(c(max(experimental_variogram$gamma),
                          median(experimental_variogram$gamma)))
  }
  else {
    initial_sill = start_vals[3]
  }
  if (!is.na(fix.values[1])) {
    fit_nugget = FALSE
    initial_nugget = fix.values[1]
  }
  else fit_nugget = TRUE
  if (!is.na(fix.values[2])) {
    fit_range = FALSE
    initial_range = fix.values[2]
  }
  else fit_range = TRUE
  if (!is.na(fix.values[3])) {
    fit_sill = FALSE
    initial_sill = fix.values[3]
  }
  else fit_sill = TRUE
  getModel = function(psill, model, range, kappa, nugget,
                      fit_range, fit_sill, fit_nugget, verbose) {
    if (verbose)
      debug.level = 1
    else debug.level = 0
    if (model == "Pow") {
      warning("Using the power model is at your own risk, read the docs of autofitVariogram for more details.")
      if (is.na(start_vals[1]))
        nugget = 0
      if (is.na(start_vals[2]))
        range = 1
      if (is.na(start_vals[3]))
        sill = 1
    }
    obj = try(fit.variogram(experimental_variogram, model = vgm(psill = psill,
                                                                model = model, range = range, nugget = nugget, kappa = kappa),
                            fit.ranges = c(fit_range), fit.sills = c(fit_nugget,
                                                                     fit_sill), debug.level = 0), TRUE)
    if ("try-error" %in% class(obj)) {
      warning("An error has occured during variogram fitting. Used:\n",
              "\tnugget:\t", nugget, "\n\tmodel:\t", model,
              "\n\tpsill:\t", psill, "\n\trange:\t", range,
              "\n\tkappa:\t", ifelse(kappa == 0, NA, kappa),
              "\n  as initial guess. This particular variogram fit is not taken into account. \nGstat error:\n",
              obj)
      return(NULL)
    }
    else return(obj)
  }
  test_models = model
  SSerr_list = c()
  vgm_list = list()
  counter = 1
  for (m in test_models) {
    if (m != "Mat" && m != "Ste") {
      model_fit = getModel(initial_sill - initial_nugget,
                           m, initial_range, kappa = 0, initial_nugget,
                           fit_range, fit_sill, fit_nugget, verbose = verbose)
      if (!is.null(model_fit)) {
        vgm_list[[counter]] = model_fit
        SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))
      }
      counter = counter + 1
    }
    else {
      for (k in kappa) {
        model_fit = getModel(initial_sill - initial_nugget,
                             m, initial_range, k, initial_nugget, fit_range,
                             fit_sill, fit_nugget, verbose = verbose)
        if (!is.null(model_fit)) {
          vgm_list[[counter]] = model_fit
          SSerr_list = c(SSerr_list, attr(model_fit,
                                          "SSErr"))
        }
        counter = counter + 1
      }
    }
  }
  strange_entries = sapply(vgm_list, function(v) any(c(v$psill,
                                                       v$range) < 0) | is.null(v))
  if (any(strange_entries)) {
    if (verbose) {
      print(vgm_list[strange_entries])
      cat("^^^ ABOVE MODELS WERE REMOVED ^^^\n\n")
    }
    warning("Some models where removed for being either NULL or having a negative sill/range/nugget, \n\tset verbose == TRUE for more information")
    SSerr_list = SSerr_list[!strange_entries]
    vgm_list = vgm_list[!strange_entries]
  }
  if (verbose) {
    cat("Selected:\n")
    print(vgm_list[[which.min(SSerr_list)]])
    cat("\nTested models, best first:\n")
    tested = data.frame(`Tested models` = sapply(vgm_list,
                                                 function(x) as.character(x[2, 1])), kappa = sapply(vgm_list,
                                                                                                    function(x) as.character(x[2, 4])), SSerror = SSerr_list)
    tested = tested[order(tested$SSerror), ]
    print(tested)
  }
  result = list(exp_var = experimental_variogram, var_model = vgm_list[[which.min(SSerr_list)]],
                sserr = min(SSerr_list))
  class(result) = c("autofitVariogram", "list")
  return(result)
}

mean_double_Station <- function(gauge = NULL, sat = NULL, longlat = TRUE) {
  point <- terra::as.points(sat)
  
  # 2. Distance from rain gauges to sat points
  distances <- function(x, ptsat = point){which.min(terra::distance(ptsat, gauge[x, ], unit="m"))}
  loc <- do.call("c", lapply(1:length(gauge), distances))
  duplicates <- loc[which(duplicated(loc))]
  
  # gauge (sp) to data.frame (gauge2p)
  gauge2p <- terra::as.data.frame(gauge, geom = 'XY')[,c(ncol(gauge)+1,ncol(gauge)+2,1:ncol(gauge))]
  
  list <- lapply(1:length(duplicates), function(i) {
    dupliStation <- which(loc == duplicates[i])
    gaugeD <- gauge2p[dupliStation, ]
    PromStation <- colMeans(gaugeD, na.rm = T)
    list(Prom = PromStation, position = dupliStation)
  })
  
  stat <- do.call("rbind", lapply(1:length(duplicates), function(x) list[[x]]$position))
  stat2 <- do.call("rbind", lapply(1:length(duplicates), function(x) list[[x]]$Prom))
  gauge2p[stat[, 1], ] <- stat2
  newG <- terra::vect(gauge2p[-stat[, 2], ], geom = c('x', 'y'), crs = 'EPSG:4326')
  return(newG)
}

RIDW <- function (gauge, cov, formula, idpR = seq(0.8, 3.5, 0.1), ...) {
  gauge <- as(gauge,'Spatial')
  cov <- as(cov, "Raster")
  sav_name <- names(gauge)
  ext <- raster::extract(cov, gauge, cellnumber = F, sp = T)
  rg_v <- which(!is.na(ext$sat))
  gauge <- gauge[rg_v,]
  station <- gauge
  linear <- dplyr::mutate_all(dplyr::mutate_all(na.omit(ext@data), as.character), as.numeric)
  station$residuals <- linear[[1]]-linear[[2]]
  point <- data.frame(rasterToPoints(cov))
  coordinates(point) <- ~x + y
  projection(point) <- projection(cov)
  idpRange <- idpR
  mse <- rep(NA, length(idpRange))
  for (i in 1:length(idpRange)) {
    mse[i] <- mean(krige.cv(residuals ~ 1, station, nfold = 10,
                            set = list(idp = idpRange[i]), verbose = F, ...)$residual^2)
  }
  poss <- which(mse %in% min(mse))
  bestparam <- idpRange[poss]
  idwError <- idw(residuals ~ 1, station, point, idp = bestparam)
  idwError <- idwError["var1.pred"]
  gridded(idwError) <- TRUE
  mapa <- raster(idwError)
  Ridw <- mapa+cov
  Ridw[Ridw < 0] <- 0
  terra::rast(Ridw)
}
