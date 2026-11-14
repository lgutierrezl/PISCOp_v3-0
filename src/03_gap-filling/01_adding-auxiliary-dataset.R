rm(list = ls())

suppressMessages(library(xts))
source("R/utils/pipe.R")

raw_data <- readRDS(
  file.path(
    "data",
    "processed",
    "point",
    "datos_PE_1981-2025_v3.0_03.RDS"
  )
)

covs_path <- file.path(file.path("data/processed/gridded/topo_vars_sa_gf_aoi.nc"))

covs <- terra::rast(covs_path)
names(covs) <- c(
  "aspectcosine", "aspectsine", "dist2coast", "dx", "dxx",
  "dy", "dyy", "eastness", "elevation", "latitude", "longitude",
  "northness", "pcurv", "roughness", "slope", "tcurv", "tpi", "tri", "vrm"
)

xy_cells <- terra::extract(
  covs[[1]],
  raw_data$xyz[, c("LON", "LAT")],
  cells = TRUE
)

xy_cells <- xy_cells$cell
xy_topo_cov <- covs[xy_cells]
raw_data$xyz <- cbind(raw_data$xyz, xy_topo_cov)

saveRDS(
  raw_data,
  file.path(
    "data",
    "processed",
    "point",
    "datos_PE_1981-2025_v3.0_03_tc.RDS"
  )
)
