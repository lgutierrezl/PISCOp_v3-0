rm(list = ls())

source("R/utils/pipe.R")
suppressMessages(library(ggplot2))
suppressMessages(library(tidyterra))
suppressMessages(library(terra))
suppressMessages(library(data.table))

raw_data <- readRDS(
  file.path("data",
            "processed",
            "point",
            "datos_PE_1981-2025_v3.0_03.RDS")
)
dem_SA_data <- rast(
  file.path("data",
            "processed",
            "gridded",
            "topo_vars_sa_gf_aoi.nc")
)[[9]]
shp_SA_data <- terra::vect(
  file.path('data',
            "raw",
            "sa_cnt.gpkg"
  )
) %>% terra::project(., 'EPSG:4326') %>% 
  terra::crop(., ext(-82,-64,-19,2))

shp_eSA_data <- vect(
  file.path("data",
            "raw",
            "sa_eco.gpkg")
) %>% 
  terra::crop(., ext(-82,-64,-19,2))

### spatial data

xyz_shp <- vect(raw_data$xyz,
                geom = c("LON", "LAT"),
                crs = terra::crs(dem_SA_data))

dem_SA_data_masked <- terra::mask(dem_SA_data, shp_SA_data)
dem_SA_data_masked <- terra::aggregate(dem_SA_data_masked, 5)
dem_SA_data_masked <- dem_SA_data_masked / 1000
dem_SA_data_masked[dem_SA_data_masked < 0] <- NA
dem_SA_data_masked <- terra::classify(
  dem_SA_data_masked,
  rcl = c(0, 0.5, 1, 2, 3, 4, 5, Inf),
  include.lowest = TRUE,
  others = NA
)

nstation_per_grid <- terra::rasterize(
  xyz_shp,
  terra::aggregate(dem_SA_data, 100),
  fun = length
)
nstation_per_grid[nstation_per_grid < 1] <- NA
nstation_per_grid <- terra::classify(
  nstation_per_grid,
  rcl = c(1, 2, 5, 10, 20, 50, Inf),
  include.lowest = TRUE,
  others = NA
)


## ecoregions per space

palette_OkabeIto <- c("#56B4E9", "#F0E442", "gray40", "#E69F00", "#009E73",
                      "#0072B2", "#D55E00", "#999999", "#CC79A7")

shp_eSA_data$nr_id <- factor(
  shp_eSA_data$nr_id,
  levels = c("NAS", "PAD", "CAS", "AOL")
)

ecoregions_plot <-
  ggplot() +
  geom_spatvector(data = shp_eSA_data,
                  aes(fill = nr_id),
                  colour = "black",
                  size = .01) +
  geom_spatvector(data = shp_SA_data,
                  fill = NA,
                  size = 0.5,
                  colour = "black") +
  scale_x_continuous(limits = c(-82, -66),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-19, 2),
                     expand = c(0, 0)) +
  scale_fill_manual("       Ecoregions",
                    values = palette_OkabeIto,
                    guide = guide_legend(override.aes = list(color = NA),
                                         title.position = "left",
                                         order = 1,
                                         barheight = .7,
                                         barwidth = .5,
                                         title.theme =
                                           element_text(
                                             size = 9,
                                             angle = 90,
                                             vjust = 0.5
                                           ))) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.text = element_text(size=8.5),
        legend.box.background = element_blank(),
        legend.key.spacing.y = unit(0, "cm"),
        legend.key.size = unit(0, "cm"),
        legend.box = "vertical",
        legend.justification = c(0, 0),
        legend.position = c(0.0, 0.0),
        legend.background = element_blank())


## elevation and countries

countries_plot <-
  ggplot() +
  geom_spatraster(data = dem_SA_data_masked) +
  geom_spatvector(data = shp_SA_data,
                  fill = NA,
                  size = 0.5,
                  colour = "black") +
  scale_fill_manual("Elevation (km asl)",
                    values = terrain.colors(7),
                    na.translate = FALSE,
                    labels = c("0-.5",
                               ".5-1",
                               "1-2",
                               "2-3",
                               "3-4",
                               "4-5",
                               "> 5"),
                    guide = guide_legend(frame.colour = "black",
                                         title.position = "left",
                                         order = 2, # display order
                                         barwidth = .5,
                                         barheight = .7,
                                         title.theme =
                                           element_text(
                                            size = 9,
                                            angle = 90,
                                            vjust = 0.5
                                           ))) +
  scale_x_continuous(limits = c(-82, -66),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-19, 2),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.box = "vertical",
        legend.justification = c(0, 0),
        legend.position = c(0.0, 0.0),
        legend.background = element_blank())

## number of stations per grid (50 km)

n_stations_plot <-
  ggplot() +
  geom_spatraster(data = nstation_per_grid) +
  geom_spatvector(data = shp_SA_data,
                  fill = NA,
                  size = 0.5,
                  colour = "black") +
  scale_fill_viridis_d("  N° of stations",
                       direction = -1,
                       na.translate = FALSE,
                       labels = c("1-2",
                                  "2-5",
                                  "5-10",
                                  "10-20",
                                  "20-50",
                                  "> 50"),
                       guide = guide_legend(frame.colour = "black",
                                            ticks.colour = "black",
                                            title.position = "left",
                                            order = 2, # display order
                                            barheight = .8,
                                            barwidth = .5,
                                            title.theme =
                                              element_text(
                                                size = 9,
                                                angle = 90,
                                                vjust = 0.5
                                              ))) +
  scale_x_continuous(limits = c(-82, -66),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-19, 2),
                     expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.box = "vertical",
        legend.justification = c(0, 0), 
        legend.position = c(0.0, 0.0),
        legend.background = element_blank())


## final plot

final_plot <-
  cowplot::plot_grid(
    countries_plot +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8)) +
      annotate("text", x = -81, y = 1, label = " a)",
               colour = "gray20", size = 4),
    ecoregions_plot +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8)) +
      annotate("text", x = -81, y = 1, label = " b)",
               colour = "gray20", size = 4),
    n_stations_plot +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8)) +
      annotate("text", x = -81, y = 1, label = " c)",
               colour = "gray20", size = 4),
    ncol = 3
  )

ggsave(
  file.path("output", "05_visualization", "study-area-sa.png"),
  plot = final_plot,
  device = "png",
  dpi = 500, scale = 1,
  width = 7, height = 4, units = "in"
)

# number of stations per time per ecoregion

time_step <- seq(as.Date("1981-01-01"), as.Date("2025-12-31"), by = "day")
nstation_per_time <-
  lapply(
    seq_len(nrow(raw_data$data)),
    function(x) {

      table(raw_data$xyz$ECOREGIONS[which(!is.na(raw_data$data[x, ]))]) %>%
        as.data.frame() %>%
        rbind(.,
              data.frame(Var1 = "Total", Freq = sum(.$Freq))) %>%
        data.frame(., time = time_step[x])

    }
  ) %>%
  do.call(rbind, .)

nstation_per_time$Var1 <- factor(
  nstation_per_time$Var1,
  levels = c("NAS",
             "PAD",
             "CAS",
             "AOL",
             "Total"),
  labels = c("NAS",
             "PAD",
             "CAS",
             "AOL",
             "Total")
)

n_stations_per_ecoregions <-
  ggplot() +
  geom_bar(data = nstation_per_time,
           aes(x = time, y = Freq),
           stat = "identity", width = 1, alpha = 1) +
  scale_x_date(breaks = "10 years", date_labels = "%Y",
               limits = c(as.Date("1981-01-01"), as.Date("2025-12-31")),
               expand = c(0, 0)) +
  xlab("") +
  ylab("Number of stations") +
  facet_wrap(~Var1, ncol = 2,
             scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

ggsave(
  file.path("output", "05_visualization", "n-stations-per-ecoregion.png"),
  plot = n_stations_per_ecoregions,
  device = "png",
  dpi = 500, scale = 1,
  width = 6, height = 6, units = "in"
)
