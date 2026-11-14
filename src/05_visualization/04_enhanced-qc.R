rm(list = ls())

source("R/utils/pipe.R")
suppressMessages(library(xts))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(tidyterra))
suppressMessages(library(terra))
aoi_ext <- ext(-82, -64, -19, 2)

shp_esa_data <- vect(
  file.path("data",
            "raw",
            "sa_cnt.gpkg")
)

# qc output data

qc <- readRDS(
  file.path(
    "output", "02_quality-control",
    "eqc-01-summary.RDS"
  )
)

qc$lenght_nona <- 1

ndata <- aggregate(lenght_nona ~ ECOREGIONS, qc, sum)

ndata_trunc <- aggregate(
  truncation ~ ECOREGIONS,
  qc,
  function(x) length(x[x > 1])
)
ndata_trunc$truncation_percent <- round(
  100 * ndata_trunc$truncation / ndata$lenght_nona, 3
)

ndata_smalgaps <- aggregate(
  small_gaps ~ ECOREGIONS,
  qc,
  function(x) length(x[x > 1])
)
ndata_smalgaps$small_gaps_percent <- round(
  100 * ndata_smalgaps$small_gaps / ndata$lenght_nona, 3
)

ndata_weekly_cycles <- aggregate(
  weekly_cycles ~ ECOREGIONS,
  qc,
  function(x) length(x[x > 1])
)
ndata_weekly_cycles$weekly_cycles_percent <- round(
  100 * ndata_weekly_cycles$weekly_cycles / ndata$lenght_nona, 3
)

ndata_precision_rounding <- aggregate(
  precision_rounding ~ ECOREGIONS,
  qc,
  function(x) length(x[x > 1])
)
ndata_precision_rounding$precision_rounding_percent <- round(
  100 * ndata_precision_rounding$precision_rounding / ndata$lenght_nona, 3
)

qc_output_size <-
list(ndata_trunc,
     ndata_smalgaps,
     ndata_weekly_cycles,
     ndata_precision_rounding) %>%
  lapply(function(x){

    out <- transpose(x)
    colnames(out) <- out[1, ]
    out <- out[-1, ]
    rownames(out) <- colnames(x)[-1]
    out

  })

qc_output_size <- do.call(rbind, qc_output_size)
qc_output_size <- qc_output_size[,
  c("NAS", "PAD", "CAS", "AOL")
]

write.csv(
  qc_output_size,
  file.path("output", "05_visualization", "eqc_output_size_level2.csv")
)


enh_qc_df <-
qc[, c("LON", "LAT", "truncation", "small_gaps", "weekly_cycles", "precision_rounding")] %>%
  reshape2::melt(
    measure.vars = c("truncation", "small_gaps", "weekly_cycles", "precision_rounding"),
    value.name = "Level"
  )
enh_qc_df$Level <- factor(enh_qc_df$Level)
enh_qc_df$variable <- factor(
  enh_qc_df$variable,
  levels = c("truncation", "small_gaps", "weekly_cycles", "precision_rounding"),
  labels = c("Truncation", "Small gaps", "Weekly cycles", "Precision and rounding")
)

enh_qc_plots <-
  ggplot() +
  geom_spatvector(data = crop(shp_esa_data, aoi_ext),
                  colour = "black",
                  fill = NA,
                  size = .1) +
  geom_point(data = enh_qc_df,
             aes(x = LON,
                 y = LAT,
                 colour = Level,
                 fill = Level,
                 size = Level),
             shape = 21) +
  scale_fill_manual(values = c("gray70", "cornflowerblue", "red2")) +
  scale_colour_manual(values = c("gray70", "cornflowerblue", "red2")) +
  scale_size_manual(values = c("0" = .25, "1" = .5, "2" = 1)) +
  scale_x_continuous(limits = c(-82, -64),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-19, 2),
                     expand = c(0, 0)) +
  facet_wrap(~ variable, ncol = 4) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.text = element_text(size = 8.5),
        legend.box.background = element_blank(),
        legend.key.spacing.y = unit(0, "cm"),
        legend.key.size = unit(0, "cm"),
        legend.box = "vertical",
        legend.justification = c(0, 0),
        legend.position = c(0.375, 0.1),
        legend.background = element_blank(),
        legend.title.position = "left",
        legend.title = element_text(size = 10, angle = 90, hjust = 0.5)) +
  theme(strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

ggsave(
  file.path("output", "05_visualization", "eqc-output-type.png"),
  plot = enh_qc_plots,
  device = "png",
  dpi = 500, scale = 1,
  width = 7, height = 3, units = "in"
)


#
qcf <- readRDS("data/processed/point/datos_PE_1981-2025_v3.0_03.RDS")$xyz
best_stations <- plyr::join(qc, qcf[,c(1,10)], by = 'ID')
names(best_stations)[16] <- 'EQCF'
best_stations$EQC <- ifelse(is.na(best_stations$EQCF), best_stations$EQC, best_stations$EQCF)
best_stations <- best_stations[, c("LON", "LAT", "EQC")]
best_stations$EQC <- factor(best_stations$EQC,
                            levels = c(0, 1),
                            labels = c("No Selected",
                                       "Selected"))
best_stations_plot <-
  ggplot() +
  geom_spatvector(data = crop(shp_esa_data, aoi_ext),
                  colour = "black",
                  fill = NA,
                  size = .1) +
  geom_point(data = best_stations,
             aes(x = LON,
                 y = LAT,
                 colour = EQC,
                 fill = EQC,
                 size = EQC),
             shape = 21) +
  scale_fill_manual(values = c("red2", "cornflowerblue")) +
  scale_colour_manual(values = c("red2", "cornflowerblue")) +
  scale_size_manual(values = c("No Selected" = .25, "Selected" = .25)) +
  scale_x_continuous(limits = c(-82, -64),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-19, 2),
                     expand = c(0, 0)) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.text = element_text(size = 8.5),
        legend.box.background = element_blank(),
        legend.key.spacing.y = unit(0, "cm"),
        legend.key.size = unit(0, "cm"),
        legend.box = "vertical",
        legend.justification = c(0, 0),
        legend.position = c(0.05, 0.05),
        legend.background = element_blank(),
        legend.title = element_text(size = 10)) +
  labs(fill = "Enhanced QC",
       colour = "Enhanced QC",
       size = "Enhanced QC") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

ggsave(
  file.path("output", "05_visualization", "eqc-output-selected.png"),
  plot = best_stations_plot,
  device = "png",
  dpi = 600, scale = 1,
  width = 80, height = 100, units = "mm"
)
