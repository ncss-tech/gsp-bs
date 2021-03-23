
setwd("D:/geodata/project_data/gsp-bs/data")


# pedon map ----

library(sf)
library(ggplot2)
library(USAboundaries)

st <- us_states() 
st <- subset(st, !state_name %in% c("Alaska", "Hawaii", "Puerto Rico"))
st <- st_transform(st, crs = 5070)

bs <- read_sf(dsn = "ldm_bs3.shp")
bs <- st_transform(bs, crs = 5070)

bs1 <- subset(bs, BS == "BS1")
bs2 <- subset(bs, BS == "BS2")
bs0 <- subset(bs, BS == "None")


gg_bs <- ggplot() +
  geom_sf(data = st, fill = NA) +
  geom_sf(data = bs0, aes(col = BS), size = 0.2) +
  geom_sf(data = bs2, aes(col = BS), size = 0.2) +
  geom_sf(data = bs1, aes(col = BS), size = 0.2) +
  scale_color_manual(values = c("Yellow", "Orange", "Navy")) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(gg_bs, file = "pedon_location.png", dpi = 350)


# BS map ----

library(tmap)
library(raster)

fp <- "D:/geodata/project_data/gsp-bs"

bs1 <- raster(file.path(fp, "bs1_dsm_mean_mask.tif"))
bs2 <- raster(file.path(fp, "bs2_dsm_mean_mask.tif"))

tm_shape(bs1) + tm_raster(breaks = 0.5)


# soil profiles ----

data <- readRDS("tdata_bs_geo.rds")
load(file = "fselection_bs.RData")

data$rn <- row.names(data)
bs1 <- na.exclude(data[c("BS1", "rn", as.character(bs1_fs$variable[1:25]))])
bs2 <- na.exclude(data[c("BS2", "rn", as.character(bs2_fs$variable[1:25]))])

# BS1
table(bs1$BS1)
idx <- bs1$rn %in% data$rn[data$rcasiteid %in% raca_h$rcasiteid]
sum(idx)  # RaCA
sum(!idx) # SCD

# SCD
table(bs2$BS2)
idx <- bs$rn %in% data$rn[data$rcasiteid %in% raca_h$rcasiteid]
sum(idx)
sum(!idx)



# methods ----

load(file = "D:/geodata/project_data/gsp-bs/data/LDM-compact_20200709.RData")
lp <- ldm_bs; rm(ldm_bs)

table(ldm$chemical$total_carbon_ncs_method)
table(ldm$chemical$cec_nh4_ph_7_method)
table(ldm$chemical$base_sat_nh4oac_ph_7)


