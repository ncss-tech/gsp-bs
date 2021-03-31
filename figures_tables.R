
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

bs1 <- readAll(raster(file.path(fp, "bs1_dsm_mean_mask.tif")))
bs2 <- readAll(raster(file.path(fp, "bs2_dsm_mean_mask.tif")))
bs1 <- projectRaster(bs1, crs = "+init=epsg:5070")
bs2 <- projectRaster(bs2, crs = "+init=epsg:5070")

bs1_pi <- readAll(raster(file.path(fp, "bs1_dsm_mean_mask.tif")))
bs2_pi <- readAll(raster(file.path(fp, "bs2_dsm_mean_mask.tif")))
bs1_pi <- projectRaster(bs1, crs = "+init=epsg:5070")
bs2_pi <- projectRaster(bs2, crs = "+init=epsg:5070")

# bs1
brks <- c(0, 0.2, 0.4, 0.5, 0.8, 1)
tm_bs1 <- tm_shape(bs1, raster.downsample = FALSE) + 
  tm_raster(breaks = brks, title = "Probability", palette = "Greys", ) + 
  tm_legend(legend.position = c("right", "bottom")) +
    # legend.outside = TRUE,
  #           legend.text.size = 0.6, 
  #           legend.title.size = 0.8, 
  #           legend.outside.size = 0.1
  #           ) +
  tm_shape(st) + tm_borders() + 
  tm_layout(frame = FALSE, 
            main.title = "Category 1", 
            main.title.size = 1,
            main.title.position = "center",
            scale = 0.8
            # inner.margins = c()
            )



# bs2
tm_bs2 <- tm_shape(bs2, raster.downsample = FALSE) + 
  tm_raster(breaks = brks, title = "Probability", palette = "Greys") + 
  tm_legend(legend.position = c("right", "bottom")) +
    # legend.outside = TRUE, legend.text.size = 0.6, legend.title.size = 0.8) +
  tm_shape(st) + tm_borders() + 
  tm_layout(frame = FALSE, 
            main.title = "Category 2", 
            main.title.position = "center",
            main.title.size = 1,
            scale = 0.8
            # inner.margins = c()
  )

test <- tmap_arrange(tm_bs1, tm_bs2, ncol = 1, asp = NA)
tmap_save(tm = test, units = "in", width = 6, height = 6, dpi = 300, outer.margins = 0, filename = "bs_test.png")


# BS table ----
brks <- c(0, 0.2, 0.4, 0.5, 0.8, 1)
brks_ <- paste(brks[-6], "to", brks[-1])
acres    <- res(bs1)[1] * res(bs1)[2] * 0.0002471
hectares <- res(bs1)[1] * res(bs1)[2] * 0.0001

bs1_tb <- values(bs1)
bs1_tb <- cut(bs1_tb, breaks = brks, labels = brks_)
bs1_tb <- as.data.frame(t(as.matrix(table(bs1_tb))))
bs2_tb <- values(bs2)
bs2_tb <- cut(bs2_tb, breaks = brks, labels = brks_)
bs2_tb <- as.data.frame(t(as.matrix(table(bs2_tb))))

bs_tb <- rbind(bs1_tb, bs2_tb)
bs_tb

bs_tb_df <- sapply(bs_tb, function(x) {
  # prettyNum(x * acres, big.mark = ",", scientific = FALSE, digits = 0)
  paste0(
    prettyNum(x*hectares/1e6, big.mark = ",", scientific = FALSE, digits = 2), 
    " ", "(", 
    prettyNum(x * acres / 1e6, big.mark = ",", scientific = FALSE, digits = 2),
    ")")
})
test <- as.data.frame(bs_tb_df)
View(test)


# BS percent ----



# CONUS
bs1_val <- values(bs1)
bs2_val <- values(bs2)
conus_bs1 <- sum(bs1_val, na.rm = TRUE) / sum(bs1_val >= 0, na.rm = TRUE) * 100
# sum(bs1_val > 0.5, na.rm = TRUE) / sum(bs1_val >= 0, na.rm = TRUE)
conus_bs2 <- sum(bs2_val, na.rm = TRUE) / sum(bs2_val >= 0, na.rm = TRUE) * 100


# AK
fp <- "D:/geodata/project_data/gsp-bs/GBS_OCONUS"
ak_r <- raster(file.path(fp, "AK_int_gbs_probability_rs_average.tif"))
ak_r <- projectRaster(ak_r, crs = "+init=epsg:3338")
ak_val <- values(ak_r)
AK <- sum(ak_val, na.rm = TRUE) / sum(ak_val >= 0, na.rm = TRUE)


# AS
as_r <- raster(file.path(fp, "AS_int_gbs_probability_rs_average.tif"))
as_val <- values(as_r)
AS <- sum(as_val, na.rm = TRUE) / sum(as_val >= 0, na.rm = TRUE)


# FM
fm_r <- raster(file.path(fp, "FM_int_gbs_probability_rs_average.tif"))
fm_val <- values(fm_r)
FM <- sum(fm_val, na.rm = TRUE) / sum(fm_val >= 0, na.rm = TRUE)


# HI
hi_r <- raster(file.path(fp, "hi_int_gbs_probability_rs_average.tif"))
hi_val <- values(hi_r)
HI <- sum(hi_val, na.rm = TRUE) / sum(hi_val >= 0, na.rm = TRUE)


# MP
mp_r <- raster(file.path(fp, "MP_int_gbs_probability_rs_average.tif"))
mp_val <- values(mp_r)
MP <- sum(mp_val, na.rm = TRUE) / sum(mp_val >= 0, na.rm = TRUE)


# PR
pr_r <- raster(file.path(fp, "PR_int_gbs_probability_rs_average.tif"))
pr_val <- values(pr_r)
PR <- sum(pr_val, na.rm = TRUE) / sum(pr_val >= 0, na.rm = TRUE)


# PW
pw_r <- raster(file.path(fp, "PW_int_gbs_probability_rs_average.tif"))
pw_val <- values(pw_r)
PW <- sum(pw_val, na.rm = TRUE) / sum(pw_val >= 0, na.rm = TRUE)


# vi
vi_r <- raster(file.path(fp, "USVI_int_gbs_probability_rs_average.tif"))
vi_val <- values(vi_r)
VI <- sum(vi_val, na.rm = TRUE) / sum(vi_val >= 0, na.rm = TRUE)


# summary
bs_pct <- round(data.frame(conus_bs1, conus_bs2, AK, AS, FM, HI, MP, PR, PW, VI))
names(bs_pct)[1:2] <- c("CONUS (BS1)", "CONUS (BS2)")
bs_pct



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



# yields ----

library(raster)

yi <- raster("C:/workspace2/corn_yields2020_1000.tif")

fp <- "D:/geodata/project_data/gsp-bs"
bs1 <- readAll(raster(file.path(fp, "bs1_dsm_mean_mask.tif")))
bs2 <- readAll(raster(file.path(fp, "bs2_dsm_mean_mask.tif")))
bs1 <- projectRaster(bs1, crs = "+init=epsg:5070")
bs2 <- projectRaster(bs2, yi, progress = TRUE)

test <- stack(yi, bs2)
test2 <- as.data.frame(values(test))
names(test2) <- c("cy2020", "BS2")
test2 <- subset(test2, complete.cases(cy2020, BS2))
idx <- sample(1:nrow(test2), 10000)



