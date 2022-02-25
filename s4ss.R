
library(sf)
library(raster)
# library(terra)

setwd("D:/geodata/project_data/gsp-bs/data")


bs  <- readRDS(file = "tdata_bs.rds")
bs2 <- subset(bs, is.na(rcasiteid) & !is.na(pedlabsampnum))
bs2$rcasiteid <- NULL
names(bs2)[names(bs2) %in% c("BS1", "BS2")] <- c("BS1_obs", "BS2_obs")

vars <- c(dsm_bs1 = "bs1_dsm_mean.tif",
          dsm_bs2 = "bs2_dsm_mean.tif"
)
fp <- "D:/geodata/project_data/gsp-bs"
rs_bs <- stack(file.path(fp, vars))
names(rs_bs) <- c("BS1_pred", "BS2_pred")

pred <- extract(rs_bs, st_coordinates(bs2))
bs3 <- cbind(st_coordinates(bs2), bs2, pred)
bs3 <- as.data.frame(st_drop_geometry(bs3))

write.csv(bs3, file = "C:/workspace2/github/ncss-tech/stats_for_soil_survey/data/gsp_bs.csv", row.names = FALSE)
