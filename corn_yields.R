library(USAboundaries)
st <- us_states()
st <- subset(st, ! state_abbr %in% c("AK", "HI", "PR"))


nassqs_auth(key = nass_key)
corn_us <- lapply(st, function(x) {
  cat("getting", x, as.character(Sys.time()), "\n")
  tryCatch({
    corn = nassqs_yields(
      list(commodity_desc = "CORN",
           agg_level_desc = "COUNTY",
           state_alpha = x
      )
    )},
    error = function(err) {
      print(paste("Error occured:  ",err))
      return(NULL)
    }
  )
})
corn_us <- do.call("rbind", corn_us)
# save(corn_us, file = "C:/workspace2/corn_us.RData")
load(file = "C:/workspace2/corn_us.RData")

corn_yield <- subset(corn_us, short_desc == "CORN, GRAIN - YIELD, MEASURED IN BU / ACRE")
corn_yield <- within(corn_yield, {
  Value      = as.numeric(Value)
  year       = as.numeric(year)
  state_name = NULL
  state      = state_alpha 
})
cy <- corn_yield
# cy <- tidyr::pivot_wider(corn_yield, names_from = year, values_from = Value)
cy$StCo_FIPSCode <- paste0(cy$state_fips_code, cy$county_ansi)
# names(cy) <- make.names(names(cy))
cy_st <- sort(unique(subset(cy, year %in% 2016:2020 & Value > 0 & !is.na(Value))$state))


library(dplyr)
library(sf)
fips <- read_sf(dsn = "D:/geodata/government_units/GOVTUNIT_NATIONAL.gdb", layer = "COUNTY")
fips <- rmapshaper::ms_simplify(fips)
fips <- st_transform(fips, crs = 5070)
# save(fips, file = "fips.RData")
load(file = "fips.RData")
fips <- filter(fips, State_FIPSCode %in% st$statefp)

test <- merge(fips, cy, by = "StCo_FIPSCode", all.x = TRUE)
test <- filter(test, Value > 0 & !is.na(Value))
filter(test, year == 2020) %>%
  ggplot() + 
  geom_sf(data = test, aes(fill = Value), col = NA) + 
  scale_fill_viridis_c(na.value = "transparent")


library(terra)
nlcd   <- rast("D:/geodata/land_use_land_cover/NLCD_2016_Land_Cover_L48_20190424.img")
crs(nlcd) <- "epsg:5070"

lapply(2016:2020, function(x) {
  temp   = filter(test, year == x)
  fips_v <- vect(temp)
  rasterize(fips_v, nlcd, field = "Value", filename = paste0("C:/workspace2/yield_", x, ".tif"), wopt = list(datatype = "INT2U"))
})


fp <- untar("2020_cdls.tar.gz", list = TRUE)
# untar("2020_cdls.tar.gz", exdir = "NASS2020", files = fp[grep(".tif$", fp)])
untar("2020_cdls.tar.gz", exdir = "NASS2020", files = fp[1:3])
fr <- list.files("C:/workspace2/NASS2020", full.names = TRUE)
fp <- list.files("C:/workspace2/NASS2020", recursive = TRUE, full.names = TRUE)
file.copy(from = fp, to = gsub("NASS_../", "", fp))

m <- matrix(c(0,   0.9, 0,
              0.9, 1.1, 1,
              1.1, 255, 0
              ),
            ncol = 3,
            byrow = TRUE
            )
yi <- rast("G:/yield_2020.tif")
co <- rast("C:/workspace2/nass2020.tif")
# co <- classify(co, rcl = m, filename = "corn2020.tif", datatype = "INT2U", overwrite = TRUE)

library(gdalUtils)
gdal_setInstallation()
gdalwarp(srcfile = "C:/workspace2/corn2020.tif",
         dstfile = "C:/workspace2/corn20202.tif",
         te = ext(yi)[c(1, 3, 2, 4)],
         r = "near",
         s_srs = crs("+init=epsg:5070"),
         t_srs = crs("+init=epsg:5070"),
         overwrite = TRUE
)
co <- rast("C:/workspace2/corn20202.tif")
test <- co * yi
writeRaster(test, filename = "corn_yields2020.tif", NAflag = 0, overwrite = TRUE)
gdalwarp(srcfile = "C:/workspace2/corn_yields2020.tif",
         dstfile = "C:/workspace2/corn_yields2020_1000.tif",
         tr = c(1000, 1000),
         r = "average",
         overwrite = TRUE
         )




# yields ----

library(raster)
library(ggplot2)

yi <- raster("C:/workspace2/corn_yields2020_1000.tif")

fp <- "D:/geodata/project_data/gsp-bs"
bs1 <- readAll(raster(file.path(fp, "bs1_dsm_mean_mask.tif")))
bs2 <- readAll(raster(file.path(fp, "bs2_dsm_mean_mask.tif")))
bs1 <- projectRaster(bs1, crs = "+init=epsg:5070")
bs2 <- projectRaster(bs2, yi, progress = TRUE)
precip <- readAll(raster("D:/geodata/climate/prism/final_MAP_mm_800m.tif"))
precip <- projectRaster(precip, yi, progress = TRUE)
temp <- readAll(raster("D:/geodata/climate/prism/final_MAAT_800m.tif"))
temp <- projectRaster(temp, yi, progress = TRUE)


test <- stack(yi, bs2, precip, temp)
test2 <- as.data.frame(values(test))
names(test2) <- c("cy2020", "BS2", "precip", "temp")
test2 <- subset(test2, complete.cases(cy2020, BS2))
idx <- sample(1:nrow(test2), 1000)

write.csv(test2, "test_yield.csv", row.names = FALSE)


gg_y_bin <- ggplot(test2, aes(x = cy2020, y = BS2)) + 
  geom_bin2d(show.legend = FALSE) +
  scale_fill_viridis_c() +
  xlab("corn yield") +
  ggtitle("2020 Corn Yield vs BS2")
gg_y_bp <- ggplot(test2, aes(x = cy2020, y = BS2 > 0.5)) + 
  geom_boxplot() +
  xlab("corn yield") +
  ggtitle("2020 Corn Yield vs BS2")
gridExtra::grid.arrange(gg_y_bin, gg_y_bp, ncol = 2)

library(quantreg)
cy_rq <- rq(cy2020 ~ BS2 + precip + temp, data = test2[idx, ])
