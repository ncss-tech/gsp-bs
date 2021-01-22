library(aqp)
library(soilDB)
library(dplyr)


# fetch SSURGO ----

f_us <- fetchGDB(dsn ="C:/geodata/soils/gNATSGO_CONUS_FY20.gdb", WHERE = "areasymbol LIKE '%'")
mu_us <- get_mapunit_from_GDB(dsn = "C:/geodata/soils/gNATSGO_CONUS_FY20.gdb", stats = TRUE)
# save(f_us, mu_us, file = "C:/Users/stephen.roecker/Nextcloud/data/gnatsgo_fy20.RData")
load(file = "C:/Users/stephen.roecker/Nextcloud/projects/2020_gsp-sas/gnatsgo_fy20.RData")

f_statsgo <- fetchGDB(dsn ="C:/geodata/soils/gNATSGO_CONUS_FY20.gdb", WHERE = "areasymbol = 'US'")

f_us <- f_us[! f_us$mukey %in% f_statsgo$mukey]
f_us <- aqp::union(list(f_us, f_statsgo))

missing_mukeys <- mu_us[is.na(mu_us$n_component), "mukey"]

# missing mukeys
f_sda <- fetchSDA(WHERE = paste0("mukey IN ('", paste0(missing_mukeys, collapse = "', '"), "')"), duplicates = TRUE)
s <- site(f_sda)[siteNames(f_us)]
h <- horizons(f_sda)[horizonNames(f_us)]
f_mis <- h
depths(f_mis) <- cokey ~ hzdept_r + hzdepb_r
site(f_mis) <- s

f_us <- aqp::union(list(f_us, f_mis))



# segment and tranform SSURGO ----

f_sub <- segment(f_us, intervals = c(0, 25))

vars <- c("cokey", "chkey", "hzdept_r", "hzdepb_r", "om_r", "cec7_r", "sumbases_r")

h <- horizons(f_sub)[vars]
names(h) <- gsub("_r", "", names(h))

h <- merge(h, site(f_sub)[c("cokey", "taxsubgrp")], by = "cokey", all.x = TRUE)

h <- transform(h, 
               BS = sumbases / cec7 * 100,
               OC = om / 1.724,
               m_chroma = ifelse(grepl("olls$|^mollic|moll", taxsubgrp), 3, 5),
               m_value  = ifelse(grepl("olls$|^mollic|moll", taxsubgrp), 3, 5),
               d_value  = ifelse(grepl("olls$|^mollic|moll", taxsubgrp), 5, 6)
               )

h2 <- allocate(object = h, hztop = "hzdept", hzbot = "hzdepb", pedonid = "cokey", BS = "BS", OC = "OC", m_chroma = "m_chroma", m_value = "m_value", d_value = "d_value", to = "FAO Black Soil")


s <- site(f_us)
s <- merge(s, h2, by = "cokey", all.x = TRUE)

s2 <- s %>%
  group_by(mukey, cokey, BS2) %>%
  summarize(pct = comppct_r) %>%
  ungroup() %>%
  arrange(mukey, -pct, -BS2) %>%
  filter(!duplicated(mukey)) %>%
  select(- cokey)

write.csv(s2, file = "ssurgo_black_soils.csv", row.names = FALSE)



# rasterize ----

library(raster)


s2 <- within(s2, {
  # mukey   = as.integer(mukey)
  ID   = as.integer(ID)
  bs      = as.integer(BS2)
})
names(s2)[1] <- "ID"

gnatsgo <- "D:/geodata/soils/gnatsgo_fy20_30m.tif"
r <- raster(gnatsgo)
gdalUtils::gdalwarp(
  srcfile = gnatsgo,
  dstfile = gsub("30m.tif", "120m.tif", gnatsgo),
  t_srs = proj4string(r), 
  te = bbox(r),
  r = "mode",
  tr = c(120, 120),
  ot = "Int32",
  verbose   = TRUE,
  overwrite = TRUE
)

r <- raster("D:/geodata/soils/gnatsgo_fy20_1km.tif")
# r2 <- r[1:100, 1:100, drop = FALSE]
# r2 <- ratify(r2)
# rat <- levels(r2)[[1]]
levels(r) <- as.data.frame(s2)

vars <- c("bs")[1]
lapply(vars, function(x) {
  cat(x, as.character(Sys.time()), "\n")
  # beginCluster(type = "SOCK")
  deratify(r, att = x, 
           filename = paste0("D:/geodata/soils/gnatsgo_fy20_1km_", x, ".tif"),
           options = c("COMPRESS=DEFLATE"), 
           overwrite = TRUE, progress = "text" 
  )
  # endCluster()
  cat(x, as.character(Sys.time()), "\n")
})


