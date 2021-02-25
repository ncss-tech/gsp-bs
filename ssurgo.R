library(aqp)
library(soilDB)
library(dplyr)


# fetch SSURGO ----

f_us <- fetchGDB(dsn ="C:/geodata/soils/gNATSGO_CONUS_FY20.gdb", WHERE = "areasymbol LIKE '%'")
mu_us <- get_mapunit_from_GDB(dsn = "C:/geodata/soils/gNATSGO_CONUS_FY20.gdb", stats = TRUE)
# save(f_us, mu_us, file = "C:/Users/stephen.roecker/Nextcloud/data/gnatsgo_fy20.RData")
load(file = "C:/Users/stephen.roecker/Nextcloud/projects/2020_gsp-sas/gnatsgo_fy20.RData")

f_statsgo <- fetchGDB(dsn ="D:/geodata/soils/gNATSGO_CONUS_FY20.gdb", WHERE = "areasymbol = 'US'")

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

f_us <- aqp::combine(list(f_us, f_mis))


# duplicated components
test <- aggregate(comppct_r ~ mukey, data = site(f_us), function(x) sum(x, na.rm = TRUE))
idx  <- subset(test, comppct_r > 100)$mukey

View(site(f_us)[site(f_us)$mukey %in% idx, ])
View(mu_us[mu_us$mukey %in% idx, ])

idx <- with(site(f_us), !duplicated(paste(mukey, compname, comppct_r, compkind, localphase, drainagecl, erocl, slope_r)))
f_us2 <- f_us[which(idx), ]


# save(f_us, file = "C:/Users/stephen.roecker/OneDrive - USDA/f_us.RData")
load(file = "C:/Users/stephen.roecker/OneDrive - USDA/f_us.RData")



# soil series colors ----
ss <- get_soilseries_from_NASIS()

ss$soilseriesname <- tolower(gsub("^ | $", "", ss$soilseriesname))
n <- length(ss$soilseriesname)
idx <- data.frame(n = 1:n, h = floor(round((1:n) / 10)) * 10)

osd <- {
  split(idx, idx$h) ->.;
  lapply(., function(x) {
    cat("getting", unique(x$h))
    temp = fetchOSD(ss$soilseriesname[x$n])
    cat(" ", length(temp), " \n")
    return(temp)
  }) ->.;
}
# saveRDS(osd, "C:/Users/stephen.roecker/OneDrive - USDA/data/osd.rds")
osd <- readRDS("C:/Users/stephen.roecker/OneDrive - USDA/data/osd.rds")


# tidy
osd2 <- lapply(osd, function(x) {
  if (!is.null(x)) {
    h <- horizons(x)
    h$hzID <- NULL
    s <- site(x)
  }
  return(list(s = s, h = h))
})

idx <- which(! sapply(osd2, function(x) ncol(x$s)) > 32)
osd2 <- osd2[-185]
s <- do.call("rbind", lapply(osd2, function(x) x$s))
h <- do.call("rbind", lapply(osd2, function(x) x$h))

osd_spc <- h
depths(osd_spc) <- id ~ top + bottom
site(osd_spc) <- s


# segment
osd_seg <- segment(osd_spc, intervals = c(0, 25))
idx_max <- aggregate(bottom ~ id, 
                     data = horizons(osd_seg), 
                     FUN = function(x) {
                       as.integer(max(x, na.rm = TRUE)) == 25
                     }, 
                     na.action = na.pass
)
sum(idx_max$bottom, na.rm = TRUE)
idx <- which(site(osd_seg)$id %in% idx_max$id)
osd_seg <- osd_seg[idx, ]

vars <- c("hue", "value", "chroma", "dry_hue", "dry_value", "dry_chroma")
new_nm <- c("m_hue", "m_value", "m_chroma", "d_hue", "d_value", "d_chroma")
horizonNames(osd_seg)[horizonNames(osd_seg) %in% vars] <- new_nm


# check colors
osd_bol <- horizons(osd_seg) %>%
  mutate(
    chroma = m_chroma <= 3,
    value  = m_value  <= 3 & d_value <= 5,
    check  = chroma & value
  ) %>%
  group_by(id) %>%
  summarize(
    check  = all(check)
  ) %>%
  ungroup() %>%
  as.data.frame()


# convert to boolean
osd_bol <- osd_bol %>%
  filter(check == TRUE) %>%
  mutate(
    id = tolower(id),
    m_chroma = 3,
    m_value  = 3,
    d_value  = 5
  )

# write.csv(osd_bol, file = "C:/Users/stephen.roecker/OneDrive - USDA/data/osd_bol.csv", row.names = FALSE)



# segment and transform SSURGO ----

osd_bol$id <- tolower(osd_bol$id)

f_sub <- segment(f_us, intervals = c(0, 25))
f_sub$compname <- tolower(f_sub$compname)

vars <- c("cokey", "chkey", "hzdept_r", "hzdepb_r", "om_r", "cec7_r", "sumbases_r")
h <- horizons(f_sub)[vars]
names(h) <- gsub("_r", "", names(h))

h <- merge(h, site(f_sub)[c("cokey", "compname")], by = "cokey", all.x = TRUE)
h <- merge(h, osd_bol, by.x = "compname", by.y = "id", all.x = TRUE)

h <- transform(h,
               BS = sumbases / cec7 * 100,
               OC = om / 1.724
               # m_chroma = ifelse(grepl("olls$|^mollic|moll", taxsubgrp), 3, 5),
               # m_value  = ifelse(grepl("olls$|^mollic|moll", taxsubgrp), 3, 5),
               # d_value  = ifelse(grepl("olls$|^mollic|moll", taxsubgrp), 5, 6)
               )

h2 <- allocate(object = h, hztop = "hzdept", hzbot = "hzdepb", pedonid = "cokey", BS = "BS", OC = "OC", m_chroma = "m_chroma", m_value = "m_value", d_value = "d_value", to = "FAO Black Soil")


s <- site(f_us)
s <- merge(s, h2, by = "cokey", all.x = TRUE)


# compute dominant condition for black soils by mukey
s2 <- s %>%
  group_by(mukey) %>%
  summarize(pct_bs = sum(comppct_r[BS2 == TRUE], na.rm = TRUE) / 100) %>%
  ungroup() %>%
  as.data.frame()


# write.csv(s2, file = "ssurgo_black_soils_v2.csv", row.names = FALSE)
s2 <- read.csv(file = "ssurgo_black_soils_v2.csv", stringsAsFactors = FALSE)



# rasterize ----

library(raster)


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

vars <- c("pct_bs")[1]
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


