library(aqp)
library(soilDB)
library(dplyr)

setwd("D:/geodata/project_data/gsp-bs/data")


# fetch SSURGO ----

f_us <- fetchGDB(dsn ="D:/geodata/soils/gNATSGO_CONUS_FY20.gdb", WHERE = "areasymbol LIKE '%'")
mu_us <- get_mapunit_from_GDB(dsn = "D:/geodata/soils/gNATSGO_CONUS_FY20.gdb", stats = TRUE)
# save(f_us, mu_us, file = "C:/Users/stephen.roecker/Nextcloud/data/gnatsgo_fy20.RData")
load(file = "C:/Users/stephen.roecker/Nextcloud/projects/2020_gsp-sas/gnatsgo_fy20.RData")

f_statsgo <- fetchGDB(dsn ="D:/geodata/soils/gNATSGO_CONUS_FY20.gdb", WHERE = "areasymbol = 'US'")
f_statsgo$mukey <- paste0("STASTGO", f_statsgo$mukey)
f_statsgo$cokey <- paste0("STASTGO", f_statsgo$cokey)

f_us <- aqp::combine(list(f_us, f_statsgo))

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
test <- aggregate(comppct_r ~ mukey, data = co_us, function(x) sum(x, na.rm = TRUE))
idx  <- subset(test, comppct_r > 100)$mukey

View(site(f_us)[site(f_us)$mukey %in% idx, ])
View(mu_us[mu_us$mukey %in% idx, ])

idx <- with(site(f_us), !duplicated(paste(mukey, compname, comppct_r, compkind, localphase, drainagecl, erocl, slope_r)))
f_us <- f_us[which(idx), ]


# save(f_us, file = "gnatsgo.RData")
load(file = "gnatsgo.RData")



# tally Mollisols ----
co_us <- get_component_from_GDB(dsn = "D:/geodata/soils/gNATSGO_CONUS_FY20.gdb", WHERE = "compname LIKE '%'", childs = FALSE, stringsAsFactors = TRUE)
mu_us <- get_mapunit_from_GDB(dsn = "D:/geodata/soils/gNATSGO_CONUS_FY20.gdb", stats = FALSE)

test <- aggregate(comppct_r ~ mukey, data = co_us, function(x) sum(x, na.rm = TRUE))
idx  <- subset(test, comppct_r > 100)$mukey
co_us <- co_us[! co_us$mukey %in% idx, ]

subgroups <- c("Mollic|Udollic|Ustollic|Xerollic|Lithic Mollic|Pachic")

co_mol <- co_us %>%
  # filter(taxorder == "Mollisols") %>%
  group_by(mukey) %>%
  summarize(
    pct_mol  = sum(comppct_r[taxorder == "Mollisols"], na.rm = TRUE),
    pct_mol2 = sum(comppct_r[grepl(subgroups, taxsubgrp) & taxorder != "Mollisols"], na.rm = TRUE)
    ) %>%
  inner_join(mu_us, by = "mukey") %>%
  mutate(acres_mol  = muacres * pct_mol  / 100,
         acres_mol2 = muacres * pct_mol2 / 100
         )
sum(co_mol$acres_mol, na.rm = TRUE) / sum(mu_us$muacres, na.rm = TRUE)
sum(co_mol$acres_mol2, na.rm = TRUE) / sum(mu_us$muacres, na.rm = TRUE)
sum(mu_us$muacres[grepl("^Prime", mu_us$farmlndcl)], na.rm = TRUE) / sum(mu_us$muacres, na.rm = TRUE)


# OCONUS
asym <- c("AS", "FM", "GU", "HI", "MH", "MP", "PR", "PW", "VI")
co_oconus <- lapply(asym, function(x) {
  cat("getting ", x, "\n")
  get_component_from_SDA(WHERE = paste0("areasymbol LIKE '", x, "%'"), childs = FALSE, duplicates = TRUE)
})
co_oconus <- do.call("rbind", co_oconus)

mu_oconus <- lapply(asym, function(x) {
  cat("getting ", x, "\n")
  get_mapunit_from_SDA(WHERE = paste0("areasymbol LIKE '", x, "%'"))
})
mu_oconus <- do.call("rbind", mu_oconus)

co_mol <- co_oconus %>%
  filter(taxorder == "Mollisols") %>%
  group_by(mukey) %>%
  summarize(pct_mol = sum(comppct_r, na.rm = TRUE)) %>%
  right_join(mu_oconus, by = "mukey") %>%
  mutate(acres_mol = muacres * pct_mol / 100)

co_mol %>%
  mutate(asym = substr(areasymbol, 1, 2)) %>%
  group_by(asym) %>%
  summarize(pct_mol = round(sum(acres_mol, na.rm = TRUE) / sum(muacres, na.rm = TRUE) * 100))



# AK
ak_mukey <- read_sf("D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_ak.shp")$MUKEY

co_ak <- get_component_from_SDA(WHERE = paste0("mukey IN ('", paste0(ak_mukey, collapse = "', '"), "')"), childs = FALSE, duplicates = TRUE)
mu_us <- get_mapunit_from_SDA(WHERE = "areasymbol LIKE 'US'")
mu_ak <- subset(mu_ak, mukey %in% ak_mukey)


co_mol <- co_ak %>%
  filter(taxorder == "Mollisols") %>%
  group_by(mukey) %>%
  summarize(pct_mol = sum(comppct_r, na.rm = TRUE)) %>%
  right_join(mu_ak, by = "mukey") %>%
  mutate(acres_mol = muacres * pct_mol / 100)

co_mol %>%
  mutate(asym = substr(areasymbol, 1, 2)) %>%
  group_by(asym) %>%
  summarize(pct_mol = round(sum(acres_mol, na.rm = TRUE) / sum(muacres, na.rm = TRUE) * 100))



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
# saveRDS(osd, "osd.rds")
osd <- readRDS("osd.rds")


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

# write.csv(osd_bol, file = "osd_bol.csv", row.names = FALSE)
osd_bol <- read.csv(file = "osd_bol.csv")


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
               BS  = sumbases / cec7 * 100,
               OC  = om / 1.724
               # m_chroma = ifelse(grepl("olls$|^mollic|moll", taxsubgrp), 3, 5),
               # m_value  = ifelse(grepl("olls$|^mollic|moll", taxsubgrp), 3, 5),
               # d_value  = ifelse(grepl("olls$|^mollic|moll", taxsubgrp), 5, 6)
               )

h2 <- allocate(object = h, hztop = "hzdept", hzbot = "hzdepb", pedonid = "cokey", BS = "BS", OC = "OC", CEC = "cec7", m_chroma = "m_chroma", m_value = "m_value", d_value = "d_value", to = "FAO Black Soil")


s <- site(f_us)
s <- merge(s, h2, by = "cokey", all.x = TRUE)
s$compname <- tolower(s$compname)
s <- merge(s, osd_bol, by.x = "compname", by.y = "id", all.x = TRUE)


# compute dominant condition for black soils by mukey
s2 <- s %>%
  group_by(mukey) %>%
  summarize(pct_bs1 = sum(comppct_r[BS1 == TRUE], na.rm = TRUE) / 100,
            pct_bs2 = sum(comppct_r[BS2 == TRUE], na.rm = TRUE) / 100,
            pct_mo = sum(comppct_r[taxorder == "Mollisols" & m_chroma <= 3 & m_value <= 3 & d_value <= 5], na.rm = TRUE) / 100
            ) %>%
  ungroup() %>%
  as.data.frame()


# write.csv(s2, file = "ssurgo_black_soils_v2.csv", row.names = FALSE)
s2 <- read.csv(file = "ssurgo_black_soils_v2.csv", stringsAsFactors = FALSE)



# rasterize ----

library(raster)


names(s2)[1] <- "ID"

gnatsgo <- "D:/geodata/soils/gnatsgo_fy20_30m.tif"
r <- raster(gnatsgo)
r2 <- clusterR(r, reclassify, args = list(rcl = s2), progress = "text")
r2 <- reclassify(r, s2, progress = "text")
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

vars <- c("pct_bs1", "pct_bs2")
lapply(vars, function(x) {
  cat(x, as.character(Sys.time()), "\n")
  # beginCluster(type = "SOCK")
  deratify(r, att = x, 
           filename = paste0("D:/geodata/project_data/gsp-bs/gnatsgo_fy20_1km_", x, ".tif"),
           options = c("COMPRESS=DEFLATE"), 
           overwrite = TRUE, progress = "text" 
  )
  # endCluster()
})



# statsgo ----

mu_statsgo <- get_mapunit_from_SDA(WHERE = "areasymbol = 'US'")

mu_statsgo$idx <- rep(1:5, length.out = nrow(mu_statsgo))  

temp <- split(mu_statsgo, mu_statsgo$idx)
temp <- lapply(temp, function(x) fetchSDA(WHERE = paste0("areasymbol = 'US' AND mukey IN ('", paste(x$mukey, collapse = "', '"), "')"), duplicates = TRUE, childs = FALSE))

f_statsgo <- aqp::combine(temp)

# saveRDS(f_statsgo, file = "statsgo.rds")               
f_statsgo <- readRDS(file = "statsgo.rds")

osd_bol <- read.csv(file = "osd_bol.csv")


# segment

osd_bol$id <- tolower(osd_bol$id)

f_sub <- segment(f_statsgo, intervals = c(0, 25))
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


s <- site(f_statsgo)
s$compname <- tolower(s$compname)
s <- merge(s, h2, by = "cokey", all.x = TRUE)
s <- merge(s, osd_bol, by.x = "compname", by.y = "id", all.x = TRUE)



# compute dominant condition for black soils by mukey
s2 <- s %>%
  group_by(mukey) %>%
  summarize(
    pct_mo = sum(comppct_r[taxorder == "Mollisols"], na.rm = TRUE) / 100,
    pct_mo_col = sum(comppct_r[taxorder == "Mollisols" & m_chroma <= 3 & m_value <= 3 & d_value <= 5], na.rm = TRUE) / 100
    # pct_bs = sum(comppct_r[BS2 == TRUE], na.rm = TRUE) / 100
  ) %>%
  ungroup() %>%
  mutate(mukey = as.character(mukey)) %>%
  as.data.frame()


# write.csv(s2, file = "statsgo_black_soils_v2.csv", row.names = FALSE)
s2 <- read.csv(file = "statsgo_black_soils_v2.csv", stringsAsFactors = FALSE)


statsgo_sf <- read_sf("D:/geodata/project_data/gsp-bs/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us.shp")
statsgo_sf <- left_join(statsgo_sf, s2, by = c("MUKEY" = "mukey"))
write_sf(statsgo_sf, dsn = "D:/geodata/project_data/gsp-bs/gsmsoilmu_a_us_bs.shp")
