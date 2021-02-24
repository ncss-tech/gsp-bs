
# load data ----

library(sf)
library(terra)


bs <- readRDS(file = "C:/Users/stephen.roecker/OneDrive - USDA/data/bs_training_data_post1997.rds")
sa <- read_sf("D:/geodata/soils/gSSURGO_CONUS_FY20_july.gdb", layer = "SAPOLYGON")
sa <- st_transform(sa, crs = 4326)

idx <- sapply(st_intersects(bs, sa), any)
bs <- bs[idx, ]


# ISRIC CONUS
path <- "D:/geodata/project_data/gsp-sas/1km covariates/ISRIC/CONUS"
lf <- list.files(path)
lf <- lf[grepl(".tif$", lf)]
isric_rs <- rast(file.path(path, lf))
names(isric_rs) <- gsub("_CONUS.tif|.tif", "", lf)
isric_vars <- names(isric_rs)


# Other CONUS
path <- "D:/geodata/project_data/gsp-sas/1km covariates/Other"
lf <- list.files(path)
lf <- lf[grepl(".tif$", lf)]
# find matching extents
# other_extent <- sapply(file.path(path, lf), function(x) paste(extent(raster(x))))
# idx <- names(sort(table(other_extent), decreasing = TRUE)[1])
# idx <- which(other_extent == idx)
# lf  <- lf[idx]
other_rs <- rast(file.path(path, lf))
names(other_rs) <- gsub(".tif", "", lf)
other_vars <- names(other_rs)


# SSURGO CONUS
path <- "D:/geodata/project_data/gsp-sas/1km covariates/SSURGO"
lf <- list.files(path)
lf <- lf[grepl(".tif$", lf)]
# # find matching extents
# other_extent <- sapply(file.path(path, lf), function(x) paste(extent(raster(x))))
# idx <- names(sort(table(other_extent), decreasing = TRUE)[1])
# idx <- which(other_extent == idx)
# lf  <- lf[idx]
ssurgo_rs <- rast(file.path(path, lf))
names(ssurgo_rs) <- gsub(".tif", "", lf)
ssurgo_vars <- names(ssurgo_rs)

rs <- c(isric_rs, other_rs, ssurgo_rs)



# Extract geodata ----

data <- extract(rs, st_coordinates(bs))
data <- cbind(bs, as.data.frame(data))
data$geometry <- NULL
data <- within(data, {
  ssurgo_pm_lu3_rs_mode_v2 = as.factor(ssurgo_pm_lu3_rs_mode_v2)
  statsgo_mlra             = as.factor(statsgo_mlra)
  statsgo_pmkind3          = as.factor(statsgo_pmkind3)
  mlra_rs                  = as.factor(mlra_rs)
  nlcd_2016_conus          = as.factor(nlcd_2016_conus)
})
data$decade <- floor(data$year * 0.1) * 10


# Split data ----

set.seed(563)
idx <- caret::createDataPartition(data$BS2, p = 0.75)$Resample1
data$train <- 1:nrow(data) %in% idx



# Feature selection ----

library(Boruta)

vars <- c(isric_vars, other_vars, ssurgo_vars[c(1, 13, 16)])

data2 <- data[names(data) %in% c("BS2", vars, "train")]
data2 <- na.exclude(data2)

idx <- data2$train == TRUE
data_fs <- Boruta(x = data2[idx, names(data2) %in% vars], y = data2[idx, "BS2"], maxRuns = 35, doTrace = 1)

# saveRDS(data_fs, file = "C:/Users/stephen.roecker/OneDrive - USDA/data/bs_training_data_fs.rds")
data_fs <- readRDS(file = "C:/Users/stephen.roecker/OneDrive - USDA/data/bs_training_data_fs.rds")


s <- attStats(data_fs)
s <- data.frame(variable = row.names(s), varImp = s$medianImp, decision = s$decision, stringsAsFactors = FALSE)
s$variable <- reorder(as.factor(s$variable), s$varImp, sum)
s <- s[order(- s$varImp), ]
row.names(s) <- 1:nrow(s)


library(ggplot2)

ggplot(s[1:25, ], aes(x = varImp, y = variable)) +
  geom_point()




# Train model ----

library(ranger)
library(caret)

vars <- c("BS2", as.character(s$variable[c(1:25)]))
td <- na.exclude(data[data$train == TRUE, names(data) %in% vars])

vd <- na.exclude(data[data$train == FALSE, names(data) %in% vars])


# ranger
bs_train <- ranger(y = td$BS, x = td[-1],
                   probability = TRUE,
                   importance = "permutation",
)

confusionMatrix(table(td$BS2, predict(bs_train, td)$predictions[, 2] > 0.4))
confusionMatrix(table(vd$BS2, predict(bs_train, vd)$predictions[, 2] > 0.5))


# caret

td$BS2 <- ifelse(td$BS2 == TRUE, "Yes", "No")
vd$BS2 <- ifelse(vd$BS2 == TRUE, "Yes", "No")

bs_train <- train(y = td$BS, x = td[-1],
                  method = "ranger", classification = TRUE,
                  importance = "permutation",
                  trControl = trainControl(
                    method = "cv", number = 10, returnResamp = "all", savePredictions = TRUE,
                    search = "random", verboseIter = FALSE
                  )
                  )

confusionMatrix(table(td$BS2, predict(bs_train, td)), positive = "Yes")
confusionMatrix(table(vd$BS2, predict(bs_train, vd)), positive = "Yes")



# Generate map

predfun <- function(model, ...) predict(model, ...)$predictions[, 2]
predfun <- function(model, ...) predict(model, ...)$predictions

rs2 <- raster::stack(rs)
names(rs2) <- gsub("_CONUS", "", names(rs2))
idx <- which(
  names(rs2) %in% bs_train$forest$independent.variable.names
  # names(rs2) %in% bs_train$finalModel$forest$independent.variable.names
)
rs2 <- rs2[[idx]]

bs_r <- predict(rs2, bs_train, fun = predfun, index = 1, progress = "text", overwrite = TRUE, filename = "C:/workspace2/test_raca.tif")
# bs_r <- predict(rs, bs_train$finalModel, fun = predfun, index = 1, progress = "text", overwrite = TRUE, filename = "C:/workspace2/test_raca.tif")


