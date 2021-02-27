
# load data ----

library(sf)
library(raster)
# library(terra)

setwd("C:/Users/stephen.roecker/OneDrive - USDA/data")

bs <- readRDS(file = "bs_training_data_post1997.rds")
sa <- read_sf("D:/geodata/soils/gSSURGO_CONUS_FY20_july.gdb", layer = "SAPOLYGON")
sa <- st_transform(sa, crs = 4326)

idx <- sapply(st_intersects(bs, sa), any)
bs <- bs[idx, ]


# ISRIC CONUS
path <- "D:/geodata/project_data/gsp-sas/1km covariates/ISRIC/CONUS"
lf <- list.files(path)
lf <- lf[grepl(".tif$", lf)]
isric_rs <- stack(file.path(path, lf))
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
other_rs <- stack(file.path(path, lf))
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
ssurgo_rs <- stack(file.path(path, lf))
names(ssurgo_rs) <- gsub(".tif", "", lf)
ssurgo_vars <- names(ssurgo_rs)

rs <- stack(isric_rs, other_rs, ssurgo_rs)



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

# saveRDS(data, file = "data.rds")
data <- readRDS("data.rds")

sum(data$BS2)
idx <- sample(data$peiid[data$BS2 == FALSE], 1489 * 3)
data_sub <- subset(data, BS2 == TRUE | peiid %in% idx)


# Feature selection ----

vars <- c(isric_vars, other_vars, ssurgo_vars[c(1, 13, 16)])

data2 <- data[names(data) %in% c("BS2", vars, "train")]
data2 <- na.exclude(data2)

idx <- data2$train == TRUE
data_fs <- Boruta::Boruta(x = data2[idx, names(data2) %in% vars], y = data2[idx, "BS2"], maxRuns = 35, doTrace = 1)

# saveRDS(data_fs, file = "C:/Users/stephen.roecker/OneDrive - USDA/data/bs_training_data_fs.rds")
data_fs <- readRDS(file = "C:/Users/stephen.roecker/OneDrive - USDA/data/bs_training_data_fs.rds")

s <- Boruta::attStats(data_fs)
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


pct <- seq(0.1, 1, 0.1)

bs_train <- lapply(pct, function(x) {
  
  # balance samples a bit
  n   <- round(nrow(data[data$BS == FALSE]) * 0.5)
  idx <- sample(data$peiid[idx], n)
  db  <- subset(data, BS2 == TRUE | peiid %in% idx)
  
  vars <- c("BS2", as.character(s$variable[c(1:25)]))
  # td <- na.exclude(data[data$train == TRUE, names(data) %in% vars])
  td <- na.exclude(db[db$train == TRUE,  names(db) %in% vars])
  vd <- na.exclude(data[data$train == FALSE, names(data) %in% vars])
  
  # ranger
  bs_train <- ranger(y = td$BS2, x = td[-1],
                   probability = TRUE,
                   importance = "permutation"
                   )
  
  td_cm <- confusionMatrix(table(obs = td$BS2, pred = predict(bs_train, td)$predictions[, 2] > 0.5), positive = 'TRUE')
  vd_cm <- confusionMatrix(table(obs = vd$BS2, pred = predict(bs_train, vd)$predictions[, 2] > 0.5), positive = 'TRUE')
  
  R2 <- cor(vd$BS2, predict(bs_train, vd)$predictions[, 2])^2
  test <- predict(bs_train, vd)$predictions[, 2]
  R22 <- mean(test[vd$BS2]) - mean(test[!vd$BS2])
  AUC <- vip::metric_auc(vd$BS2, test)
  # R22 <- mean(vd$probabilit[trueI,]) - mean(vd$probability[falseI,])
  
  return(list(data = data_sub, rf = bs_train, td_cm = td_cm, vd_cm = vd_cm, R2 = R2, R22 = R22, AUC = AUC))
  })

cm <- sapply(bs_train, function(x) round(x$vd_cm$byClass, 2))
cm <- as.data.frame(t(cm))
cm$pct <- pct

ggplot(cm, aes(x = pct)) +
  geom_line(aes(y = `Sensitivity`), col = "red") + 
  geom_line(aes(y = `Precision`), col = "blue") + 
  ylim(0, 1)

ov <- sapply(bs_train, function(x) round(x$vd_cm$overall, 2))
ov <- as.data.frame(t(ov))
ov$pct <- pct

ggplot(ov, aes(x = pct)) +
  geom_line(aes(y = Accuracy), col = "red") + 
  geom_line(aes(y = Kappa), col = "blue") + 
  
  ylim(0, 1)


# caret

vars <- c("BS2", as.character(s$variable[c(1:25)]))
data_sub     <- na.exclude(data[names(data) %in% vars])
data_sub$BS2 <- ifelse(data_sub$BS2 == TRUE, "Yes", "No")

bs_train <- train(y = data_sub$BS, x = data_sub[-1],
                  method = "ranger",
                  importance = "permutation",
                  trControl = trainControl(
                    method = "cv", number = 10, returnResamp = "all", savePredictions = TRUE,
                    search = "random", verboseIter = FALSE
                  )
                  )

confusionMatrix(table(obs = td$BS2, pred = predict(bs_train, td)), positive = "Yes")
confusionMatrix(table(obs = vd$BS2, pred = predict(bs_train, vd)), positive = "Yes")



# Calculate uncertainty maps ----

i <- 1:50
bs_unc <- lapply(i, function(x) {
  
  set.seed(x)
  
  # balance samples a bit
  n   <- round(nrow(data[data$BS == FALSE]) * 0.5)
  idx <- sample(data[data$BS2 == FALSE, "peiid"], n)
  db  <- subset(data, BS2 == TRUE | peiid %in% idx)
  
  
  # select samples
  n   <- round(nrow(db) * 0.75)
  idx <- sample(db$peiid, n)
  td  <- subset(db,  peiid %in% idx)
  vd  <- subset(db, !peiid %in% idx)
  
  
  # select variables
  vars <- c("BS2", as.character(s$variable[c(1:25)]))
  td <- na.exclude(td[, names(td) %in% vars])
  vd <- na.exclude(vd[, names(vd) %in% vars])
  
  
  # ranger
  bs_train <- ranger(y = td$BS2, x = td[-1],
                     probability = TRUE,
                     importance = "permutation"
  )
  
  td_cm <- confusionMatrix(table(obs = td$BS2, pred = predict(bs_train, td)$predictions[, 2] > 0.5), positive = 'TRUE')
  vd_cm <- confusionMatrix(table(obs = vd$BS2, pred = predict(bs_train, vd)$predictions[, 2] > 0.5), positive = 'TRUE')
  
  vd$pred <- predict(bs_train, vd)$predictions[, 2]
  R2  <- cor(vd$BS2, vd$pred)^2
  # R22 <- mean(test[vd$BS2]) - mean(test[!vd$BS2])
  AUC <- vip::metric_auc(vd$BS2, vd$pred)
  
  
  return(list(
    td    = td,
    vd    = vd,
    rf    = bs_train, 
    td_cm = td_cm, 
    vd_cm = vd_cm, 
    R2    = R2, 
    AUC   = AUC
    ))
})
# saveRDS(bs_unc, file = "bs_unc.rds")
bs_unc <- readRDS(file = "bs_unc.rds")


rs2 <- rs
names(rs2) <- gsub("_CONUS", "", names(rs2))
idx <- which(
  names(rs2) %in% c("BS2", as.character(s$variable[c(1:25)]))
  # names(rs2) %in% bs_train$finalModel$forest$independent.variable.names
)
rs2 <- rs2[[idx]]

i <- 26:50
lapply(i, function(x) {
  predfun <- function(model, ...) predict(model, ...)$predictions[, 2]
  bs_r <- predict(rs2, bs_unc[[x]]$rf, fun = predfun, index = 1, progress = "text", overwrite = TRUE, filename = paste0("D:/geodata/project_data/gsp-bs/bs_dsm_", x, ".tif"))
})



# Generate map ----

predfun <- function(model, ...) predict(model, ...)$predictions[, 2]
predfun <- function(model, ...) predict(model, ...)$predictions

rs2 <- raster::stack(rs)
names(rs2) <- gsub("_CONUS", "", names(rs2))
idx <- which(
  names(rs2) %in% c("BS2", as.character(s$variable[c(1:25)]))
  # names(rs2) %in% bs_train$finalModel$forest$independent.variable.names
)
rs2 <- rs2[[idx]]

bs_r <- predict(rs2, bs_train[[10]]$rf, fun = predfun, index = 1, progress = "text", overwrite = TRUE, filename = "D:/geodata/project_data/gsp-bs/test_raca_1.tif")
# bs_r <- predict(rs, bs_train$finalModel, fun = predfun, index = 1, progress = "text", overwrite = TRUE, filename = "C:/workspace2/test_raca.tif")

image(bs_r > 0.5)
tb_dsm <- table(values(bs_r) > 0.5)
prop.table(tb_dsm)


# Prediction interval