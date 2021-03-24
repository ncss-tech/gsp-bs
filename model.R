
# load data ----

library(sf)
library(raster)
# library(terra)

setwd("D:/geodata/project_data/gsp-bs/data")

bs <- readRDS(file = "tdata_bs.rds")
bs_nr <- readRDS(file = "tdata_bs_noRaCA.rds")
sa <- read_sf("D:/geodata/soils/gSSURGO_CONUS_FY20_july.gdb", layer = "SAPOLYGON")
sa <- st_transform(sa, crs = 4326)

idx <- sapply(st_intersects(bs, sa), any)
bs <- bs[idx, ]

idx <- sapply(st_intersects(bs_nr, sa), any)
bs_nr <- bs_nr[idx, ]

bs2   <- subset(bs, !is.na(BS1) | !is.na(BS2))
bs$BS <- ifelse(bs$BS2 == TRUE, "BS2", "None")
bs$BS <- ifelse(bs$BS1 == TRUE, "BS1", bs$BS)
write_sf(bs2, dsn = "ldm_bs3.shp")

bs2 <- subset(bs, BS2 == TRUE & !is.na(BS1))
mapview::mapview(bs2, zcol = "BS1")



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

# set.seed(563)
# idx <- caret::createDataPartition(data$BS2, p = 0.75)$Resample1
# data$train <- 1:nrow(data) %in% idx

# saveRDS(data, file = "tdata_bs_geo.rds")
data <- readRDS("tdata_bs_geo.rds")

sum(data$BS2)
idx <- sample(data$peiid[data$BS2 == FALSE], 1489 * 3)
data_sub <- subset(data, BS2 == TRUE | peiid %in% idx)



# Feature selection ----

vars <- c(isric_vars, other_vars, ssurgo_vars[c(1, 13, 16)])

data1 <- data[names(data) %in% c("BS1", vars)]
data1 <- na.exclude(data1)

bs1_fs_bor <- Boruta::Boruta(x = data1[, names(data1) %in% vars], y = data1[, "BS1"], maxRuns = 35, doTrace = 1)


data2 <- data[names(data) %in% c("BS2", vars, "train")]
data2 <- na.exclude(data2)
# idx   <- data2$train == TRUE

bs2_fs_bor <- Boruta::Boruta(x = data2[, names(data2) %in% vars], y = data2[, "BS2"], maxRuns = 35, doTrace = 1)


bs1_fs <- Boruta::attStats(bs1_fs_bor)
bs1_fs <- with(bs1_fs, data.frame(
  variable = row.names(bs1_fs), 
  varImp   = medianImp, 
  decision = decision, 
  stringsAsFactors = FALSE
  ))
bs1_fs$variable <- with(bs1_fs, reorder(as.factor(variable), varImp, sum))
bs1_fs <- bs1_fs[order(- bs1_fs$varImp), ]
row.names(bs1_fs) <- 1:nrow(bs1_fs)


bs2_fs <- Boruta::attStats(bs2_fs_bor)
bs2_fs <- with(bs2_fs, data.frame(
  variable = row.names(bs2_fs), 
  varImp   = medianImp, 
  decision = decision, 
  stringsAsFactors = FALSE
))
bs2_fs$variable <- with(bs2_fs, reorder(as.factor(variable), varImp, sum))
bs2_fs <- bs2_fs[order(- bs2_fs$varImp), ]
row.names(bs2_fs) <- 1:nrow(bs2_fs)

# save(bs1_fs_bor, bs2_fs_bor, bs1_fs, bs2_fs, file = "fselection_bs.RData")
load(file = "fselection_bs.RData")


cbind(as.character(bs1_fs$variable[1:25]), as.character(bs2_fs$variable[1:25]))


library(ggplot2)

ggplot(bs1_fs[1:25, ], aes(x = varImp, y = variable)) +
  geom_point()

ggplot(bs2_fs[1:25, ], aes(x = varImp, y = variable)) +
  geom_point()



# test down sampling balance model ----

library(ranger)
library(caret)


# modify BS class to switch
data$BS <- ifelse(data$BS2 == TRUE, "Yes", "No")
vars <- c("BS", as.character(bs2_fs$variable[c(1:25)]))

pct <- seq(0.1, 1, 0.1)
bs_train <- lapply(pct, function(x) {
  
  set.seed(534)
  
  # selet variables
  db   <- na.exclude(data[, names(data) %in% vars])
  
  # select samples
  n1   <- round(nrow(db) * 0.75)
  # idx1 <- sample(1:nrow(db), n1)
  idx1 <- caret::createDataPartition(db$BS, p = 0.75)$Resample1
  
  td  <- db[ idx1, ]
  vd  <- db[-idx1, ]
  # td <- db
  # vd <- db
  
  
  # balance samples a bit
  idx2 <- td$BS == "No"
  n2   <- round(nrow(td[idx2, ]) * x)
  rn  <- sample(row.names(td[idx2, ]), n2)
  td  <- td[td$BS == "Yes" | row.names(td) %in% rn, ]
  
  # ranger
  bs_train <- ranger(y = td$BS == "Yes", x = td[names(td) != "BS"],
                  probability = TRUE, # splitrule = "hellinger",
                  importance = "permutation"
                  )
  td$pred <- predict(bs_train, td)$predictions[, 2]
  # td$pred <- predict(bs_train, td)$predictions
  td_cm <- confusionMatrix(table(pred = td$pred > 0.5, obs = td$BS == "Yes"), positive = 'TRUE')
  
  pred    <- as.data.frame(predict(bs_train, vd)$predictions)
  names(pred) <- c("No", "Yes")
  pred$BS <- vd$BS
  vd$pred <- pred[, 2]
  # vd$pred <- predict(bs_train, vd)$predictions
  vd_cm <- confusionMatrix(table(pred = vd$pred > 0.5, obs = vd$BS == "Yes"), positive = 'TRUE')

  # # caret
  # bs_train <- train(y = td$BS, x = td[, names(td) != "BS"],
  #                   method = "ranger", # probability = TRUE,
  #                   importance = "permutation",
  #                   tuneGrid = expand.grid(
  #                     min.node.size = 10,
  #                     splitrule     = "gini",
  #                     mtry          = floor(sqrt(25))
  #                     ),
  #                   trControl = trainControl(
  #                     method = "cv", number = 10, returnResamp = "all", savePredictions = TRUE,
  #                     classProbs = TRUE,
  #                     summaryFunction = twoClassSummary,
  #                     # search = "random",
  #                     verboseIter = FALSE
  #                   )
  # )
  # td_cm <- confusionMatrix(table(pred = predict(bs_train$finalModel, td)$predictions[, 2] > 0.5, obs = td$BS == "Yes"), positive = 'TRUE')
  # vd_cm <- confusionMatrix(table(pred = predict(bs_train$finalModel, vd)$predictions[, 2] > 0.5, obs = vd$BS == "Yes"), positive = 'TRUE')
  # vd$pred <- predict(bs_train$finalModel, vd)$predictions[, 2]

  m1 <- mean(vd$pred[which(vd$BS == "Yes")], na.rm = TRUE)
  m0 <- mean(vd$pred[which(vd$BS == "No")], na.rm = TRUE)
  tjur_d <- abs(m1 - m0)
  R2 <- cor(vd$BS == "Yes", vd$pred)^2
  # R22 <- mean(test[vd$BS]) - mean(test[!vd$BS])
  
  dev  <- sigr::calcDeviance(vd$pred, vd$BS == "Yes")
  dev0 <- sigr::calcDeviance(mean(vd$BS == "Yes"), vd$BS == "Yes")
  D2   <- 1 - dev / dev0
  
  
  # BS2 = mean(vd$pred - as.integer(vd$BS == "Yes")^2)
  BS = aqp::brierScore(pred, c("No", "Yes"), "BS")
  # SE = aqp::shannonEntropy(pred$Yes)
  
  AUC <- vip::metric_auc(vd$BS == "Yes", vd$pred)
  # R22 <- mean(vd$probabilit[trueI,]) - mean(vd$probability[falseI,])
  
  return(list(data = td, rf = bs_train, td_cm = td_cm, vd_cm = vd_cm, R2 = R2, D2 = D2, tjur_d = tjur_d, AUC = AUC, BS = BS))
  })
# saveRDS(bs_train, file = "downsamp_ranger_bs1.rds")
# saveRDS(bs_train, file = "downsamp_caret_bs1.rds")
# saveRDS(bs_train, file = "downsamp_ranger_bs2.rds")
# saveRDS(bs_train, file = "downsamp_caret_bs2.rds")
# bs_train <- readRDS(file = "downsamp_ranger_bs1.rds")
# saveRDS(bs_train, file = "downsamp_caret_bs1.rds")
# saveRDS(bs_train, file = "downsamp_ranger_bs2.rds")
# saveRDS(bs_train, file = "downsamp_caret_bs2.rds")


cm <- sapply(bs_train, function(x) round(x$vd_cm$byClass, 2))
cm <- as.data.frame(t(cm))
cm$pct <- seq(0.1, 1, 0.1)
cm


library(ggplot2)

ggplot(cm, aes(x = pct)) +
  geom_line(aes(y = `Precision`), col = "red") + 
  geom_line(aes(y = `Recall`), col = "blue") + 
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(0, 1, 0.1))

ov <- sapply(bs_train, function(x) round(x$vd_cm$overall, 2))
ov <- as.data.frame(t(ov))
ov$pct <- pct

ggplot(ov, aes(x = pct)) +
  geom_line(aes(y = Accuracy), col = "red") + 
  geom_line(aes(y = Kappa), col = "blue") + 
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(0, 1, 0.1))


vars <- as.character(s$variable[c(1:25)])
vd   <- na.exclude(data[data$train == FALSE, c("BS", vars)])
vd$pred <- predict(bs_train[[3]]$rf, vd, type = "prob")[, 2]
pred <- prediction(
  prediction = vd$pred, 
  labels     = vd$BS
)
plot(performance(pred, "t"))
abline(0, 1)



# caret ----

data$BS <- data$BS1

# balance samples a bit
idx <- data$BS == FALSE
n   <- round(nrow(data[idx, ]) * 0.3)
idx <- sample(data$peiid[idx], n)
db  <- subset(data, BS == TRUE | peiid %in% idx)

vars <- c("BS", as.character(bs1_fs$variable[c(1:25)]))
# td <- na.exclude(data[data$train == TRUE, names(data) %in% vars])
td <- na.exclude(db[,   names(db) %in% vars])
vd <- na.exclude(data[, names(data) %in% vars])

td$BS <- ifelse(td$BS == TRUE, "Yes", "No")
bs_train <- train(y = td$BS, x = td[, names(td) != "BS"],
                  method = "ranger", # probability = TRUE,
                  importance = "permutation",
                  tuneGrid = expand.grid(
                    min.node.size = 10,
                    splitrule     = "gini",
                    mtry          = floor(sqrt(25))
                  ),
                  trControl = trainControl(
                    method = "cv", number = 10, returnResamp = "all", savePredictions = TRUE,
                    # classProbs = TRUE,
                    summaryFunction = twoClassSummary,
                    # search = "random",
                    verboseIter = FALSE
                    )
                  )


# save(bs_train, file = "bs_train30_ranger.rds")
# save(bs_train, file = "bs_train30_caret1010.rds")

confusionMatrix(table(pred = predict(bs_train, td), obs = td$BS), positive = "Yes")
vd_cm1 <- confusionMatrix(table(obs = ifelse(vd$BS == TRUE, "Yes", "No"), pred = predict(bs_train, vd)), positive = "Yes")
vd_cm2 <- confusionMatrix(table(obs = vd$BS, pred = predict(bs_train, vd, type = "prob")$Yes > 0.5), positive = "TRUE")



# Generate map ----

predfun <- function(model, ...) predict(model, ...)$predictions[, 2]
predfun <- function(model, ...) predict(model, ...)$predictions

# rs2 <- raster::stack(rs)
names(rs) <- gsub("_CONUS", "", names(rs))
idx <- which(
  names(rs) %in% c("BS2", as.character(s$variable[c(1:25)]))
  # names(rs2) %in% bs_train$finalModel$forest$independent.variable.names
)
rs <- rs[[idx]]

bs_r <- predict(rs, bs_train[[6]]$rf, fun = predfun, index = 1, progress = "text", overwrite = TRUE, filename = "D:/geodata/project_data/gsp-bs/bs_dsm_60pct.tif")
# bs_r <- predict(rs, bs_train$finalModel, fun = predfun, index = 1, progress = "text", overwrite = TRUE, filename = "D:/geodata/project_data/gsp-bs/test_raca.tif")

image(bs_r > 0.5)
tb_dsm <- table(values(bs_r) > 0.5)
prop.table(tb_dsm)



# Calculate uncertainty maps ----

# modify BS class to switch, and manually set the down sampling threshold (0.20 for BS1, 0.40 for BS2)
data$BS <- data$BS2
vars <- c("BS", as.character(bs2_fs$variable[c(1:25)]))
pct     <- 0.6

i       <- 1:50
bs_unc <- lapply(i, function(x) {
  
  set.seed(x)
  
  # select variables
  db   <- na.exclude(data[, names(data) %in% vars])
  
  
  # select samples
  n   <- round(nrow(db) * 0.75)
  # idx <- sample(1:nrow(db), n)
  idx <- caret::createDataPartition(db$BS, p = 0.75)$Resample1
  td  <- db[ idx, ]
  vd  <- db[-idx, ]
  
  
  # balance samples a bit
  idx <- td$BS == FALSE
  n   <- round(nrow(td[idx, ]) * pct)
  rn  <- sample(row.names(td[idx, ]), n)
  td  <- td[td$BS == TRUE | row.names(td) %in% rn, ]
  
  
  # ranger
  bs_train <- ranger(y = td$BS, x = td[names(td) != "BS"],
                     probability = TRUE,
                     importance = "permutation"
  )
  
  td$pred <- predict(bs_train, td)$predictions[, 2]
  td_cm <- confusionMatrix(table(pred = td$pred > 0.5, obs = td$BS), positive = 'TRUE')
  
  pred    <- as.data.frame(predict(bs_train, vd)$predictions)
  names(pred) <- c("No", "Yes")
  pred$BS <- vd$BS
  vd$pred <- predict(bs_train, vd)$predictions[, 2]
  vd_cm <- confusionMatrix(table(pred = vd$pred > 0.5, obs = vd$BS), positive = 'TRUE')
  
  m1 <- mean(vd$pred[which(vd$BS == TRUE)], na.rm = TRUE)
  m0 <- mean(vd$pred[which(vd$BS == FALSE)], na.rm = TRUE)
  tjur_d <- abs(m1 - m0)
  R2  <- cor(vd$BS, vd$pred)^2
  
  dev  <- sigr::calcDeviance(vd$pred, vd$BS)
  dev0 <- sigr::calcDeviance(mean(vd$BS), vd$BS)
  D2   <- 1 - dev / dev0
  
  BS = mean(vd$pred - as.integer(vd$BS == "Yes")^2)
  # BS = aqp::brierScore(pred, c("No", "Yes"), "BS")
  
  # R22 <- mean(test[vd$BS]) - mean(test[!vd$BS])
  AUC <- vip::metric_auc(vd$BS, vd$pred)
  
  
  return(list(
    td     = td,
    vd     = vd,
    rf     = bs_train, 
    td_cm  = td_cm, 
    vd_cm  = vd_cm, 
    R2     = R2, 
    D2     = D2,
    tjur_d = tjur_d,
    AUC    = AUC,
    BS     = BS
    ))
})
# saveRDS(bs_unc, file = "bs1_unc_60.rds")
# saveRDS(bs_unc, file = "bs2_unc_60.rds")
# bs_unc <- readRDS(file = "bs1_unc_60.rds")
# bs_unc <- readRDS(file = "bs2_unc_60.rds")


summary(t(sapply(bs_unc, function(x) x$vd_cm$byClass)))
summary(t(sapply(bs_unc, function(x) x$vd_cm$overall)))
summary(sapply(bs_unc, function(x) x$AUC))
summary(sapply(bs_unc, function(x) x$D2))
summary(sapply(bs_unc, function(x) x$BS))
summary(sapply(bs_unc, function(x) round(x$tjur_d, 2)))



# change BS
bs1_unc_60 <- readRDS(file = "bs1_unc_60.rds")
bs2_unc_60 <- readRDS(file = "bs2_unc_60.rds")

bs_unc <- bs1_unc_60
vars <- as.character(bs1_fs$variable[c(1:25)])

rs2 <- rs
names(rs2) <- gsub("_CONUS", "", names(rs2))
idx <- which(
  names(rs2) %in% vars
  # names(rs2) %in% bs_train$finalModel$forest$independent.variable.names
)
rs2 <- readAll(rs2[[idx]])

i <- 1:50
lapply(i, function(x) {
  predfun <- function(model, ...) predict(model, ...)$predictions[, 2]
  bs_r <- predict(rs2, bs_unc[[x]]$rf, fun = predfun, index = 1, progress = "text", overwrite = TRUE, filename = paste0("D:/geodata/project_data/gsp-bs/bs1_dsm_", x, ".tif"))
})



# BS 1
fp <- "D:/geodata/project_data/gsp-bs"
lf <- list.files(fp)
lf <- lf[grepl("bs1_dsm_[1-9]", lf)]
bs_rs <- readAll(stack(file.path(fp, lf)))

bs_mean <- calc(bs_rs, mean, filename = file.path(fp, "bs1_dsm_mean.tif"), progress = "text", overwrite = TRUE)
bs_sd   <- calc(bs_rs, sd, filename = file.path(fp, "bs1_dsm_sd.tif"), progress = "text", overwrite = TRUE)
bs_pi   <- calc(bs_rs, function(x) {
  bs_q <- quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  bs_q[[2]] - bs_q[[1]]
  }, 
  filename = file.path(fp, "bs1_dsm_pi.tif"), progress = "text", overwrite = TRUE
  )
writeRaster(bs_rs, file = file.path(fp, "bs1_dsm_iterations_ds60.tif"), progress = "text", overwrite = TRUE)
file.remove(file.path(fp, lf))





# BS 2
fp <- "D:/geodata/project_data/gsp-bs"
lf <- list.files(fp)
lf <- lf[grepl("bs2_dsm_[1-9]", lf)]
bs_rs <- readAll(stack(file.path(fp, lf)))


bs_mean <- calc(bs_rs, mean, filename = file.path(fp, "bs2_dsm_mean.tif"), progress = "text", overwrite = TRUE)
bs_sd   <- calc(bs_rs, sd, filename = file.path(fp, "bs2_dsm_sd.tif"), progress = "text", overwrite = TRUE)
bs_pi   <- calc(bs_rs, function(x) {
  bs_q <- quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  bs_q[[2]] - bs_q[[1]]
}, 
filename = file.path(fp, "bs2_dsm_pi.tif"), progress = "text", overwrite = TRUE
)
writeRaster(bs_rs, file = file.path(fp, "bs2_dsm_iteration_ds60.tif"), progress = "text", overwrite = TRUE)
file.remove(file.path(fp, lf))



# mask
lc <- raster("D:/geodata/project_data/gsp-sas/1km covariates/ISRIC/CONUS_redundant/LCEE10__CONUS.tif")
m <- c(
  0,   209, 1, 
  209, 210, NA,
  210, 300, 1
)
m <- matrix(m, ncol = 3, byrow = TRUE)
lc <- reclassify(lc, m, filename = "D:/geodata/project_data/gsp-bs/mask.tif", overwrite = TRUE, progress = "text")
lc <- raster("D:/geodata/project_data/gsp-bs/mask.tif")

bs1 <- raster(file.path(fp, "bs1_dsm_mean.tif"))
bs2 <- raster(file.path(fp, "bs2_dsm_mean.tif"))

bs1 <- mask(bs1, lc, filename = file.path(fp, "bs1_dsm_mean_mask.tif"), overwrite = TRUE, progress = "text")
bs2 <- mask(bs2, lc, filename = file.path(fp, "bs2_dsm_mean_mask.tif"), overwrite = TRUE, progress = "text")


# classify ----
m <- c(
  0,   0.5, 0, 
  0.5, 1,   1
)
m <- matrix(m, nrow = 2, byrow = TRUE)
bs1_class <- reclassify(bs1, m, filename = "D:/geodata/project_data/gsp-bs/bs1_dsm_class_mask.tif", overwrite = TRUE, progress = "text")

m <- c(
  0,   0.5, 0, 
  0.5, 1,   2
)
m <- matrix(m, nrow = 2, byrow = TRUE)
bs2_class <- reclassify(bs2, m, filename = "D:/geodata/project_data/gsp-bs/bs2_dsm_class_mask.tif", overwrite = TRUE, progress = "text")





# compute final confusion matrix ----
vars <- c(dsm_bs1 = "bs1_dsm_mean.tif",
          dsm_bs2 = "bs2_dsm_mean.tif",
          ssurgo_bs1 = "gnatsgo_fy20_1km_pct_bs1.tif",
          ssurgo_bs2 = "gnatsgo_fy20_1km_pct_bs2.tif"
          )
fp <- "D:/geodata/project_data/gsp-bs"
rs_bs <- stack(file.path(fp, vars))
names(rs_bs) <- c("dsm_bs1", "dsm_bs2", "ssurgo_bs1", "ssurgo_bs2")

pred <- extract(rs_bs, st_coordinates(bs))
pred <- cbind(bs[c("BS1", "BS2")], pred)

cm_dsm_bs1    <- confusionMatrix(table(pred = pred$dsm_bs1    > 0.5, obs = pred$BS2), positive = "TRUE")
cm_dsm_bs2    <- confusionMatrix(table(pred = pred$dsm_bs2    > 0.5, obs = pred$BS2), positive = "TRUE")
cm_ssurgo_bs1 <- confusionMatrix(table(pred = pred$ssurgo_bs1 > 0.5, obs = pred$BS2), positive = "TRUE")
cm_ssurgo_bs2 <- confusionMatrix(table(pred = pred$ssurgo_bs2 > 0.5, obs = pred$BS2), positive = "TRUE")

D2 <- function(pred, obs) {
  dev  = sigr::calcDeviance(pred, obs, na.rm = TRUE)
  dev0 = sigr::calcDeviance(mean(obs, na.rm = TRUE), obs, na.rm = TRUE)
  D2   = 1 - dev / dev0
  return(D2)
}
BS <- function(pred, obs) mean(pred - as.integer(obs)^2, na.rm = TRUE)

D2 <- rbind(
  D2(pred$dsm_bs1,    pred$BS1),
  D2(pred$dsm_bs2,    pred$BS2),
  D2(pred$ssurgo_bs1, pred$BS1),
  D2(pred$ssurgo_bs2, pred$BS2)
)
brier <- rbind(
  BS(pred$dsm_bs1,    pred$BS1),
  BS(pred$dsm_bs2,    pred$BS2),
  BS(pred$ssurgo_bs1, pred$BS1),
  BS(pred$ssurgo_bs2, pred$BS2)
)

data.frame(source   = c("DSM", "DSM", "SSURGO", "SSURGO"),
           category = c("BS1", "BS2", "BS1", "BS2"),
           D2, brier
           )

test <- data.frame(
  No  = 1 - pred$dsm_bs1, 
  Yes = pred$dsm_bs1, 
  BS  = ifelse(pred$BS1 > 0.5, "Yes", "No"), 
  stringsAsFactors = FALSE
  )
aqp::brierScore(na.exclude(test),  c("No", "Yes"), "BS")
