
# load data ----

library(sf)
library(raster)
# library(terra)

setwd("C:/Users/stephen.roecker/OneDrive - USDA/data")

bs <- readRDS(file = "bs_training_data_all.rds")
sa <- read_sf("D:/geodata/soils/gSSURGO_CONUS_FY20_july.gdb", layer = "SAPOLYGON")
sa <- st_transform(sa, crs = 4326)

idx <- sapply(st_intersects(bs, sa), any)
bs <- bs[idx, ]

mapview::mapview(bs, zcol = "BS2")



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

# idx <- data2$train == TRUE
# data_fs <- Boruta::Boruta(x = data2[idx, names(data2) %in% vars], y = data2[idx, "BS2"], maxRuns = 35, doTrace = 1)
data_fs <- Boruta::Boruta(x = data2[names(data2) %in% vars], y = data2[, "BS2"], maxRuns = 35, doTrace = 1)

# saveRDS(data_fs, file = "C:/Users/stephen.roecker/OneDrive - USDA/data/bs_training_data_fs.rds")
# saveRDS(data_fs, file = "C:/Users/stephen.roecker/OneDrive - USDA/data/bs_training_data_fs_v2.rds")
data_fs <- readRDS(file = "C:/Users/stephen.roecker/OneDrive - USDA/data/bs_training_data_fs.rds")
data_fs <- readRDS(file = "C:/Users/stephen.roecker/OneDrive - USDA/data/bs_training_data_fs_v2.rds")

s <- Boruta::attStats(data_fs)
s <- data.frame(variable = row.names(s), varImp = s$medianImp, decision = s$decision, stringsAsFactors = FALSE)
s$variable <- reorder(as.factor(s$variable), s$varImp, sum)
s <- s[order(- s$varImp), ]
row.names(s) <- 1:nrow(s)

s2 <- s

library(ggplot2)

ggplot(s[1:25, ], aes(x = varImp, y = variable)) +
  geom_point()




# test down sampling balance model ----

library(ranger)
library(caret)


pct <- seq(0.1, 1, 0.1)

bs_train <- lapply(pct, function(x) {
  
  # balance samples a bit
  idx <- data$BS2 == FALSE
  n   <- round(nrow(data[idx, ]) * x)
  idx <- sample(data$peiid[idx], n)
  db  <- subset(data, BS2 == TRUE | peiid %in% idx)
  
  vars <- c("BS2", as.character(s$variable[c(1:25)]))
  # td <- na.exclude(data[data$train == TRUE, names(data) %in% vars])
  # n   <- nrow(db)
  # idx <- sample(1:n, n * 0.75)
  td  <- na.exclude(db[  data$train == TRUE,  names(db) %in% vars])
  vd  <- na.exclude(data[data$train == FALSE, names(data) %in% vars])
  
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
  
  return(list(data = td, rf = bs_train, td_cm = td_cm, vd_cm = vd_cm, R2 = R2, R22 = R22, AUC = AUC))
  })

cm <- sapply(bs_train, function(x) round(x$vd_cm$byClass, 2))
cm <- as.data.frame(t(cm))
cm$pct <- pct

ggplot(cm, aes(x = pct)) +
  geom_line(aes(y = `Sensitivity`), col = "red") + 
  geom_line(aes(y = `Precision`), col = "blue") + 
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



# caret

vars <- c("BS2", as.character(s$variable[c(1:25)]))
data_sub     <- na.exclude(data[names(data) %in% vars])

# balance samples a bit
idx <- data$BS2 == FALSE
n   <- round(nrow(data[idx, ]) * 0.5)
idx <- sample(data$peiid[idx], n)
db  <- subset(data, BS2 == TRUE | peiid %in% idx)

vars <- c("BS2", as.character(s$variable[c(1:25)]))
# td <- na.exclude(data[data$train == TRUE, names(data) %in% vars])
td <- na.exclude(db[,   names(db) %in% vars])
vd <- na.exclude(data[, names(data) %in% vars])
td$BS2 <- ifelse(td$BS2 == TRUE, "Yes", "No")

sf <- function(data, lev, model) {
  
  
}

bs_train <- train(y = td$BS2, x = td[-1],
                  method = "ranger",
                  importance = "permutation",
                  trControl = trainControl(
                    method = "cv", number = 10, returnResamp = "all", savePredictions = TRUE,
                    classProbs = TRUE,
                    search = "random", verboseIter = FALSE
                  )
                  )

confusionMatrix(table(obs = td$BS2, pred = predict(bs_train, td)), positive = "Yes")
confusionMatrix(table(obs = vd$BS2, pred = predict(bs_train, vd)), positive = "Yes")



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
# bs_r <- predict(rs, bs_train$finalModel, fun = predfun, index = 1, progress = "text", overwrite = TRUE, filename = "C:/workspace2/test_raca.tif")

image(bs_r > 0.5)
tb_dsm <- table(values(bs_r) > 0.5)
prop.table(tb_dsm)



# Calculate uncertainty maps ----

i <- 1:50
bs_unc <- lapply(i, function(x) {
  
  set.seed(x)
  
  # selet variables
  vars <- c("BS2", as.character(s$variable[c(1:25)]))
  db   <- na.exclude(data[, names(data) %in% vars])
  
  
  # select samples
  n   <- round(nrow(db) * 0.75)
  idx <- sample(1:nrow(db), n)
  td  <- db[ idx, ]
  vd  <- db[-idx, ]
  
  
  # balance samples a bit
  idx <- td$BS2 == FALSE
  n   <- round(nrow(td[idx, ]) * 0.6)
  rn  <- sample(row.names(td[idx, ]), n)
  td  <- td[td$BS2 == TRUE | row.names(td) %in% rn, ]
  
  
  # ranger
  bs_train <- ranger(y = td$BS2, x = td[-1],
                     probability = TRUE,
                     importance = "permutation"
  )
  
  td$pred <- predict(bs_train, td)$predictions[, 2]
  td_cm <- confusionMatrix(table(obs = td$BS2, pred = td$pred > 0.5), positive = 'TRUE')
  
  vd$pred <- predict(bs_train, vd)$predictions[, 2]
  vd_cm <- confusionMatrix(table(obs = vd$BS2, pred = vd$pred > 0.5), positive = 'TRUE')
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

summary(t(sapply(bs_unc, function(x) x$vd_cm$byClass)))
summary(t(sapply(bs_unc, function(x) x$vd_cm$overall)))
summary((sapply(bs_unc, function(x) x$AUC)))
summary((sapply(bs_unc, function(x) x$R2)))



rs2 <- rs
names(rs2) <- gsub("_CONUS", "", names(rs2))
idx <- which(
  names(rs2) %in% c("BS2", as.character(s$variable[c(1:25)]))
  # names(rs2) %in% bs_train$finalModel$forest$independent.variable.names
)
rs2 <- rs2[[idx]]

i <- 1:50
lapply(i, function(x) {
  predfun <- function(model, ...) predict(model, ...)$predictions[, 2]
  bs_r <- predict(rs2, bs_unc[[x]]$rf, fun = predfun, index = 1, progress = "text", overwrite = TRUE, filename = paste0("D:/geodata/project_data/gsp-bs/bs_dsm_", x, ".tif"))
})

fp <- "D:/geodata/project_data/gsp-bs"
lf <- list.files(fp)
lf <- lf[grepl("bs_dsm_[1-9]", lf)]
bs_rs <- stack(file.path(fp, lf))
bs_mean <- calc(bs_rs, mean, filename = file.path(fp, "bs_dsm_mean_v2.tif"), progress = "text", overwrite = TRUE)
bs_sd   <- calc(bs_rs, sd, filename = file.path(fp, "bs_dsm_sd.tif"), progress = "text", overwrite = TRUE)
bs_pi   <- calc(bs_rs, function(x) {
  quantile(x, probs = 0.975, na.rm = TRUE) - quantile(x, probs = 0.025, na.rm = TRUE)
  }, 
  filename = file.path(fp, "bs_dsm_pi.tif"), progress = "text", overwrite = TRUE
  )



# compute final confusion matrix ----
test <- extract(bs_mean, st_coordinates(bs))
ssurgo <- extract(raster(file.path(fp, "gnatsgo_fy20_1km_pct_bs.tif")), st_coordinates(bs))
test <- data.frame(BS2 = bs$BS2, dsm = test, ssurgo = ssurgo)

cm_dsm    <- confusionMatrix(table(pred = test$dsm    > 0.5, obs = test$BS2), positive = "TRUE")
cm_ssurgo <- confusionMatrix(table(pred = test$ssurgo > 0.5, obs = test$BS2), positive = "TRUE")


