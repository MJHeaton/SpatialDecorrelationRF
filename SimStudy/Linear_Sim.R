## Libraries
library(tidyverse)
library(viridis)
library(spmodel)
library(ranger)
library(patchwork)
library(pdp)
source("TransformFunctions.R")

scale_to_01 <- function(vis){
  return((vis-min(vis))/(max(vis)-min(vis)))
}

## Parameters
beta0 <- 0
beta1 <- 1
beta2 <- 1
beta3 <- 0
beta4 <- 0
sigma2_x1 <- 1 # these are variance, not st dev
sigma2_x2 <- 1
sigma2_x3 <- 1
sigma2_x4 <- 1
sigma2_W <- 1
sigma2_E <- 1
phix1 <- 1.5
phix3 <- 1.5
phiw <- 1.5
rhox1 <- .01
rhox3 <- .01
rhow <- 0
n <- 1000
split_prop <- 0.8

## Spatial Settings
de_X1 <- (sigma2_x1)*(1-rhox1)
de_X3 <- (sigma2_x2)*(1-rhox3)
de_W <- sigma2_W*(1-rhow)
ie_X1 <- sigma2_x1*rhox1
ie_X3 <- sigma2_x3*rhox3
ie_W <- sigma2_W*rhow
range <- 0.5*sqrt(2)/-log(0.05)

x1_cov_params <- spcov_params("exponential",
                             de=de_X1,
                             ie=ie_X1,
                             range=range)
x3_cov_params <- spcov_params("exponential",
                             de=de_X3,
                             ie=ie_X3,
                             range=range)
w_cov_params <- spcov_params("exponential",
                             de=de_W,
                             ie=ie_W,
                             range=range)

## Function
results_list <- list()
get_vip <- function(idx){
  ## White noise covariates
  X2 <- rnorm(n,0,sqrt(sigma2_x2))
  X4 <- rnorm(n,0,sqrt(sigma2_x4))
  epsilon <- rnorm(n,0,sqrt(sigma2_E))
  
  ## Simulate Locations
  locs <- cbind(runif(n,0,1), runif(n,0,1))
  
  # Data Frame
  df <- data.frame(X2=X2,
                   X4=X4,
                   xc=locs[,1],
                   yc=locs[,2])
  
  ## Spatial Covariates
  X1 <- sprnorm(x1_cov_params,
                data=df,
                xcoord=xc,
                ycoord=yc)
  X3 <- sprnorm(x3_cov_params,
                data=df,
                xcoord=xc,
                ycoord=yc)
  W <- sprnorm(w_cov_params,
               data=df,
               xcoord=xc,
               ycoord=yc)
  
  ## View Spatial Data
  # ggplot(data=df, aes(x=xc,y=yc,color=X1)) +
  #   geom_point() +
  #   scale_color_viridis()
  
  ## Formulas
  loc_formula <- y~X1+X2+X3+X4+xc+yc
  no_loc_formula <- y~X1+X2+X3+X4
  
  ## Y
  y <- beta0+beta1*X1+beta2*X2+W+epsilon
  
  ## Full Data frame
  df <- data.frame(y=y,
                   X1=X1,
                   X2=X2,
                   X3=X3,
                   X4=X4,
                   xc=locs[,1],
                   yc=locs[,2])
  
  ## Train Test Split
  train_set_ind <- sample.int(n, size=split_prop*n)
  trainSet <- df[train_set_ind,]
  testSet <- df[-train_set_ind,]
  
  ## Ind RF Without Locs
  ind_RF_no_locs <- ranger(no_loc_formula, data=trainSet, importance="impurity")
  ind_RF_no_locs$variable.importance <- scale_to_01(ind_RF_no_locs$variable.importance)
  
  # pdp plots
  # features <- c('X1','X2')
  # ind_RF_no_locs_plots <- lapply(features, function(f) {
  #   partial(ind_RF_no_locs, pred.var = f) %>%
  #     autoplot()
  # })
  # wrap_plots(ind_RF_no_locs_plots) +
  #   plot_annotation(title = "Ind RF No Locs PDP")
  
  ind_RF_no_locs_preds <- predict(ind_RF_no_locs, data=testSet)$predictions
  ind_RF_no_locs_RMSE <- sqrt(mean((testSet$y-ind_RF_no_locs_preds)^2))
  ind_RF_no_locs_R2 <- cor(testSet$y,ind_RF_no_locs_preds)^2
  ind_RF_no_locs_X1 <- importance(ind_RF_no_locs)[1]
  ind_RF_no_locs_X2 <- importance(ind_RF_no_locs)[2]
  ind_RF_no_locs_X3 <- importance(ind_RF_no_locs)[3]
  ind_RF_no_locs_X4 <- importance(ind_RF_no_locs)[4]
  ind_RF_no_locs_xc <- NA
  ind_RF_no_locs_yc <- NA
  
  ## Ind RF With Locs
  ind_RF_locs <- ranger(loc_formula, data=trainSet, importance="impurity")
  ind_RF_locs$variable.importance <- scale_to_01(ind_RF_locs$variable.importance)
  
  ind_RF_locs_preds <- predict(ind_RF_locs, data=testSet)$predictions
  ind_RF_locs_RMSE <- sqrt(mean((testSet$y-ind_RF_locs_preds)^2))
  ind_RF_locs_R2 <- cor(testSet$y,ind_RF_locs_preds)^2
  ind_RF_locs_X1 <- importance(ind_RF_locs)[1]
  ind_RF_locs_X2 <- importance(ind_RF_locs)[2]
  ind_RF_locs_X3 <- importance(ind_RF_locs)[3]
  ind_RF_locs_X4 <- importance(ind_RF_locs)[4]
  ind_RF_locs_xc <- importance(ind_RF_locs)[5]
  ind_RF_locs_yc <- importance(ind_RF_locs)[6]
  
  ## spRF Without Locs
  sp_RF_no_locs <- splmRF(no_loc_formula, data=trainSet, xcoord=xc, ycoord=yc, spcov_type = "exponential")
  sp_RF_no_locs_preds <- predict(sp_RF_no_locs, newdata=testSet)
  spRF_no_locs_RMSE <- sqrt(mean((testSet$y-sp_RF_no_locs_preds)^2))
  spRF_no_locs_R2 <- cor(testSet$y,sp_RF_no_locs_preds)^2
  
  ## spRF With Locs
  sp_RF_locs <- splmRF(loc_formula, data=trainSet, xcoord=xc, ycoord=yc, spcov_type = "exponential")
  sp_RF_locs_preds <- predict(sp_RF_locs, newdata=testSet)
  spRF_locs_RMSE <- sqrt(mean((testSet$y-sp_RF_locs_preds)^2))
  spRF_locs_R2 <- cor(testSet$y,sp_RF_locs_preds)^2
  
  ## Decorr Without Locs
  tune_grid <- expand.grid(range=-seq(0.001, sqrt(2), length=10)/log(0.05),
                           nugget=seq(0.05, 0.95, length=10))
  get_RMSE_no_locs <- function(idx){
    
    indData <- transform_to_ind(no_loc_formula,
                                trainData=trainSet,
                                trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                testData=testSet,
                                testLocs=as.matrix(testSet[,c('xc','yc')]),
                                M=30,
                                MaternParams=as.numeric(tune_grid[idx,]),
                                ncores=1)
    
    ## Fit RF and Backtransform
    spatial_RF <- ranger(y~., data=indData$trainData)
    spatial_preds <- predict(spatial_RF, data=indData$testData)$predictions %>%
      back_transform_to_spatial(indData)
    
    ## Calculate and return RMSE
    return(sqrt(mean((spatial_preds-testSet$y)^2)))
    
  }
  tune_rmse_no_locs <- sapply(1:nrow(tune_grid), FUN=get_RMSE_no_locs)
  decorr_range_no_locs <- tune_grid$range[which.min(tune_rmse_no_locs)]
  decorr_nugget_no_locs <- tune_grid$nugget[which.min(tune_rmse_no_locs)]
  decorr_data_no_locs <- transform_to_ind(no_loc_formula,
                                          trainData=trainSet,
                                          trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                          testData=testSet,
                                          testLocs=as.matrix(testSet[,c('xc','yc')]),
                                          M=30,
                                          MaternParams=c(decorr_range_no_locs, decorr_nugget_no_locs),
                                          ncores=1)
  deccor_RF_no_locs <- ranger(y~., data=decorr_data_no_locs$trainData, importance='impurity')
  deccor_RF_no_locs$variable.importance <- scale_to_01(deccor_RF_no_locs$variable.importance)
  decorr_RF_no_locs_preds <- predict(deccor_RF_no_locs, data=decorr_data_no_locs$testData)$predictions %>%
    back_transform_to_spatial(decorr_data_no_locs)
  decorr_no_locs_RMSE <- sqrt(mean((testSet$y-decorr_RF_no_locs_preds)^2))
  decorr_no_locs_R2 <- cor(testSet$y,decorr_RF_no_locs_preds)^2
  decorr_no_locs_X1 <- importance(deccor_RF_no_locs)[2]
  decorr_no_locs_X2 <- importance(deccor_RF_no_locs)[3]
  decorr_no_locs_X3 <- importance(deccor_RF_no_locs)[4]
  decorr_no_locs_X4 <- importance(deccor_RF_no_locs)[5]
  decorr_no_locs_xc <- NA
  decorr_no_locs_yc <- NA
  
  ## Decorr With Locs
  get_RMSE_locs <- function(idx){
    
    indData <- transform_to_ind(loc_formula,
                                trainData=trainSet,
                                trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                testData=testSet,
                                testLocs=as.matrix(testSet[,c('xc','yc')]),
                                M=30,
                                MaternParams=as.numeric(tune_grid[idx,]),
                                ncores=1)
    
    ## Fit RF and Backtransform
    spatial_RF <- ranger(y~., data=indData$trainData)
    spatial_preds <- predict(spatial_RF, data=indData$testData)$predictions %>%
      back_transform_to_spatial(indData)
    
    ## Calculate and return RMSE
    return(sqrt(mean((spatial_preds-testSet$y)^2)))
    
  }
  tune_rmse_locs <- sapply(1:nrow(tune_grid), FUN=get_RMSE_locs)
  decorr_range_locs <- tune_grid$range[which.min(tune_rmse_locs)]
  decorr_nugget_locs <- tune_grid$nugget[which.min(tune_rmse_locs)]
  decorr_data_locs <- transform_to_ind(loc_formula,
                                       trainData=trainSet,
                                       trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                       testData=testSet,
                                       testLocs=as.matrix(testSet[,c('xc','yc')]),
                                       M=30,
                                       MaternParams=c(decorr_range_locs, decorr_nugget_locs),
                                       ncores=1)
  deccor_RF_locs <- ranger(y~., data=decorr_data_locs$trainData, importance='impurity')
  deccor_RF_locs$variable.importance <- scale_to_01(deccor_RF_locs$variable.importance)
  decorr_RF_locs_preds <- predict(deccor_RF_locs, data=decorr_data_locs$testData)$predictions %>%
    back_transform_to_spatial(decorr_data_locs)
  decorr_locs_RMSE <- sqrt(mean((testSet$y-decorr_RF_locs_preds)^2))
  decorr_locs_R2 <- cor(testSet$y,decorr_RF_locs_preds)^2
  decorr_locs_X1 <- importance(deccor_RF_locs)[2]
  decorr_locs_X2 <- importance(deccor_RF_locs)[3]
  decorr_locs_X3 <- importance(deccor_RF_locs)[4]
  decorr_locs_X4 <- importance(deccor_RF_locs)[5]
  decorr_locs_xc <- importance(deccor_RF_locs)[6]
  decorr_locs_yc <- importance(deccor_RF_locs)[7]
  
  # Put Results in Data Frame
  results <- data.frame(Method=c("Ind RF Loc", "Ind RF No Loc", "spRF Loc", "spRF No Loc", "Decorr Loc", "Decorr No Loc"),
                        RMSE=c(ind_RF_locs_RMSE, ind_RF_no_locs_RMSE, spRF_locs_RMSE, spRF_no_locs_RMSE, decorr_locs_RMSE, decorr_no_locs_RMSE),
                        R2=c(ind_RF_locs_R2, ind_RF_no_locs_R2, spRF_locs_R2, spRF_no_locs_R2, decorr_locs_R2, decorr_no_locs_R2),
                        X1=c(ind_RF_locs_X1, ind_RF_no_locs_X1, ind_RF_locs_X1, ind_RF_no_locs_X1, decorr_locs_X1, decorr_no_locs_X1),
                        X2=c(ind_RF_locs_X2, ind_RF_no_locs_X2, ind_RF_locs_X2, ind_RF_no_locs_X2, decorr_locs_X2, decorr_no_locs_X2),
                        X3=c(ind_RF_locs_X3, ind_RF_no_locs_X3, ind_RF_locs_X3, ind_RF_no_locs_X3, decorr_locs_X3, decorr_no_locs_X3),
                        X4=c(ind_RF_locs_X4, ind_RF_no_locs_X4, ind_RF_locs_X4, ind_RF_no_locs_X4, decorr_locs_X4, decorr_no_locs_X4),
                        xc=c(ind_RF_locs_xc, ind_RF_no_locs_xc, ind_RF_locs_xc, ind_RF_no_locs_xc, decorr_locs_xc, decorr_no_locs_xc),
                        yc=c(ind_RF_locs_yc, ind_RF_no_locs_yc, ind_RF_locs_yc, ind_RF_no_locs_yc, decorr_locs_yc, decorr_no_locs_yc))
  return(results)
}
n_iters <- 3
results_list <- mclapply(1:100, FUN=get_vip, mc.cores=5)

final_results <- bind_rows(results_list) %>%
  mutate(Method = factor(Method, levels = unique(Method))) %>% # preserve order
  group_by(Method) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))  # Compute mean for numeric columns

# write.csv(final_results, 'SimulateDataV3.csv')



# pdp plots
features <- c('X1', 'X2')
ind_RF_locs_plots <- lapply(features, function(f) {
  partial(ind_RF_locs, pred.var = f) %>%
    autoplot()
})
ind_RF_no_locs_plots <- lapply(features, function(f) {
  partial(ind_RF_no_locs, pred.var = f) %>%
    autoplot()
})
sp_RF_locs_plots <- lapply(features, function(f) {
  partial(sp_RF_locs, pred.var = f) %>%
    autoplot()
})
sp_RF_no_locs_plots <- lapply(features, function(f) {
  partial(sp_RF_no_locs, pred.var = f) %>%
    autoplot()
})
decorr_RF_locs_plots <- lapply(features, function(f) {
  partial(deccor_RF_locs, pred.var = f) %>%
    autoplot()
})
decorr_RF_no_locs_plots <- lapply(features, function(f) {
  partial(deccor_RF_no_locs, pred.var = f) %>%
    autoplot()
})

wrap_plots(c(ind_RF_locs_plots,
             sp_RF_locs_plots,
             decorr_RF_locs_plots)) +
  plot_annotation(
    title = "Partial Dependence Plots",
    tag_levels = "A"  # adds row/col labels automatically
  ) &
  theme(plot.tag.position = "left") 

wrap_plots(ind_RF_locs_plots)/
  wrap_plots(sp_RF_locs_plots)/
  wrap_plots(decorr_RF_locs_plots)




