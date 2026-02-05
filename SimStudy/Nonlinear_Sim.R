##
## Simulate Data for Spatial ML Study
##

## Libraries
library(tidyverse)
library(viridis)
library(spmodel)
library(ranger)
library(pdp)
library(patchwork)
scale_to_01 <- function(vis){
  return((vis-min(vis))/(max(vis)-min(vis)))
}
get_perm_pdp_decorr <- function(){
  ftr_vals <- trainSet[["X1"]]
  ftr_seq <- seq(min(ftr_vals), max(ftr_vals), length=nrow(testSet))
  for(x2val in c(0,1)){
    pdp_val <- rep(NA, length(ftr_seq))
    for(v in ftr_seq){
      data_copy <- testSet
      data_copy[["X1"]] <- v
      data_copy$X2 <- x2val
      ind_data_copy <- transform_to_ind(our_formula,
                                        trainData=trainSet,
                                        trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                        testData=data_copy,
                                        testLocs=as.matrix(data_copy[,c('xc','yc')]),
                                        M=30,
                                        MaternParams=c(decorr_range, decorr_nugget),
                                        ncores=1)
      pdp_val[v==ftr_seq] <- mean(predict(deccor_RF, data=ind_data_copy$testData)$predictions %>%
                                    back_transform_to_spatial(ind_data_copy))
    }
    if(x2val==0){
      pdp_out <- data.frame(X1=ftr_seq, PDP=pdp_val, X2=x2val)
    } else {
      pdp_out <- bind_rows(pdp_out,
                           data.frame(X1=ftr_seq, PDP=pdp_val, X2=x2val))
    }
  }
  return(pdp_out)
}
get_perm_pdp_ind <- function(){
  ftr_vals <- trainSet[["X1"]]
  ftr_seq <- seq(min(ftr_vals), max(ftr_vals), length=nrow(testSet))
  for(x2val in c(0,1)){
    pdp_val <- rep(NA, length(ftr_seq))
    for(v in ftr_seq){
      data_copy <- testSet
      data_copy[["X1"]] <- v
      data_copy$X2 <- x2val
      pdp_val[v==ftr_seq] <- mean(predict(ind_RF, data=data_copy)$predictions)
    }
    if(x2val==0){
      pdp_out <- data.frame(X1=ftr_seq, PDP=pdp_val, X2=x2val)
    } else {
      pdp_out <- bind_rows(pdp_out,
                           data.frame(X1=ftr_seq, PDP=pdp_val, X2=x2val))
    }
  }
  return(pdp_out)
}
source("TransformFunctions.R")

## Settings
sample_size <- 1000
split_prop <- 0.8

## Spatial Settings
de <- 10*0.95
de_X <- 10*0.5
ie <- 10*0.05
ie_X <- 0#10*0.5
true_range <- 0.5*sqrt(2)/-log(0.05)

## Formula
our_formula <- y~X1+X2+X3+X4+X5

## Covariates
X1 <- rnorm(sample_size, 0, 1)
X2 <- rbinom(sample_size, size=1,prob=0.5)
X4 <- rnorm(sample_size)
X5 <- rbinom(sample_size, size=1,prob=0.5)

## Spatial Locations
xc <- runif(sample_size)
yc <- runif(sample_size)

## Data frame
df <- data.frame(X1=X1,
                 X2=X2,
                 X4=X4,
                 X5=X5,
                 xc=xc,
                 yc=yc)

## Generate Spatial Covariate
x_cov_params <- spcov_params("exponential",
                             de=de_X,
                             ie=ie_X,
                             range=true_range)
X3 <- sprnorm(x_cov_params,
              data=df,
              xcoord=xc,
              ycoord=yc)
df$X3 <- X3
ggplot(data=df, aes(x=xc,y=yc,color=X3)) +
  geom_point() +
  scale_color_viridis()

## Data frame
df <- data.frame(X1=X1,
                 X2=X2,
                 X3=X3,
                 X4=X4,
                 X5=X5,
                 xc=xc,
                 yc=yc)

## Mean Structure
mn <- sin(X1) + X2 * cos(X1) + 0.1 * X1^2
ggplot() + geom_point(aes(x=xc, y=yc, color=mn)) +
  scale_color_viridis()

## Simulate Error
cov_params <- spcov_params("exponential",
                           de=de,
                           ie=ie,
                           range=true_range)
sp_error <- sprnorm(spcov_params=cov_params, data=df,
        xcoord=xc, ycoord=yc)
ggplot() + geom_point(aes(x=xc, y=yc, color=sp_error)) +
  scale_color_viridis()

## Response
df$y <- mn+sp_error
ggplot(data=df, aes(x=xc, y=yc, color=y)) + geom_point() +
  scale_color_viridis()

## Split test-train
train_set_ind <- sample.int(sample_size, size=split_prop*sample_size)
trainSet <- df[train_set_ind,]
testSet <- df[-train_set_ind,]

x1 <- seq(min(testSet$X1),max(testSet$X1), length.out = 500)
x2 <- c(0, 1)
grid <- expand.grid(X1 = x1, X2 = x2)
grid <- bind_cols(grid, X3=X3, X4=X4, X5=X5, xc=xc, yc=yc)
grid$Z <- with(grid, sin(X1) + X2 * cos(X1) + 0.1 * X1^2)
ggplot(grid, aes(x = X1, y = Z, color = factor(X2))) +
  geom_line(linewidth = 1) +
  labs(x = expression(X[1]), y = expression(f(X[1],X[2])), color = expression(X[2])) +
  theme_minimal()

# list of variable names
features <- setdiff(names(df), c('xc', 'yc', 'y'))

## Fit IND RF
ind_RF <- ranger(our_formula, data=trainSet, importance="impurity")
ind_RF$variable.importance <- scale_to_01(ind_RF$variable.importance)
vip::vip(ind_RF)
# pdp plots
ind_RF_plots <- lapply(features, function(f) {
  partial(ind_RF, pred.var = f) %>%
    autoplot()
})
wrap_plots(ind_RF_plots) + 
  plot_annotation(title = "Ind RF Partial Dependence Plots")
ind_preds <- predict(ind_RF, data=testSet)$predictions
ind_RMSE <- sqrt(mean((testSet$y-ind_preds)^2))
ind_R2 <- cor(testSet$y,ind_preds)^2
ind_X1 <- importance(ind_RF)[1]
ind_X2 <- importance(ind_RF)[2]
ind_X3 <- importance(ind_RF)[3]
ind_X4 <- importance(ind_RF)[4]
ind_X5 <- importance(ind_RF)[5]
ind_X6 <- importance(ind_RF)[6]
ind_preds <- bind_cols(testSet,Ind=ind_preds) %>%
  arrange(X1)
ggplot() +
  geom_line(data=grid, 
            mapping=aes(x = X1, y = Z, color = factor(X2))) +
  labs(x = expression(X[1]), y = expression(f(X[1],X[2])), color = expression(X[2])) +
  theme_minimal() +
  geom_point(data=ind_preds, aes(x=X1, y=Ind, color=factor(X2)))
  

## Fit spRF
sp_RF <- splmRF(our_formula, data=trainSet, xcoord=xc, ycoord=yc)
sp_RF_preds <- predict(sp_RF, newdata=testSet)
spRF_RMSE <- sqrt(mean((testSet$y-sp_RF_preds)^2))
spRF_R2 <- cor(testSet$y,sp_RF_preds)^2
# pdp plots
sp_RF_plots <- lapply(features, function(f) {
  partial(sp_RF, pred.var = f, type='regression') %>%
    autoplot()
})
wrap_plots(sp_RF_plots) + 
  plot_annotation(title = "Sp RF Partial Dependence Plots")

# Compute importance scores
spRF_X1 <- ind_X1
spRF_X2 <- ind_X2
spRF_X3 <- ind_X3
spRF_X4 <- ind_X4
spRF_X5 <- ind_X5
spRF_X6 <- ind_X6

## Fit Decorr
tune_grid <- expand.grid(range=seq(0.001, sqrt(2), length=10),
                         nugget=seq(0.05, 0.95, length=10))
get_RMSE <- function(idx){
  
  indData <- transform_to_ind(our_formula,
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
tune_rmse <- sapply(1:nrow(tune_grid), FUN=get_RMSE)
decorr_range <- tune_grid$range[which.min(tune_rmse)]
decorr_nugget <- tune_grid$nugget[which.min(tune_rmse)]
decorr_data <- transform_to_ind(our_formula,
                 trainData=trainSet,
                 trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                 testData=testSet,
                 testLocs=as.matrix(testSet[,c('xc','yc')]),
                 M=30,
                 MaternParams=c(decorr_range, decorr_nugget),
                 ncores=1)
deccor_RF <- ranger(y~., data=decorr_data$trainData, importance='impurity')
deccor_RF$variable.importance <- scale_to_01(deccor_RF$variable.importance)
vip::vip(deccor_RF)
# pdp plots
decorr_plots <- lapply(features, function(f) {
  partial(deccor_RF, pred.var = f) %>%
    autoplot()
})
wrap_plots(decorr_plots) + 
  plot_annotation(title = "Decorr Partial Dependence Plots")
decorr_RF_preds <- predict(deccor_RF, data=decorr_data$testData)$predictions %>%
  back_transform_to_spatial(decorr_data)
decorr_RMSE <- sqrt(mean((testSet$y-decorr_RF_preds)^2))
decorr_R2 <- cor(testSet$y,decorr_RF_preds)^2
decorr_X1 <- importance(deccor_RF)[2]
decorr_X2 <- importance(deccor_RF)[3]
decorr_X3 <- importance(deccor_RF)[4]
decorr_X4 <- importance(deccor_RF)[5]
decorr_X5 <- importance(deccor_RF)[6]
decorr_grid <- transform_to_ind(our_formula,
                                trainData=trainSet,
                                trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                testData=grid,
                                testLocs=as.matrix(grid[,c('xc','yc')]),
                                M=30,
                                MaternParams=c(decorr_range, decorr_nugget),
                                ncores=1)
decorr_grid <- bind_cols(grid,
                             P=predict(deccor_RF, data=decorr_grid$testData)$predictions) %>%
  arrange(X1)
ggplot() +
  geom_line(data=grid, 
            mapping=aes(x = X1, y = Z, color = factor(X2))) +
  labs(x = expression(X[1]), y = expression(f(X[1],X[2])), color = expression(X[2])) +
  theme_minimal() +
  geom_point(data=decorr_grid, aes(x=X1, y=P, color=factor(X2)))

## Compare PDPs between features
decorr_RF_pdp <- get_perm_pdp_decorr()
ind_RF_pdp <- get_perm_pdp_ind()
PDP_plots <- bind_rows(decorRF=decorr_RF_pdp,
                       indRF=ind_RF_pdp,
                       .id="Method") %>%
  group_by(X2, Method) %>%
  mutate(PDP=PDP-mean(PDP))
x1pdp <- ggplot() +
  geom_line(data=PDP_plots, aes(x=X1, y=PDP, 
                                color=Method, 
                                linetype = factor(X2)),
            linewidth=0.75) +
  labs(linetype = expression(X[2])) +
  ylab(expression(hat(y)))

## Compare PDP for X3
get_X3_pdp_ind <- function(){
  ftr_vals <- trainSet[["X3"]]
  ftr_seq <- seq(min(ftr_vals), max(ftr_vals), length=nrow(testSet))
  pdp_val <- rep(NA, length(ftr_seq))
  for(v in ftr_seq){
    data_copy <- testSet
    data_copy[["X3"]] <- v
    pdp_val[v==ftr_seq] <- mean(predict(ind_RF, data=data_copy)$predictions)
  }
  pdp_out <- data.frame(X1=ftr_seq, PDP=pdp_val)
  return(pdp_out)
}
get_X3_pdp_decorr <- function(){
  ftr_vals <- trainSet[["X3"]]
  ftr_seq <- seq(min(ftr_vals), max(ftr_vals), length=nrow(testSet))
  pdp_val <- rep(NA, length(ftr_seq))
  for(v in ftr_seq){
    data_copy <- testSet
    data_copy[["X3"]] <- v
    ind_data_copy <- transform_to_ind(our_formula,
                                      trainData=trainSet,
                                      trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                      testData=data_copy,
                                      testLocs=as.matrix(data_copy[,c('xc','yc')]),
                                      M=30,
                                      MaternParams=c(decorr_range, decorr_nugget),
                                      ncores=1)
    pdp_val[v==ftr_seq] <- mean(predict(deccor_RF, data=ind_data_copy$testData)$predictions %>%
                                  back_transform_to_spatial(ind_data_copy))
  }
  pdp_out <- data.frame(X1=ftr_seq, PDP=pdp_val)
  return(pdp_out)
}
X3_pdp_plots <- bind_rows(decorRF=get_X3_pdp_decorr(),
                          indRF=get_X3_pdp_ind(),
                          .id="Method") %>%
  group_by(Method) %>%
  mutate(PDP=PDP-mean(PDP))
x3pdp <- ggplot(data=X3_pdp_plots, aes(x=X1, y=PDP, color=Method)) +
  geom_line() +
  xlab(expression(X[3])) +
  ylab(expression(hat(y)))
x1pdp + x3pdp

## Compare PDPs evaluated on original vs transformed
get_perm_pdp_decorr_os_margX2 <- function(){
  ftr_vals <- trainSet[["X1"]]
  ftr_seq <- seq(min(ftr_vals), max(ftr_vals), length=nrow(testSet))
  pdp_val <- rep(NA, length(ftr_seq))
  for(v in ftr_seq){
    data_copy <- testSet
    data_copy[["X1"]] <- v
    ind_data_copy <- transform_to_ind(our_formula,
                                      trainData=trainSet,
                                      trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                      testData=data_copy,
                                      testLocs=as.matrix(data_copy[,c('xc','yc')]),
                                      M=30,
                                      MaternParams=c(decorr_range, decorr_nugget),
                                      ncores=1)
    pdp_val[v==ftr_seq] <- mean(predict(deccor_RF, data=ind_data_copy$testData)$predictions %>%
                                  back_transform_to_spatial(ind_data_copy))
  }
  pdp_out <- data.frame(X1=ftr_seq, PDP=pdp_val)
  return(pdp_out)
}
get_perm_pdp_decorr_ts_margX2 <- function(){
  ind_data_copy <- transform_to_ind(our_formula,
                                    trainData=trainSet,
                                    trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                    testData=testSet,
                                    testLocs=as.matrix(testSet[,c('xc','yc')]),
                                    M=30,
                                    MaternParams=c(decorr_range, decorr_nugget),
                                    ncores=1)
  ftr_vals <- ind_data_copy$trainData[["X1"]]
  ftr_seq <- seq(min(ftr_vals), max(ftr_vals), length=nrow(testSet))
  pdp_val <- rep(NA, length(ftr_seq))
  for(v in ftr_seq){
    ind_data_copy$testData[["X1"]] <- v
    pdp_val[v==ftr_seq] <- mean(predict(deccor_RF, data=ind_data_copy$testData)$predictions %>%
                                  back_transform_to_spatial(ind_data_copy))
  }
  pdp_out <- data.frame(X1=ftr_seq, PDP=pdp_val)
  return(pdp_out)
}
pdp_comp <- bind_rows(
  bind_cols(data.frame(Type=rep("OS", nrow(testSet))),get_perm_pdp_decorr_os_margX2()),
  bind_cols(data.frame(Type=rep("TS", nrow(testSet))),get_perm_pdp_decorr_ts_margX2()))
 ggplot(data=pdp_comp, aes(x=X1, y=PDP, color=Type)) +
   geom_line()
  

## Fit decorr_spRF
decorr_spRF_range <- coef(sp_RF$splm, type="spcov")[3]
decorr_spRF_nugget <- coef(sp_RF$splm, type="spcov")[2]/
  sum(coef(sp_RF$splm, type="spcov")[1:2])
decorr_spRF_data <- transform_to_ind(our_formula,
                                     trainData=trainSet,
                                     trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                     testData=testSet,
                                     testLocs=as.matrix(testSet[,c('xc','yc')]),
                                     M=30,
                                     MaternParams=c(decorr_spRF_range, decorr_spRF_nugget),
                                     ncores=1)
deccor_spRF_RF <- ranger(y~., data=decorr_spRF_data$trainData, importance='impurity')
deccor_spRF_RF$variable.importance <- scale_to_01(deccor_spRF_RF$variable.importance)
vip::vip(deccor_spRF_RF)
# pdp plots
deccor_spRF_RF_plots <- lapply(features, function(f) {
  partial(deccor_spRF_RF, pred.var = f) %>%
    autoplot()
})
wrap_plots(deccor_spRF_RF_plots) + 
  plot_annotation(title = "Decorr SpRF Partial Dependence Plots")
decorr_spRF_RF_preds <- predict(deccor_spRF_RF, data=decorr_spRF_data$testData)$predictions %>%
  back_transform_to_spatial(decorr_spRF_data)
decorr_spRF_RMSE <- sqrt(mean((testSet$y-decorr_spRF_RF_preds)^2))
decorr_spRF_R2 <- cor(testSet$y,decorr_spRF_RF_preds)^2
decorr_spRF_X1 <- importance(deccor_spRF_RF)[2]
decorr_spRF_X2 <- importance(deccor_spRF_RF)[3]
decorr_spRF_X3 <- importance(deccor_spRF_RF)[4]
decorr_spRF_X4 <- importance(deccor_spRF_RF)[5]
decorr_spRF_X5 <- importance(deccor_spRF_RF)[6]

## Decorrelate Y's but not X's
trainLocs <- trainSet[,c('xc','yc')]
testLocs <- testSet[,c('xc','yc')]
y_train <- trainSet$y
X_train <- trainSet[,c('X1', 'X2', 'X3', 'X4', 'X5')]
y_test <- testSet$y
X_test <- testSet[,c('X1', 'X2', 'X3', 'X4', 'X5')]
K <- 30
tune_grid <- expand.grid(range=seq(0.001, sqrt(2), length=10),
                         nugget=seq(0.05, 0.95, length=10))
get_rmse <- function(idx){
  phi <- tune_grid[idx,1]
  omega <- tune_grid[idx,2]
  
  # decorrelate training data
  y_uncor <- matrix(NA, nrow = nrow(trainSet), ncol = 1)
  y_uncor[1] <- y_train[1]
  for(i in 2:nrow(trainSet)){
    distances <- fields::rdist(trainLocs[i, ,drop=F],trainLocs[1:(i-1),])
    neighbors <- order(distances)[1:min(i-1,K)]
    D <- fields::rdist(trainLocs[c(i,neighbors),])
    R <- (1-omega)*fields::Exponential(D, aRange=phi)+omega*diag(nrow(D))
    R_inv <- solve(R[-1, -1])
    denom <- as.numeric(1-(R[1,-1] %*% R_inv %*% R[-1,1]))
    y_uncor[i] <- y_train[i] - R[1,-1] %*% R_inv %*% y_train[neighbors]
    y_uncor[i] <- y_uncor[i]/sqrt(denom)
  }
  
  # fit model
  data_ind <- data.frame(y_uncor = y_uncor, X_train)
  
  rf <- ranger(y_uncor ~ ., data = data_ind)
  preds_ind <- predict(rf, data=X_test)$predictions
  
  # backtransform predictions
  preds <- matrix(NA, nrow=nrow(testSet),ncol=1)
  for (i in 1:length(preds_ind)){
    distances <- fields::rdist(testLocs[i, ,drop=F], trainLocs)
    neighbors <- order(distances)[1:K]
    D <- fields::rdist(rbind(testLocs[i,], trainLocs[neighbors,]))
    R <- (1-omega)*fields::Exponential(D, aRange=phi)+omega*diag(nrow(D))
    R_inv <- solve(R[-1, -1])
    mult <- sqrt(as.numeric(1-R[1,-1] %*% R_inv %*% t(t(R[1,-1]))))
    preds[i] <- preds_ind[i]*mult + R[1,-1] %*% R_inv %*% y_train[neighbors]
  }
  
  # find and return rmse
  predictions <- data.frame(truth = testSet$y,
                            estimate = preds)
  rf_rmse <- yardstick::rmse(data=predictions, truth=truth, estimate=estimate)$.estimate
  
  return(rf_rmse)
}

# save tuned parameters
tune_rmse <- mclapply(1:nrow(tune_grid), FUN=get_rmse, mc.cores=1)
best_range <- tune_grid$range[which.min(tune_rmse)]
best_nugget <- tune_grid$nugget[which.min(tune_rmse)]
phi <- best_range
omega <- best_nugget
# now with best hyperparameters
# decorrelate training data
y_uncor <- matrix(NA, nrow = nrow(trainSet), ncol = 1)
y_uncor[1] <- y_train[1]
for(i in 2:nrow(trainSet)){
  distances <- fields::rdist(trainLocs[i, ,drop=F],trainLocs[1:(i-1),])
  neighbors <- order(distances)[1:min(i-1,K)]
  D <- fields::rdist(trainLocs[c(i,neighbors),])
  R <- (1-omega)*fields::Exponential(D, aRange=phi)+omega*diag(nrow(D))
  R_inv <- solve(R[-1, -1])
  denom <- as.numeric(1-(R[1,-1] %*% R_inv %*% R[-1,1]))
  y_uncor[i] <- y_train[i] - R[1,-1] %*% R_inv %*% y_train[neighbors]
  y_uncor[i] <- y_uncor[i]/sqrt(denom)
}

# fit model
data_ind <- data.frame(y_uncor = y_uncor, X_train)

decorr_y_RF <- ranger(y_uncor ~ ., data = data_ind, importance = 'impurity')
decorr_y_RF$variable.importance <- scale_to_01(decorr_y_RF$variable.importance)
preds_ind <- predict(decorr_y_RF, data=X_test)$predictions

# backtransform predictions
preds <- matrix(NA, nrow=nrow(testSet),ncol=1)
for (i in 1:length(preds_ind)){
  distances <- fields::rdist(testLocs[i, ,drop=F], trainLocs)
  neighbors <- order(distances)[1:K]
  D <- fields::rdist(rbind(testLocs[i,], trainLocs[neighbors,]))
  R <- (1-omega)*fields::Exponential(D, aRange=phi)+omega*diag(nrow(D))
  R_inv <- solve(R[-1, -1])
  mult <- sqrt(as.numeric(1-R[1,-1] %*% R_inv %*% t(t(R[1,-1]))))
  preds[i] <- preds_ind[i]*mult + R[1,-1] %*% R_inv %*% y_train[neighbors]
}

vip::vip(decorr_y_RF)
# pdp plots
decorr_y_plots <- lapply(features, function(f) {
  partial(decorr_y_RF, pred.var = f) %>%
    autoplot()
})
wrap_plots(decorr_y_plots) + 
  plot_annotation(title = "Y Decorr Partial Dependence Plots")
decorr_y_RF_preds <- preds
decorr_y_RMSE <- sqrt(mean((testSet$y-decorr_y_RF_preds)^2))
decorr_y_R2 <- cor(testSet$y,decorr_y_RF_preds)^2
decorr_y_X1 <- importance(decorr_y_RF)[1]
decorr_y_X2 <- importance(decorr_y_RF)[2]
decorr_y_X3 <- importance(decorr_y_RF)[3]
decorr_y_X4 <- importance(decorr_y_RF)[4]
decorr_y_X5 <- importance(decorr_y_RF)[5]




## Only Use Intercept to Estimate Range and Nugget
base_mod <- splm(y~1, spcov_type = 'exponential', data=trainSet, xcoord=xc, ycoord=yc)
base_mod_range <- coef(base_mod$splm, type="spcov")[3]
base_mod_nugget <- coef(base_mod$splm, type="spcov")[2]/
  sum(coef(sp_RF$splm, type="spcov")[1:2])
base_mod_data <- transform_to_ind(our_formula,
                                     trainData=trainSet,
                                     trainLocs=as.matrix(trainSet[,c('xc','yc')]),
                                     testData=testSet,
                                     testLocs=as.matrix(testSet[,c('xc','yc')]),
                                     M=30,
                                     MaternParams=c(decorr_spRF_range, decorr_spRF_nugget),
                                     ncores=1)
base_mod_RF <- ranger(y~., data=base_mod_data$trainData, importance='impurity')
base_mod_RF$variable.importance <- scale_to_01(base_mod_RF$variable.importance)
vip::vip(base_mod_RF)
# pdp plots
base_mod_plots <- lapply(features, function(f) {
  partial(base_mod_RF, pred.var = f) %>%
    autoplot()
})
wrap_plots(base_mod_plots) + 
  plot_annotation(title = "Base Mod Partial Dependence Plots")
base_mod_RF_preds <- predict(base_mod_RF, data=base_mod_data$testData)$predictions %>%
  back_transform_to_spatial(base_mod_data)
base_mod_RMSE <- sqrt(mean((testSet$y-base_mod_RF_preds)^2))
base_mod_R2 <- cor(testSet$y,base_mod_RF_preds)^2
base_mod_X1 <- importance(base_mod_RF)[2]
base_mod_X2 <- importance(base_mod_RF)[3]
base_mod_X3 <- importance(base_mod_RF)[4]
base_mod_X4 <- importance(base_mod_RF)[5]
base_mod_X5 <- importance(base_mod_RF)[6]


## Put Results in Data Frame
results <- data.frame(Method=c("Ind", "splmRF", "Decorr", "spRF_Decorr", "y_Decorr", "Base_Mod"),
                      RMSE=c(ind_RMSE, spRF_RMSE, decorr_RMSE, decorr_spRF_RMSE, decorr_y_RMSE, base_mod_RMSE),
                      R2=c(ind_R2, spRF_R2, decorr_R2, decorr_spRF_R2, decorr_y_R2, base_mod_R2),
                      X1=c(ind_X1, spRF_X1, decorr_X1, decorr_spRF_X1, decorr_y_X1, base_mod_X1),
                      X2=c(ind_X2, spRF_X2, decorr_X2, decorr_spRF_X2, decorr_y_X2, base_mod_X2),
                      X3=c(ind_X3, spRF_X3, decorr_X3, decorr_spRF_X3, decorr_y_X3, base_mod_X3),
                      X4=c(ind_X4, spRF_X4, decorr_X4, decorr_spRF_X4, decorr_y_X4, base_mod_X4),
                      X5=c(ind_X5, spRF_X5, decorr_X5, decorr_spRF_X5, decorr_y_X5, base_mod_X5))
print(results)

