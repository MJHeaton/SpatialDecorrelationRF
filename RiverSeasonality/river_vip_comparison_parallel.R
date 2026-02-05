library(tidyverse)
library(fields)
library(ranger)
library(vroom)
library(viridis)
library(maps)
library(doParallel)
library(foreach)
library(yardstick)
library(spmodel)
library(patchwork)
library(pdp)

sp_empirical_acf <- function(formula,
                             data,
                             xcoord,
                             ycoord,
                             ...,
                             var_from = c("residual", "response"),
                             plot = TRUE) {
  # formula: e.g. z ~ 1 or z ~ covariates (esv uses residuals of this model)
  # data: data frame
  # xcoord, ycoord: quoted column names with coordinates
  # ...: extra args to esv() such as bins, cutoff, cloud
  # var_from:
  #   "residual" → use variance of residuals from lm(formula)
  #   "response" → use variance of the response directly
  # plot: TRUE returns a ggplot object
  
  var_from <- match.arg(var_from)
  
  # 1. Empirical semivariogram
  ev <- esv(
    formula = formula,
    data    = data,
    xcoord  = xcoord,
    ycoord  = ycoord,
    ...
  )
  
  # 2. Estimate total variance (sigma^2 + tau^2)
  if (var_from == "residual") {
    lm_fit <- lm(formula, data = data)
    z_vals <- residuals(lm_fit)
  } else {
    mf <- model.frame(formula, data = data)
    z_vals <- model.response(mf)
  }
  
  sigma2_hat <- var(z_vals, na.rm = TRUE)
  
  # 3. Convert semivariogram → ACF
  ev$acf <- 1 - ev$gamma / sigma2_hat
  ev$variance_hat <- sigma2_hat
  
  # 4. Plot ACF if desired
  if (plot) {
    p <- ggplot(ev, aes(x = dist, y = acf)) +
      geom_point(size = 2, alpha = 0.8) +
      geom_line(alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      labs(
        x = "Distance h",
        y = "Empirical ACF  ρ(h)",
        title = "Empirical Spatial Autocorrelation Function"
      ) +
      ylim(-1, 1) +
      theme_bw()
    
    print(p)
  } else {
    p <- NULL
  }
  
  # return both data and plot
  list(
    acf_table = ev,
    plot = p
  )
}

# choose number of cores & folds
nfolds <- 10
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

load('./PCA.RData')

## Filter to just USA and Canada
na_indices <- which(
  locs[,2] >= 25 & locs[,2] <= 83 &
    locs[,1] >= -170 & locs[,1] <= -52
)
locs <- locs[na_indices,] %>%
  drop_na()

## Use PCA value as response
y <- CombinedMetrics[,2][na_indices] 

## Get X
X <- as.matrix(CompX[,-c(9,96)][na_indices,]) # BIOME column(96) not numeric. took out gord, weird variable

## Use these variables
keep_variables <- c(#"strmOrder",
                    "area_sqkm",
                    "drain_den",
                    "MeanTempAnn",# or bio1 (2 options of temperature available)
                    "MeanPrecAnn",# or CumPrecTotal (not sure if one will be better than the other, but leaning toward using mean)
                    "cls7", # (% of watershed in agricultural land cover)
                    "cls9", # (% of watershed in urban land cover)
                    "HydroLakes_Area_sqkm", # (was important to the story in the paper already and could indicate overall hydrology of the system)
                    "meanPercentDC_Imperfectly",
                    "meanPercentDC_ModeratelyWell",
                    "meanPercentDC_Poor",
                    "meanPercentDC_SomewhatExcessive",
                    "meanPercentDC_VeryPoor",
                    "meanPercentDC_Well",
                    "Dam_Count")
X <- cbind(X[,keep_variables], clsSum=rowSums(X[,c("cls1", "cls2", "cls3", "cls4")]))

pretty_labels <- c(#"strmOrder" = "Strahler Stream Order",
                   "area_sqkm" = "Total Catchment Area",
                   "drain_den" = "Drainage Density",
                   "MeanTempAnn" = "Mean Annual Temperature",# or bio1 (2 options of temperature available)
                   "MeanPrecAnn" = "Mean Annual Precip",# or CumPrecTotal (not sure if one will be better than the other, but leaning toward using mean)
                   "cls7" = "Cultivated Vegetation", # (% of watershed in agricultural land cover)
                   "cls9" = "Urban", # (% of watershed in urban land cover)
                   "HydroLakes_Area_sqkm" = "Lake Area", # (was important to the story in the paper already and could indicate overall hydrology of the system)
                   "meanPercentDC_Imperfectly" = "Pct. Moist Soil",
                   "meanPercentDC_ModeratelyWell" = "Pct. Slightly Moist Soil",
                   "meanPercentDC_Poor" = "Pct. Wet Soil",
                   "meanPercentDC_SomewhatExcessive" = "Pct Very Dry Soil",
                   "meanPercentDC_VeryPoor" = "Pct. Very Wet Soil",
                   "meanPercentDC_Well" = "Pct. Dry Soil",
                   "Dam_Count" = "Dam Count",
                   "clsSum" = "Tree Cover")

## Plot Data
na_countries <- c("USA", "Canada", "Mexico")
na_map <- subset(map_data("world"), region %in% na_countries)

ggplot() +
  geom_polygon(data = na_map, aes(x = long, y = lat, group = group),
               fill = "white", color="gray70") +
  geom_point(mapping=aes(x=locs[,1], y=locs[,2], color=y), size=0.5) +
  scale_color_viridis() +
  coord_quickmap(xlim = range(locs[,1]), ylim = range(locs[,2])) +
  labs(color=expression(PC[2]), x='Longitude', y='Latitude') +
  theme_minimal()

## Empirical ACF
dat <- data.frame(
  y      = y,
  xcoord = locs[, 1],
  ycoord = locs[, 2],
  X      # columns from your cbind above (area_sqkm, drain_den, etc.)
)
result_resid <- sp_empirical_acf(
  y ~ area_sqkm + drain_den + MeanTempAnn + MeanPrecAnn +
    cls7 + cls9 + HydroLakes_Area_sqkm +
    meanPercentDC_Imperfectly + meanPercentDC_ModeratelyWell +
    meanPercentDC_Poor + meanPercentDC_SomewhatExcessive +
    meanPercentDC_VeryPoor + meanPercentDC_Well +
    Dam_Count + clsSum,
  data   = dat,
  xcoord = "xcoord",
  ycoord = "ycoord",
  bins   = 20
)

## Raw Data
# tune mtry
mtry_grid <- 1:ncol(X)
m_rmses <- matrix(NA, nrow = length(mtry_grid), ncol = 2)
colnames(m_rmses) <- c("mtry", "avg_rmse")

for (m in seq_along(mtry_grid)) {
  mtry <- mtry_grid[m]
  
  # parallel over folds
  fold_rmse <- foreach(
    fld = 1:nfolds,
    .combine  = "c",
    .packages = c("ranger", "yardstick")
  ) %dopar% {
    
    # random 80/20 split
    train_indices <- sample(seq_len(length(y)), size = round(0.8 * length(y)))
    y_train <- y[train_indices]
    y_val   <- y[-train_indices]
    X_train <- X[train_indices, ]
    X_val   <- X[-train_indices, ]
    
    # ranger training data
    data_train <- data.frame(y = y_train, X_train)
    
    raw_rf <- ranger::ranger(
      y ~ .,
      data       = data_train,
      mtry       = mtry,
      importance = "impurity"
    )
    
    preds <- predict(raw_rf, data = X_val)$predictions
    
    predictions <- data.frame(
      truth      = y_val,
      prediction = preds
    )
    
    yardstick::rmse(
      data     = predictions,
      truth    = truth,
      estimate = prediction
    )$.estimate
  }
  
  # average over folds for this mtry
  m_rmses[m, ] <- c(mtry, mean(fold_rmse))
}
# get best mtry and us it in random forest
plot(mtry_grid, m_rmses[,2], type="l")
best_mtry <- m_rmses[which.min(m_rmses[,2]),1]
m_rmses[which.min(m_rmses[,2]),2]
data <- cbind(y,X)
raw_rf <- ranger(y~., data=data, mtry=best_mtry, importance='impurity')
vip::vip(raw_rf)

ggplot() +
  geom_polygon(data = na_map, aes(x = long, y = lat, group = group),
               fill = "white", color="gray70") +
  geom_point(mapping=aes(x=locs[,1], y=locs[,2], color=X[,colnames(X)=="MeanPrecAnn"]), size=0.5) +
  scale_color_viridis() +
  coord_quickmap(xlim = range(locs[,1]), ylim = range(locs[,2]))


## Decorrelated Data
# tune
D <- fields::rdist(locs)
D <- fields::rdist.earth(locs, miles=FALSE)
min_range <- Matern.cor.to.range(min(D[D > 0]), nu = 1/2, cor.target = 0.05)
max_range <- Matern.cor.to.range(max(D[D > 0]), nu = 1/2, cor.target = 0.05)

tune_grid <- expand.grid(
  range  = exp(seq(log(min_range), log(max_range), length = 5)),
  nugget = seq(0.05, 0.95, length = 5),
  mtry   = seq(1, 15, length = 15)
)

K      <- 30
rmses <- foreach(
  j = 1:nrow(tune_grid),
  .combine  = rbind,
  .packages = c("fields", "ranger", "yardstick"),
  .export   = c("Exponential")  # if Exponential is your own function
) %dopar% {

  phi   <- tune_grid[j, 1]
  omega <- tune_grid[j, 2]
  mtry  <- tune_grid[j, 3]

  fold_rmse <- numeric(nfolds)

  for (fld in 1:nfolds) {
    train_indices <- sample(seq_len(length(y)), size = round(0.8 * length(y)))
    y_train <- y[train_indices]
    y_val   <- y[-train_indices]
    X_train <- X[train_indices, ]
    X_val   <- X[-train_indices, ]
    locs_train <- locs[train_indices, ]
    locs_val   <- locs[-train_indices, ]

    # decorrelate train data
    y_uncor <- matrix(NA, nrow = nrow(X_train), ncol = 1)
    x_uncor <- matrix(NA, nrow = nrow(X_train), ncol = ncol(X_train))
    y_uncor[1]  <- y_train[1]
    x_uncor[1,] <- X_train[1,]

    for (i in 2:nrow(locs_train)) {
      distances <- fields::rdist(locs_train[i, , drop = FALSE], locs_train[1:(i - 1), ])
      nbrs <- order(distances)[1:min(i - 1, K)]
      Dloc <- fields::rdist(locs_train[c(i, nbrs), ])
      R <- (1 - omega) * Exponential(Dloc, aRange = phi) + omega * diag(nrow(Dloc))
      R_inv <- solve(R[-1, -1])
      denom <- as.numeric(1 - (R[1, -1] %*% R_inv %*% R[-1, 1]))
      y_uncor[i]  <- (y_train[i]  - R[1, -1] %*% R_inv %*% y_train[nbrs])  / sqrt(denom)
      x_uncor[i,] <- (X_train[i,] - R[1, -1] %*% R_inv %*% X_train[nbrs,]) / sqrt(denom)
    }

    # decorrelate validation data
    x_uncor_val <- matrix(NA, nrow = nrow(X_val), ncol = ncol(X_val))
    for (i in 1:nrow(X_val)) {
      distances <- fields::rdist(locs_val[i, , drop = FALSE], locs_train)
      nbrs <- order(distances)[1:K]
      Dloc <- fields::rdist(rbind(locs_val[i,], locs_train[nbrs, ]))
      R <- (1 - omega) * Exponential(Dloc, aRange = phi) + omega * diag(nrow(Dloc))
      R_inv <- solve(R[-1, -1])
      denom <- as.numeric(1 - R[1, -1] %*% R_inv %*% R[1, -1])
      x_uncor_val[i,] <- (X_val[i,] - R[1, -1] %*% R_inv %*% X_train[nbrs,]) / sqrt(denom)
    }

    colnames(x_uncor_val) <- colnames(x_uncor) <- colnames(X)
    data_uncor <- data.frame(y = y_uncor, x_uncor)

    uncor_rf <- ranger::ranger(
      y ~ .,
      data       = data_uncor,
      mtry       = mtry,
      importance = "impurity",
      num.threads = 1   # good idea when you’re already parallelizing outside
    )
    preds_uncor <- predict(uncor_rf, x_uncor_val)$predictions

    # backtransform
    preds <- matrix(NA, nrow = length(preds_uncor), ncol = 1)
    preds[1] <- preds_uncor[1]
    for (i in 1:length(preds_uncor)) {
      distances <- fields::rdist(locs_val[i, , drop = FALSE], locs_train)
      nbrs <- order(distances)[1:K]
      Dloc <- fields::rdist(rbind(locs_val[i,], locs_train[nbrs, ]))
      R <- (1 - omega) * Exponential(Dloc, aRange = phi) + omega * diag(nrow(Dloc))
      R_inv <- solve(R[-1, -1])
      denom <- as.numeric(1 - R[1, -1] %*% R_inv %*% R[1, -1])
      preds[i] <- (preds_uncor[i] * sqrt(denom)) + (R[1, -1] %*% R_inv %*% y_train[nbrs])
    }

    predictions <- data.frame(truth = y_val, prediction = preds)

    fold_rmse[fld] <- yardstick::rmse(
      data     = predictions,
      truth    = truth,
      estimate = prediction
    )$.estimate
  }

  avg.rmse <- mean(fold_rmse)

  c(
    avg_rmse = avg.rmse,
    range    = phi,
    nugget   = omega,
    mtry     = mtry
  )
}

stopCluster(cl)
registerDoSEQ()

# put into a nice data.frame
rmses <- as.data.frame(rmses)

# find best
best_idx  <- order(rmses$avg_rmse)[1]
best_rmse <- rmses[best_idx, ]

phi   <- best_rmse$range
omega <- best_rmse$nugget
mtry  <- best_rmse$mtry

# decorrelate
y_uncor <- matrix(NA, nrow = nrow(X), ncol = 1)
x_uncor <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
y_uncor[1] <- y[1]
x_uncor[1,] <- X[1,]
for(i in 2:nrow(locs)){
  distances <- rdist(locs[i, ,drop=F],locs[1:(i-1),])
  nbrs <- order(distances)[1:min(i-1,K)]
  D <- rdist(locs[c(i,nbrs),])
  R <- (1-omega)*Exponential(D, aRange=phi)+omega*diag(nrow(D))
  R_inv <- solve(R[-1, -1])
  denom <- as.numeric(1-(R[1,-1] %*% R_inv %*% R[-1,1]))
  print(denom)
  y_uncor[i] <- (y[i] - R[1,-1] %*% R_inv %*% y[nbrs]) / sqrt(denom)
  x_uncor[i,] <- (X[i,] - R[1,-1] %*% R_inv %*% X[nbrs,]) / sqrt(denom)
}
colnames(x_uncor) <- colnames(X)
data_uncor <- data.frame(y=y_uncor,x_uncor)
uncor_rf <- ranger(y~., data=data_uncor, mtry=mtry, importance='impurity')
vip::vip(uncor_rf)
my_partial <- pdp::partial(
  object = uncor_rf,          # your ranger model
  pred.var = colnames(X)[11],
  train = as.data.frame(x_uncor), # MUST supply original training data
  grid.resolution = 50
)

uncor_imp <- vip::vi(uncor_rf)
cor_imp <- vip::vi(raw_rf)
bind_cols(uncor_imp, cor_imp$Importance)

# Scale importances to [0, 1] within model
uncor_scaled <- uncor_imp %>%
  mutate(Importance = Importance / max(Importance, na.rm = TRUE)) %>%
  rename(Decor = Importance)

cor_scaled <- cor_imp %>%
  mutate(Importance = Importance / max(Importance, na.rm = TRUE)) %>%
  rename(Ind = Importance)

# Merge and compute difference
compare_imp <- inner_join(uncor_scaled, cor_scaled, by = "Variable") %>%
  mutate(Diff = Decor - Ind)

# Optional: keep top variables by absolute difference
compare_top <- compare_imp %>%
  arrange(desc(abs(Diff)))

compare_top <- compare_top %>%
  mutate(Variable = recode(Variable, !!!pretty_labels))

# Plot difference
diff_plot <- ggplot(compare_top, aes(x = reorder(Variable, Diff), y = Diff,
                        fill = Diff > 0)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(limits = c(-0.77, 0.77)) +
  theme_minimal(base_size = 8) +
  labs(
    x = NULL,
    y = "Difference in VI (DecorRF – IndRF)",
  ) + 
  theme(legend.position = 'none')

df <- compare_top   # or substitute your object name

# Reshape to long format for ggplot
df_long <- df %>%
  select(Variable, Decor, Ind) %>%
  pivot_longer(cols = c("Decor", "Ind"),
               names_to = "Model",
               values_to = "Importance")

# Order variables by Decor importance (or any metric you want)
df_long$Variable <- factor(df_long$Variable,
                           levels = df %>% arrange(Diff) %>% pull(Variable))

# Plot
vi_plot <- ggplot(df_long,
       aes(x = Importance, y = Variable, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("Ind" = "grey60",
                               "Decor" = "steelblue")) +
  theme_minimal(base_size = 8) +
  labs(
    x = "Scaled Variable Importance",
    y = NULL,
    fill = "",
  ) +
  theme(
    axis.text.y = element_text(size = 9)
  )

vi_plot + diff_plot


vars_to_plot <- c("meanPercentDC_Well", "MeanPrecAnn")  

## ------------------------------------------------------------------
## 1. Helper to get a PDP for a single model & single variable
## ------------------------------------------------------------------

get_pdp <- function(model, train_df, var, model_name, grid.resolution = 50) {
  
  # pdp::partial sometimes conflicts with objects in workspace, so force clean input
  var <- var[1]  
  
  pd <- pdp::partial(
    object          = model,
    pred.var        = var,
    train           = train_df,
    grid.resolution = grid.resolution,
    plot            = FALSE
  )
  
  pd %>%
    rename(x = !!var) %>%
    mutate(
      Variable = var,
      Model    = model_name
    ) %>%
    select(Variable, Model, x, yhat)
}

## ------------------------------------------------------------------
## 2. Choose variables & build PDPs for both models
## ------------------------------------------------------------------

# indRF PDPs: ranger trained on original X
pdp_ind_list <- lapply(
  vars_to_plot,
  function(v) get_pdp(
    model      = raw_rf,
    train_df   = as.data.frame(X),
    var        = v,
    model_name = "indRF"
  )
)

# decorRF PDPs: ranger trained on decorrelated predictors x_uncor
pdp_dec_list <- lapply(
  vars_to_plot,
  function(v) get_pdp(
    model      = uncor_rf,
    train_df   = as.data.frame(x_uncor),
    var        = v,
    model_name = "decorRF"
  )
)

# Combine and center within variable-model groups
pdp_all <- bind_rows(pdp_ind_list, pdp_dec_list) %>%
  group_by(Variable, Model) %>%
  mutate(yhat_centered = yhat - mean(yhat)) %>%
  ungroup()

pdp_all <- pdp_all %>%
  mutate(
    VarLabel = ifelse(
      Variable %in% names(pretty_labels),
      pretty_labels[Variable],
      Variable   # fallback if not found in pretty_labels
    )
  )

## ------------------------------------------------------------------
## 3. 2×2 panel: rows = variables, columns = models
## ------------------------------------------------------------------
ggplot(pdp_all, aes(x = x, y = yhat_centered)) +
  geom_line(linewidth = 0.9) +
  facet_grid(VarLabel ~ Model, scales = "free_x") +
  labs(
    x = NULL,
    y = expression("Centered Partial Dependence"~hat(y))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_line(color = "grey85")
  )

