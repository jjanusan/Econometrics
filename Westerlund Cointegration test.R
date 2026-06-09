new_panel <- new_panel %>%
  rename(id = geo, time = year)

# Newey-West variance estimator function
NW_var <- function(v, L = 1) {
  T <- length(v)
  v_bar <- mean(v)
  gamma0 <- sum((v - v_bar)^2) / T
  sum_cov <- 0
  if (L > 0) {
    for (l in 1:L) {
      # compute lag-l autocovariance
      gamma_l <- sum((v[(l+1):T] - v_bar) * (v[1:(T-l)] - v_bar)) / T
      weight <- 1 - l/(L+1)
      sum_cov <- sum_cov + 2 * weight * gamma_l
    }
  }
  gamma0 + sum_cov
}

# Function to perform the cointegration test for one cross-sectional unit
compute_unit_test <- function(data_unit, NW_lag = 1) {
  # Ensure data is sorted by time
  data_unit <- data_unit[order(data_unit$time), ]
  T_unit <- nrow(data_unit)
  if (T_unit < 2) stop("Not enough time observations for unit ", unique(data_unit$id))
  
  # Create differenced y and lagged values for y and x
  data_unit <- data_unit %>%
    mutate(dy = c(NA, diff(y)),
           y_lag = lag(y, 1),
           x_lag = lag(x, 1))
  
  # Remove the first observation 
  data_reg <- na.omit(data_unit)
  if(nrow(data_reg) < 5) stop("Too few observations after differencing for unit ", unique(data_unit$id))
  
  # Estimate the error correction model:
  #   dy = δ + α_i * y_lag + λ_i * x_lag + error
  model <- lm(dy ~ y_lag + x_lag, data = data_reg)
  
  # Extract the coefficient and its standard error for y_lag (the adjustment parameter α_i)
  coef_alpha <- coef(model)["y_lag"]
  se_alpha   <- summary(model)$coefficients["y_lag", "Std. Error"]
  
  # Obtain residuals from the model and the differenced y vector
  res <- residuals(model)
  dy_vec <- data_reg$dy
  
  # Compute Newey-West long-run variance estimates
  omega_u <- NW_var(res, L = NW_lag)
  omega_y <- NW_var(dy_vec, L = NW_lag)
  
  # Compute adjustment parameter for standardization: α_i(1)
  alpha1 <- omega_u / omega_y
  
  return(list(alpha = coef_alpha, se = se_alpha, alpha1 = alpha1, T = nrow(data_reg)))
}

# Main function: loop over all cross-sectional units and compute group-mean test statistics
cointegration_test <- function(df, NW_lag = 1) {
  # df must contain: id, time, y, x
  # Split data by cross-sectional unit
  units <- split(df, df$id)
  
  # Run the unit-level test and store the results for each unit
  results <- lapply(units, compute_unit_test, NW_lag = NW_lag)
  
  # Convert results to a data frame
  res_df <- do.call(rbind, lapply(names(results), function(id) {
    r <- results[[id]]
    data.frame(id = id, alpha = r$alpha, se = r$se, alpha1 = r$alpha1, T = r$T)
  }))
  
  # Number of cross-sectional units
  N <- nrow(res_df)
  
  # Compute group-mean test statistics:
  #   G_tau = (1/N)*sum(alpha_i / SE(alpha_i))
  G_tau <- mean(res_df$alpha / res_df$se)
  
  #   G_alpha = (1/N)*sum(T_i * alpha_i / alpha1_i)
  G_alpha <- mean(res_df$T * res_df$alpha / res_df$alpha1)
  
  return(list(results_by_unit = res_df, G_tau = G_tau, G_alpha = G_alpha))
}

# List of independent variable names to test against gdp_ldiff
indep_vars <- c("RnDexp_ldiff", "REexp_ldiff", "invst_ldiff", "lbr_ldiff")
results_list <- list()

# Loop through each independent variable, prepare the data, and run the cointegration test
for (var in indep_vars) {
  # Create a temporary dataset with the required column names:
  # y: dependent variable (gdp_ldiff)
  # x: independent variable (one of the variables in indep_vars)
  pdata_temp <- pdata %>%
    select(id, time, y = gdp_ldiff, x = all_of(var))
  
  cat("Testing cointegration between gdp_ldiff and", var, "\n")
  
  # Run the cointegration test for the current independent variable
  results <- cointegration_test(pdata_temp, NW_lag = 1)
  results_list[[var]] <- results
  
  print(results$results_by_unit)
  cat("Group-mean test statistic G_tau =", results$G_tau, "\n")
  cat("Group-mean test statistic G_alpha =", results$G_alpha, "\n\n")
}


