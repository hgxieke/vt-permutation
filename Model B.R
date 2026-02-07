# ============================================================
# Model B Simulation: power under Correlated Errors
# ============================================================

setwd("C:/Users/Ke/OneDrive - Washington University in St. Louis/Desktop/python/Model22/cor_error/ModelB/")

library(dplyr)
library(tibble)
library(randomForest)
library(rpart)
library(foreach)
library(doParallel)
library(ggplot2)
library(mclust)
library(MASS)

# Simulation Parameters 
nrun <- 100
nx <- 4
B <- 200
beta <- c(2, 2, 2, 2, 1, 1)
CUT0 <- 0.5
sigma <- 1
alpha <- 0.05
rho_vals <- c(0.5, 1)
n_vals <- c(300, 600)
ite_methods <- c("ite1", "ite2")

################################################################################
# ---- Data Generation and Supporting Functions ----
################################################################################

generate_data <- function(N, nx, CUT0, beta, sigma, rho, ite_method = "ite2") {
  X_df <- as.data.frame(matrix(runif(N * nx), N, nx))
  names(X_df) <- paste0("X", 1:nx)
  Z1 <- as.integer(X_df$X1 <= CUT0)
  Z2 <- as.integer(X_df$X2 <= CUT0)
  trt <- rbinom(N, 1, 0.5)
  
  # Correlated errors
  # eps0 and eps1 have mean 0, variance sigma^2, and correlation rho
  Sigma_eps <- matrix(c(1, rho, rho, 1), ncol = 2) * sigma^2
  eps_pair <- MASS::mvrnorm(n = N, mu = c(0, 0), Sigma = Sigma_eps)
  eps0_full <- eps_pair[, 1]
  eps1_full <- eps_pair[, 2]
  
  mu <- beta[1] +
    beta[2] * trt +
    beta[3] * Z1 +
    beta[4] * Z2 +
    beta[5] * trt *Z1 +
    beta[6] * (trt * Z1 * Z2)
  
  y <- numeric(N)
  y[trt == 0] <- mu[trt == 0] + eps0_full[trt == 0]
  y[trt == 1] <- mu[trt == 1] + eps1_full[trt == 1]
  
  df <- bind_cols(
    tibble(id = 1:N, trt = trt, Z1 = Z1, Z2 = Z2, y = y),
    X_df
  )
  
  # Random forest fits
  rf0 <- randomForest(y ~ X1 + X2 + X3 + X4, data = filter(df, trt == 0))
  rf1 <- randomForest(y ~ X1 + X2 + X3 + X4, data = filter(df, trt == 1))
  
  df <- df %>%
    mutate(pred0 = predict(rf0, newdata = df),
           pred1 = predict(rf1, newdata = df))
  
  if (ite_method == "ite1") {
    df <- df %>%
      mutate(ITE = ifelse(trt == 1, y - pred0, pred1 - y))
  } else if (ite_method == "ite2") {
    df <- df %>%
      mutate(ITE = pred1 - pred0)
  }
  return(df)
}


permute_data <- function(data, ite_method = "ite2") {
  zbar <- mean(data$y[data$trt == 1]) - mean(data$y[data$trt == 0])
  data_p <- data
  data_p$y[data_p$trt == 1] <- data_p$y[data_p$trt == 1] - zbar
  data_p$trt <- sample(data_p$trt, size = nrow(data_p))
  rf0_p <- randomForest(y ~ X1 + X2 + X3 + X4, data = filter(data_p, trt == 0))
  rf1_p <- randomForest(y ~ X1 + X2 + X3 + X4, data = filter(data_p, trt == 1))
  data_p$pred0 <- predict(rf0_p, newdata = data_p)
  data_p$pred1 <- predict(rf1_p, newdata = data_p)
  
  if (ite_method == "ite1") {
    data_p <- data_p %>% mutate(ITE_p = ifelse(trt == 1, y - pred0, pred1 - y))
  } else if (ite_method == "ite2") {
    data_p <- data_p %>% mutate(ITE_p = pred1 - pred0)
  }
  return(data_p)
}


get_mnse <- function(permuted_data) {
  full_tree <- rpart(ITE_p ~ X1 + X2 + X3 + X4,
                     data = permuted_data,
                     control = rpart.control(cp = 0, xval = 10))
  cp_tab <- full_tree$cptable
  opt_idx <- which.min(cp_tab[, "xerror"])
  xerror_opt <- cp_tab[opt_idx, "xerror"]
  xstd_opt <- cp_tab[opt_idx, "xstd"]
  xerror_chosen <- cp_tab[1, "xerror"]
  if (xstd_opt == 0) return(NA)
  (xerror_chosen - xerror_opt) / xstd_opt
}


prune_cv <- function(formula, data, gamma) {
  full <- rpart(formula, data = data, control = rpart.control(cp = 0, xval = 10, minsplit = 30))
  cp_tab <- full$cptable
  opt <- which.min(cp_tab[, "xerror"])
  err <- cp_tab[, "xerror"]; std <- cp_tab[, "xstd"]
  cutoff <- min(err) + gamma * std[opt]
  pos <- min(which(err <= cutoff))
  cp1 <- cp_tab[pos, "CP"]
  prune(full, cp = cp1)
}

################################################################################
# ---- Simulation Loop (rho, N, ite_method) ----
################################################################################

results_summary <- data.frame()

for (rho in rho_vals) {
  for (N in n_vals) {
    for (ite_method_to_use in ite_methods) {
      
      cat("Running simulation for rho =", rho, ", N =", N, ", method =", ite_method_to_use, "\n")
      
      # Prepare storage
      mnse_matrix <- matrix(NA, nrow = nrun, ncol = B)
      ari_vector <- numeric(nrun)
      gamma_vec <- numeric(nrun)
      tree_sizes <- numeric(nrun)
      
      n_cores <- parallel::detectCores() - 1
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      
      start_time <- Sys.time()
      for (sim in 1:nrun) {
        set.seed(sim)
        original_data <- generate_data(N, nx, CUT0, beta, sigma, rho, ite_method = ite_method_to_use)
        mnse_distribution <- foreach(i = 1:B, .combine = c, .packages = c("dplyr","randomForest","rpart","tibble")) %dopar% {
          permuted_df <- permute_data(original_data, ite_method = ite_method_to_use)
          get_mnse(permuted_df)
        }
        mnse_matrix[sim, ] <- mnse_distribution
        gamma_threshold <- quantile(na.omit(mnse_distribution), probs = 1 - alpha, na.rm = TRUE)
        gamma_vec[sim] <- gamma_threshold
        pruned_tree <- prune_cv(ITE ~ X1 + X2 + X3 + X4, data = original_data, gamma = gamma_threshold)
        predicted_groups <- pruned_tree$where

        true_groups <- ifelse(original_data$Z1 == 0, 1,
                              ifelse(original_data$Z2 == 0, 2, 3))
        tree_sizes[sim] <- length(unique(predicted_groups)) - 1
        ari_vector[sim] <- if (length(unique(predicted_groups)) > 1) adjustedRandIndex(predicted_groups, true_groups) else 0
      }
      stopCluster(cl)
      
      results_df <- data.frame(
        Adjusted_Rand_Index = ari_vector,
        Gamma = gamma_vec,
        Tree_Size = tree_sizes
      )
      
      # Save each combination’s raw results
      fname <- paste0("rho", rho, "_N", N, "_", ite_method_to_use, ".csv")
      write.csv(results_df, fname, row.names = FALSE)
      
      # Compute means ± SDs
      mean_gamma <- mean(results_df$Gamma, na.rm = TRUE)
      sd_gamma <- sd(results_df$Gamma, na.rm = TRUE)
      mean_tree <- mean(results_df$Tree_Size, na.rm = TRUE)
      mean_ari <- mean(results_df$Adjusted_Rand_Index, na.rm = TRUE)
      sd_ari <- sd(results_df$Adjusted_Rand_Index, na.rm = TRUE)
      
      # Add summary to overall table
      results_summary <- rbind(
        results_summary,
        data.frame(
          Method = ite_method_to_use,
          rho = rho,
          N = N,
          Gamma = sprintf("%.3f ± %.3f", mean_gamma, sd_gamma),
          Tree_Size = round(mean_tree, 2),
          ARI = sprintf("%.3f ± %.3f", mean_ari, sd_ari)
        )
      )
    }
  }
}

# Save final summary table
write.csv(results_summary, "ModelB.csv", row.names = FALSE)
print(results_summary)
