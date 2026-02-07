# ============================================================
# Model D Simulation: Type I Error under Correlated Errors
# ============================================================

setwd("C:/Users/Ke/OneDrive - Washington University in St. Louis/Desktop/python/Model22/cor_error/ModelD/")

library(dplyr)
library(tibble)
library(randomForest)
library(rpart)
library(foreach)
library(doParallel)
library(MASS)

# Simulation Parameters 
nrun <- 500        # simulations per setting
nx <- 4           # number of covariates
B <- 200          # permutations per simulation
CUT0 <- 0.5
beta <- c(2, 2)   # null model (no interaction)
sigma <- 1
alpha <- 0.05
rho_vals <- c(0, 0.5, 1)
n_vals <- c(300, 600)
ite_methods <- c("ite1", "ite2")

################################################################################
# Functions 
################################################################################
generate_data_null <- function(N, nx, CUT0, beta, sigma, rho, ite_method = "ite2") {
  X_df <- as.data.frame(matrix(runif(N * nx), N, nx))
  names(X_df) <- paste0("X", 1:nx)
  Z1 <- as.integer(X_df$X1 <= CUT0)
  Z2 <- as.integer(X_df$X2 <= CUT0)
  trt <- rbinom(N, 1, 0.5)
  
  cov_matrix <- matrix(c(1, rho, rho, 1), ncol = 2) * sigma^2
  eps_pair <- MASS::mvrnorm(n = N, mu = c(0, 0), Sigma = cov_matrix)
  eps0_full <- eps_pair[, 1]
  eps1_full <- eps_pair[, 2]
  
  # Compute mu_i under the NULL hypothesis (no interaction term)
  mu <- beta[1] +  beta[2] * exp((X_df$X1 - 0.5)^2 + (X_df$X2 - 0.5)^2)  
  
  y <- numeric(N)
  y[trt==0] <- mu[trt==0] + eps0_full[trt==0]
  y[trt==1] <- mu[trt==1] + eps1_full[trt==1]
  
  df <- bind_cols(tibble(id = 1:N, trt = trt, Z1 = Z1, Z2 = Z2, y = y), X_df)
  
  rf0 <- randomForest(y ~ X1 + X2 + X3 + X4, data = filter(df, trt == 0))
  rf1 <- randomForest(y ~ X1 + X2 + X3 + X4, data = filter(df, trt == 1))
  
  df <- df %>%
    mutate(pred0 = predict(rf0, newdata = df),
           pred1 = predict(rf1, newdata = df))
  
  if (ite_method == "ite1") {
    df <- df %>% mutate(ITE = ifelse(trt == 1, y - pred0, pred1 - y))
  } else if (ite_method == "ite2") {
    df <- df %>% mutate(ITE = pred1 - pred0)
  } else stop("Invalid ite_method.")
  
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

calculate_mnse_multiplier <- function(permuted_data) {
  formula_rpart <- as.formula(paste("ITE_p ~", paste0("X", 1:nx, collapse = "+")))
  full_tree <- rpart::rpart(formula = formula_rpart,
                            data = permuted_data,
                            control = rpart.control(cp = 0, xval = 10))
  cp_tab <- full_tree$cptable
  opt_idx <- which.min(cp_tab[, "xerror"])
  xerror_opt <- cp_tab[opt_idx, "xerror"]
  xstd_opt <- cp_tab[opt_idx, "xstd"]
  xerror_chosen <- cp_tab[1, "xerror"]
  
  if (xstd_opt == 0) return(NA)
  gamma <- (xerror_chosen - xerror_opt) / xstd_opt
  return(gamma)
}

prune_cv <- function(formula, data, gamma) {
  full <- rpart(formula, data = data,
                control = rpart.control(cp = 0, xval = 10, minsplit = 30))
  cp_tab <- full$cptable
  opt <- which.min(cp_tab[, "xerror"])
  err <- cp_tab[, "xerror"]
  std <- cp_tab[, "xstd"]
  cutoff <- min(err) + gamma * std[opt]
  pos <- min(which(err <= cutoff))
  cp1 <- cp_tab[pos, "CP"]
  prune(full, cp = cp1)
}

# Simulation loop
summary_results <- list()

for (N in n_vals) {
  for (rho in rho_vals) {
    for (ite_method_to_use in ite_methods) {
      
      cat("Running combination: N=", N, ", rho=", rho, ", method=", ite_method_to_use, "\n")
      
      gamma_thresholds_per_run <- numeric(nrun)
      did_split <- logical(nrun)
      
      n_cores <- parallel::detectCores() - 1
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      
      for (sim in 1:nrun) {
        set.seed(sim)
        original_data <- generate_data_null(N, nx, CUT0, beta, sigma, rho, ite_method_to_use)
        
        gamma_distribution <- foreach(i = 1:B, .combine = c,
                                      .packages = c("dplyr", "randomForest", "rpart", "tibble")) %dopar% {
                                        permuted_df <- permute_data(original_data, ite_method_to_use)
                                        calculate_mnse_multiplier(permuted_df)
                                      }
        
        gamma_distribution <- na.omit(gamma_distribution)
        gamma_threshold <- if (length(gamma_distribution) > 0) {
          quantile(gamma_distribution, probs = 1 - alpha, na.rm = TRUE)
        } else { NA }
        
        gamma_thresholds_per_run[sim] <- gamma_threshold
        
        pruned_tree <- prune_cv(formula = ITE ~ X1 + X2 + X3 + X4,
                                data = original_data, gamma = gamma_threshold)
        did_split[sim] <- (nrow(pruned_tree$frame) > 1)
        
        cat(paste0("Sim ", sim, " done. Split=", did_split[sim], "\n"))
      }
      
      stopCluster(cl)
      
      # Summary stats
      false_positive_rate <- mean(did_split, na.rm = TRUE)
      mean_gamma <- mean(gamma_thresholds_per_run, na.rm = TRUE)
      sd_gamma <- sd(gamma_thresholds_per_run, na.rm = TRUE)
      
      result_df <- data.frame(
        sim = 1:nrun,
        Gamma = gamma_thresholds_per_run,
        DidSplit = as.integer(did_split)
      )
      
      # Save per combination
      fname <- paste0("rho", rho, "_n", N, "_", ite_method_to_use, ".csv")
      write.csv(result_df, fname, row.names = FALSE)
      
      summary_results[[paste(rho, N, ite_method_to_use, sep = "_")]] <- 
        data.frame(Method = ite_method_to_use,
                   Rho = rho,
                   N = N,
                   MeanGamma = round(mean_gamma, 3),
                   SDGamma = round(sd_gamma, 3),
                   FPR = round(false_positive_rate, 3))
    }
  }
}

# summary table 
summary <- bind_rows(summary_results)
write.csv(summary, "ModelD_beta2.csv", row.names = FALSE)
summary
