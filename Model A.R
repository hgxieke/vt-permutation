# ============================================================
# Model A Simulation: Type I Error under Correlated Errors
# ============================================================

setwd("C:/Users/Ke/OneDrive - Washington University in St. Louis/Desktop/python/Model22/cor_error/ModelA/")

library(dplyr)
library(tibble)
library(randomForest)
library(rpart)
library(foreach)
library(doParallel)
library(MASS)

# Simulation Parameters 
nrun <- 500        # simulations runs
nx <- 4           # number of covariates
B <- 200          # permutations per simulation run
CUT0 <- 0.5
beta <- c(2, 2, 2, 2)  # null model (no interaction)
sigma <- 1
alpha <- 0.05
rho_vals <- c(0, 0.5, 1)
n_vals <- c(300, 600)
ite_methods <- c("ite1", "ite2")  # ite1 for o-ITE method, ite2 for m-ITE method

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
  
  mu <- beta[1] + beta[2] * trt + beta[3] * Z1 + beta[4] * Z2
  
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

get_mnse <- function(permuted_data) {
  formula_rpart <- as.formula(paste("ITE_p ~", paste0("X", 1:nx, collapse = "+")))
  full_tree <- rpart::rpart(formula = formula_rpart,
                            data = permuted_data,
                            control = rpart.control(cp = 0, xval = 10))
  cp_tab <- full_tree$cptable
  # Get the minimum cross-validation error and its standard deviation
  opt_idx <- which.min(cp_tab[, "xerror"])
  xerror_opt <- cp_tab[opt_idx, "xerror"]
  xstd_opt <- cp_tab[opt_idx, "xstd"]
  # Get the cross-validation error for the root-only tree (nsplit=0)
  xerror_chosen <- cp_tab[1, "xerror"]
  
  if (xstd_opt == 0) return(NA)
  gamma <- (xerror_chosen - xerror_opt) / xstd_opt
  return(gamma)
}

prune_cv <- function(formula, data, minsplit = 30, rule = c("0-SE", "gamma-SE"), gamma) {
  
  rule <- match.arg(rule)
  # Grow a full rpart tree
  full <- rpart(formula, data = data,
                control = rpart.control(cp = 0, xval = 10, minsplit = minsplit))
  # Complexity parameter table
  cp_tab <- full$cptable
  # Row with the minimum cross-validation error
  opt <- which.min(cp_tab[, "xerror"])
  
  # If the "0-SE" rule is chosen, return the tree with the lowest error
  if (rule == "0-SE") {
    cp0 <- cp_tab[opt, "CP"]
    return(prune(full, cp = cp0))
  }
  # If the "gamma-SE" rule is chosen, proceed with the custom pruning
  err <- cp_tab[, "xerror"]
  std <- cp_tab[, "xstd"]
  
  # Calculate the cutoff using the provided gamma value
  cutoff <- min(err) + gamma * std[opt]
  # Find the position of the simplest tree (minimum cp) whose error is below the cutoff
  pos <- min(which(err <= cutoff))
  # Get the corresponding complexity parameter
  cp1 <- cp_tab[pos, "CP"]
  # Prune the full tree and return the result
  prune(full, cp = cp1)
}


# Simulation loop
summary_table <- list()

for (N in n_vals) {
  for (rho in rho_vals) {
    for (ite_method_to_use in ite_methods) {
      
      cat("N=", N, ", rho=", rho, ", method=", ite_method_to_use, "\n")
      
      gamma_thresholds <- numeric(nrun)
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
                                        get_mnse(permuted_df)
                                      }
        
        gamma_distribution <- na.omit(gamma_distribution)
        gamma_threshold <- if (length(gamma_distribution) > 0) {
          quantile(gamma_distribution, probs = 1 - alpha, na.rm = TRUE)
        } else { NA }
        
        gamma_thresholds[sim] <- gamma_threshold
        
        pruned_tree <- prune_cv(
          formula = ITE ~ X1 + X2 + X3 + X4, 
          data = original_data, 
          rule = "gamma-SE", 
          gamma = gamma_threshold
        )
        
        did_split[sim] <- (nrow(pruned_tree$frame) > 1)
        
        cat(paste0("Sim ", sim, " done. Split=", did_split[sim], "\n"))
      }
      
      stopCluster(cl)
      
      false_positive_rate <- mean(did_split, na.rm = TRUE)
      mean_gamma <- mean(gamma_thresholds, na.rm = TRUE)
      sd_gamma <- sd(gamma_thresholds, na.rm = TRUE)
      
      result_df <- data.frame(
        sim = 1:nrun,
        Gamma = gamma_thresholds,
        DidSplit = as.integer(did_split)
      )
      
      # Save results per combination
      fname <- paste0("rho", rho, "_n", N, "_", ite_method_to_use, ".csv")
      write.csv(result_df, fname, row.names = FALSE)
      
      summary_table[[paste(rho, N, ite_method_to_use, sep = "_")]] <- 
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
summary <- bind_rows(summary_table)
write.csv(summary, "ModelA_beta2.csv", row.names = FALSE)
summary
