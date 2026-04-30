#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(glmnet)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/benchmark_changepoint_fused_vs_glmnet.R")
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  root <- normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(root, "R", "l1_fusion.R"))
source(file.path(root, "R", "l1_fusion_new_utils.R"))
source(file.path(root, "R", "l1_fusion_operator_new.R"))

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

predict_from_beta <- function(beta, X, groups) {
  yhat <- numeric(nrow(X))
  group_levels <- sort(unique(groups))
  for (g in group_levels) {
    idx <- which(groups == g)
    yhat[idx] <- as.numeric(X[idx, , drop = FALSE] %*% beta[, g, drop = FALSE])
  }
  yhat
}

make_block_design <- function(X, groups, k) {
  n <- nrow(X)
  p <- ncol(X)
  row_idx <- rep(seq_len(n), each = p)
  col_idx <- as.vector(t((groups - 1L) * p + matrix(seq_len(p), nrow = n, ncol = p, byrow = TRUE)))
  Matrix::sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = as.vector(t(X)),
    dims = c(n, p * k)
  )
}

generate_piecewise_data <- function(seed = 20260430, k = 60, p = 20,
                                    n_group_fit = 25, n_group_val = 10,
                                    n_group_test = 80, sigma = 0.7) {
  set.seed(seed)
  breakpoints <- c(20L, 40L)
  segment <- cut(seq_len(k), breaks = c(0L, breakpoints, k), labels = FALSE)

  beta <- matrix(0, nrow = p, ncol = k)
  beta[1, ] <- c(1.6, -1.1, 1.2)[segment]
  beta[2, ] <- c(0.0, 1.3, -1.0)[segment]
  beta[3, ] <- c(-1.2, -1.2, 0.8)[segment]
  beta[4, ] <- 0.8
  beta[5, ] <- c(0.5, 0.0, 0.0)[segment]

  make_split <- function(n_per_group) {
    groups <- rep(seq_len(k), each = n_per_group)
    X <- matrix(rnorm(length(groups) * p), ncol = p)
    mean_y <- rowSums(X * t(beta[, groups, drop = FALSE]))
    y <- mean_y + rnorm(length(groups), sd = sigma)
    list(X = X, y = y, groups = groups)
  }

  fit <- make_split(n_group_fit)
  val <- make_split(n_group_val)
  test <- make_split(n_group_test)

  x_center <- colMeans(fit$X)
  x_scale <- apply(fit$X, 2, sd)
  x_scale[x_scale == 0] <- 1
  standardize <- function(X) sweep(sweep(X, 2, x_center, "-"), 2, x_scale, "/")

  fit$X <- standardize(fit$X)
  val$X <- standardize(val$X)
  test$X <- standardize(test$X)

  G <- matrix(0, k, k)
  for (i in seq_len(k - 1L)) {
    G[i, i + 1L] <- 1
    G[i + 1L, i] <- 1
  }

  list(fit = fit, val = val, test = test, beta = beta, breakpoints = breakpoints, G = G)
}

fit_pooled_glmnet <- function(d) {
  cv <- glmnet::cv.glmnet(d$fit$X, d$fit$y, alpha = 1, standardize = FALSE, nfolds = 5)
  yhat_val <- as.numeric(predict(cv, d$val$X, s = "lambda.min"))
  yhat_test <- as.numeric(predict(cv, d$test$X, s = "lambda.min"))
  coef_vec <- as.numeric(coef(cv, s = "lambda.min"))[-1]
  beta <- matrix(rep(coef_vec, ncol(d$beta)), nrow = nrow(d$beta))
  list(beta = beta, val_rmse = rmse(d$val$y, yhat_val), test_rmse = rmse(d$test$y, yhat_test))
}

fit_group_glmnet <- function(d) {
  k <- ncol(d$beta)
  Z_fit <- make_block_design(d$fit$X, d$fit$groups, k)
  Z_val <- make_block_design(d$val$X, d$val$groups, k)
  Z_test <- make_block_design(d$test$X, d$test$groups, k)
  cv <- glmnet::cv.glmnet(Z_fit, d$fit$y, alpha = 1, standardize = FALSE, nfolds = 5)
  yhat_val <- as.numeric(predict(cv, Z_val, s = "lambda.min"))
  yhat_test <- as.numeric(predict(cv, Z_test, s = "lambda.min"))
  coef_vec <- as.numeric(coef(cv, s = "lambda.min"))[-1]
  beta <- matrix(coef_vec, nrow = nrow(d$beta), ncol = k)
  list(beta = beta, val_rmse = rmse(d$val$y, yhat_val), test_rmse = rmse(d$test$y, yhat_test))
}

fit_fused_l1 <- function(d) {
  grid <- expand.grid(
    lambda = c(0, 1e-4, 1e-3, 1e-2),
    gamma = c(0.01, 0.05, 0.1, 0.25, 0.5),
    KEEP.OUT.ATTRS = FALSE
  )
  fits <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    beta <- suppressWarnings(fusedLassoProximalNewOperator(
      d$fit$X, d$fit$y, d$fit$groups,
      lambda = grid$lambda[i],
      gamma = grid$gamma[i],
      G = d$G,
      mu = 1e-3,
      tol = 1e-4,
      num.it = 2000,
      intercept = FALSE,
      scaling = TRUE,
      edge.block = 128L
    ))
    yhat_val <- predict_from_beta(beta, d$val$X, d$val$groups)
    fits[[i]] <- list(beta = beta, val_rmse = rmse(d$val$y, yhat_val))
  }
  best <- which.min(vapply(fits, `[[`, numeric(1), "val_rmse"))
  beta <- fits[[best]]$beta
  yhat_test <- predict_from_beta(beta, d$test$X, d$test$groups)
  list(
    beta = beta,
    val_rmse = fits[[best]]$val_rmse,
    test_rmse = rmse(d$test$y, yhat_test),
    lambda = grid$lambda[best],
    gamma = grid$gamma[best]
  )
}

top_jumps <- function(beta, n_cp) {
  score <- sqrt(colSums((beta[, -1, drop = FALSE] - beta[, -ncol(beta), drop = FALSE])^2))
  sort(order(score, decreasing = TRUE)[seq_len(n_cp)])
}

jump_count <- function(beta, threshold = 0.25) {
  score <- sqrt(colSums((beta[, -1, drop = FALSE] - beta[, -ncol(beta), drop = FALSE])^2))
  sum(score > threshold)
}

cp_mean_distance <- function(estimated, truth) {
  mean(vapply(truth, function(cp) min(abs(estimated - cp)), numeric(1)))
}

summarize_fit <- function(name, fit, d) {
  estimated_cp <- top_jumps(fit$beta, length(d$breakpoints))
  data.frame(
    method = name,
    test_rmse = fit$test_rmse,
    beta_rel_error = sqrt(sum((fit$beta - d$beta)^2) / sum(d$beta^2)),
    cp_mean_distance = cp_mean_distance(estimated_cp, d$breakpoints),
    jump_count_025 = jump_count(fit$beta, threshold = 0.25),
    estimated_changepoints = paste(estimated_cp, collapse = ","),
    lambda = fit$lambda %||% NA_real_,
    gamma = fit$gamma %||% NA_real_
  )
}

plot_paths <- function(d, fits, out_png) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    png(out_png, width = 1800, height = 900, res = 180)
    matplot(t(d$beta[1:3, ]), type = "l", lty = 1, lwd = 2,
            main = "True coefficient paths", xlab = "Group", ylab = "Coefficient")
    abline(v = d$breakpoints, lty = 2)
    dev.off()
    return(invisible(NULL))
  }
  df <- do.call(rbind, lapply(names(fits), function(method) {
    beta <- fits[[method]]$beta
    do.call(rbind, lapply(seq_len(3), function(feature) {
      data.frame(method = method, group = seq_len(ncol(beta)),
                 feature = paste0("x", feature), beta = beta[feature, ])
    }))
  }))
  true_df <- do.call(rbind, lapply(seq_len(3), function(feature) {
    data.frame(method = "truth", group = seq_len(ncol(d$beta)),
               feature = paste0("x", feature), beta = d$beta[feature, ])
  }))
  df <- rbind(true_df, df)
  df$method <- factor(df$method, levels = c("truth", "pooled glmnet", "group glmnet", "fused L1"))

  p <- ggplot2::ggplot(df, ggplot2::aes(group, beta, color = method, linewidth = method)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = d$breakpoints, linetype = "dashed", color = "grey45") +
    ggplot2::facet_wrap(~ feature, ncol = 1, scales = "free_y") +
    ggplot2::scale_linewidth_manual(values = c("truth" = 1.2, "pooled glmnet" = 0.7,
                                               "group glmnet" = 0.7, "fused L1" = 0.9)) +
    ggplot2::labs(x = "Ordered group", y = "Coefficient",
                  color = NULL, linewidth = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(out_png, p, width = 7.0, height = 4.8, dpi = 300)
}

main <- function() {
  d <- generate_piecewise_data()
  fits <- list(
    "pooled glmnet" = fit_pooled_glmnet(d),
    "group glmnet" = fit_group_glmnet(d),
    "fused L1" = fit_fused_l1(d)
  )
  summary <- do.call(rbind, Map(summarize_fit, names(fits), fits, MoreArgs = list(d = d)))

  out_dir <- file.path(root, "build", "changepoint")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(summary, file.path(out_dir, "fused_vs_glmnet_summary.csv"), row.names = FALSE)

  fig_path <- file.path(root, "publication", "paper1", "figures", "changepoint_fused_vs_glmnet.png")
  plot_paths(d, fits, fig_path)
  print(summary)
}

main()
