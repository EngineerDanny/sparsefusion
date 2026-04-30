#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(glmnet)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "scripts/benchmark_spatial_hotspot_fused_vs_glmnet.R")
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

grid_edges <- function(side) {
  edge_list <- list()
  add_edge <- function(a, b) {
    edge_list[[length(edge_list) + 1L]] <<- c(a, b)
  }
  id <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) {
    for (c in seq_len(side)) {
      if (r < side) add_edge(id(r, c), id(r + 1L, c))
      if (c < side) add_edge(id(r, c), id(r, c + 1L))
    }
  }
  do.call(rbind, edge_list)
}

generate_hotspot_data <- function(seed = 20260430, side = 8, p = 16,
                                  n_group_fit = 25, n_group_val = 10,
                                  n_group_test = 70, sigma = 0.75) {
  set.seed(seed)
  k <- side * side
  coords <- expand.grid(row = seq_len(side), col = seq_len(side))
  hotspot1 <- coords$row %in% 2:4 & coords$col %in% 2:4
  hotspot2 <- coords$row %in% 6:7 & coords$col %in% 6:8

  region <- rep("background", k)
  region[hotspot1] <- "hotspot1"
  region[hotspot2] <- "hotspot2"

  beta <- matrix(0, nrow = p, ncol = k)
  beta[1, ] <- c(background = 0.4, hotspot1 = 1.8, hotspot2 = -1.1)[region]
  beta[2, ] <- c(background = -0.5, hotspot1 = 0.9, hotspot2 = 1.4)[region]
  beta[3, ] <- c(background = 0.0, hotspot1 = -1.2, hotspot2 = 0.6)[region]
  beta[4, ] <- c(background = 0.7, hotspot1 = 0.7, hotspot2 = -0.6)[region]
  beta[5, ] <- c(background = 0.0, hotspot1 = 0.5, hotspot2 = 0.0)[region]

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

  edges <- grid_edges(side)
  G <- matrix(0, k, k)
  for (i in seq_len(nrow(edges))) {
    G[edges[i, 1], edges[i, 2]] <- 1
    G[edges[i, 2], edges[i, 1]] <- 1
  }

  list(fit = fit, val = val, test = test, beta = beta, G = G,
       edges = edges, coords = coords, region = region, side = side)
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
    gamma = c(0.05, 0.1, 0.25, 0.5, 1.0),
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
      num.it = 1800,
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

edge_scores <- function(beta, edges) {
  sqrt(colSums((beta[, edges[, 1], drop = FALSE] - beta[, edges[, 2], drop = FALSE])^2))
}

boundary_metrics <- function(beta, d, threshold = 0.45) {
  truth <- d$region[d$edges[, 1]] != d$region[d$edges[, 2]]
  detected <- edge_scores(beta, d$edges) > threshold
  tp <- sum(truth & detected)
  fp <- sum(!truth & detected)
  fn <- sum(truth & !detected)
  precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  f1 <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  c(boundary_f1 = f1, boundary_precision = precision, boundary_recall = recall,
    large_boundaries = sum(detected))
}

summarize_fit <- function(name, fit, d) {
  bm <- boundary_metrics(fit$beta, d)
  data.frame(
    method = name,
    test_rmse = fit$test_rmse,
    beta_rel_error = sqrt(sum((fit$beta - d$beta)^2) / sum(d$beta^2)),
    boundary_f1 = bm[["boundary_f1"]],
    large_boundaries = bm[["large_boundaries"]],
    lambda = fit$lambda %||% NA_real_,
    gamma = fit$gamma %||% NA_real_
  )
}

plot_hotspots <- function(d, fits, out_png) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    png(out_png, width = 1200, height = 800, res = 160)
    image(matrix(d$beta[1, ], nrow = d$side, byrow = TRUE), main = "True x1 coefficient")
    dev.off()
    return(invisible(NULL))
  }
  df <- do.call(rbind, lapply(names(fits), function(method) {
    beta <- fits[[method]]$beta
    data.frame(method = method, d$coords, beta = beta[1, ])
  }))
  truth_df <- data.frame(method = "truth", d$coords, beta = d$beta[1, ])
  df <- rbind(truth_df, df)
  df$method <- factor(df$method, levels = c("truth", "pooled glmnet", "group glmnet", "fused L1"))

  p <- ggplot2::ggplot(df, ggplot2::aes(col, row, fill = beta)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_y_reverse(breaks = seq_len(d$side)) +
    ggplot2::scale_x_continuous(breaks = seq_len(d$side)) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
    ggplot2::coord_equal() +
    ggplot2::facet_wrap(~ method, nrow = 1) +
    ggplot2::labs(x = "Grid column", y = "Grid row", fill = "x1 coef.") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom")
  ggplot2::ggsave(out_png, p, width = 7.0, height = 2.7, dpi = 300)
}

main <- function() {
  d <- generate_hotspot_data()
  fits <- list(
    "pooled glmnet" = fit_pooled_glmnet(d),
    "group glmnet" = fit_group_glmnet(d),
    "fused L1" = fit_fused_l1(d)
  )
  summary <- do.call(rbind, Map(summarize_fit, names(fits), fits, MoreArgs = list(d = d)))

  out_dir <- file.path(root, "build", "spatial_hotspot")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(summary, file.path(out_dir, "fused_vs_glmnet_summary.csv"), row.names = FALSE)

  fig_path <- file.path(root, "publication", "paper1", "figures", "spatial_hotspot_fused_vs_glmnet.png")
  plot_hotspots(d, fits, fig_path)
  print(summary)
}

main()
