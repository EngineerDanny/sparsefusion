library(sparsefusion)

make_sparse_fusion_data <- function(seed = 7L, k = 3L, n_group = 6L, p = 4L) {
  set.seed(seed)
  groups <- rep(seq_len(k), each = n_group)
  X <- matrix(stats::rnorm(length(groups) * p), nrow = length(groups), ncol = p)
  beta <- matrix(0, nrow = p, ncol = k)
  beta[1, ] <- c(1, 1, -1)
  beta[2, ] <- c(0.5, -0.5, -0.5)
  y <- rowSums(X * t(beta[, groups, drop = FALSE])) + stats::rnorm(length(groups), sd = 0.05)
  G <- matrix(0, k, k)
  for (i in seq_len(k - 1L)) {
    G[i, i + 1L] <- 1
    G[i + 1L, i] <- 1
  }
  list(X = X, y = y, groups = groups, G = G)
}

expect_coef_shape <- function(fit, d) {
  expect_true(is.matrix(fit))
  expect_equal(dim(fit), c(ncol(d$X), length(unique(d$groups))))
}

suppress_max_iter_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("Reached max iterations without convergence", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

test_that("sparse_fusion dispatches L1 solver modes", {
  d <- make_sparse_fusion_data()

  fit_active <- suppress_max_iter_warning(sparse_fusion(
    d$X, d$y, d$groups, d$G,
    lambda = 1e-3, gamma = 1e-2,
    fusion = "l1", solver = "active_edge",
    mu = 1e-4, tol = 1e-3, num.it = 50
  ))
  expect_coef_shape(fit_active, d)

  fit_full <- suppress_max_iter_warning(sparse_fusion(
    d$X, d$y, d$groups, d$G,
    lambda = 1e-3, gamma = 1e-2,
    fusion = "l1", solver = "full_pairwise",
    mu = 1e-4, tol = 1e-3, num.it = 50
  ))
  expect_coef_shape(fit_full, d)

  fit_chain <- suppress_max_iter_warning(sparse_fusion(
    d$X, d$y, d$groups, d$G,
    lambda = 1e-3, gamma = 1e-2,
    fusion = "l1", solver = "chain_approx",
    mu = 1e-4, tol = 1e-3, num.it = 50
  ))
  expect_coef_shape(fit_chain, d)
})

test_that("sparse_fusion dispatches L2 solver modes", {
  d <- make_sparse_fusion_data()

  fit_active <- sparse_fusion(
    d$X, d$y, d$groups, d$G,
    lambda = 1e-3, gamma = 1e-2,
    fusion = "l2", solver = "active_edge"
  )
  expect_coef_shape(fit_active, d)

  fit_full <- sparse_fusion(
    d$X, d$y, d$groups, d$G,
    lambda = 1e-3, gamma = 1e-2,
    fusion = "l2", solver = "full_pairwise"
  )
  expect_coef_shape(fit_full, d)

  expect_error(
    sparse_fusion(
      d$X, d$y, d$groups, d$G,
      lambda = 1e-3, gamma = 1e-2,
      fusion = "l2", solver = "chain_approx"
    ),
    "only available"
  )
})
