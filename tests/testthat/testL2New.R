library(fuserplus)

make_tiny_l2_data <- function(seed = 42L, k = 3L, n_group = 10L, p = 5L, sigma = 0.1) {
  set.seed(seed)
  groups <- rep(seq_len(k), each = n_group)
  X <- matrix(rnorm(length(groups) * p), nrow = length(groups), ncol = p)
  beta <- matrix(rnorm(p * k), nrow = p, ncol = k)
  y <- rowSums(X * t(beta[, groups, drop = FALSE])) + stats::rnorm(length(groups), sd = sigma)
  list(X = X, y = y, groups = groups)
}

test_that("L2-new validates scalar lambda and gamma", {
  d <- make_tiny_l2_data()
  G <- matrix(1, 3, 3)
  diag(G) <- 0

  expect_error(
    fuserplus:::fusedL2DescentGLMNetNew(d$X, d$y, d$groups, lambda = NULL, G = G, gamma = 1e-3),
    "lambda"
  )
  expect_error(
    fuserplus:::fusedL2DescentGLMNetNew(d$X, d$y, d$groups, lambda = c(1e-3, 1e-2), G = G, gamma = 1e-3),
    "lambda"
  )
  expect_error(
    fuserplus:::fusedL2DescentGLMNetNew(d$X, d$y, d$groups, lambda = 1e-3, G = G, gamma = NA_real_),
    "gamma"
  )
})

test_that("L2-new validates G and symmetrizes asymmetric G", {
  d <- make_tiny_l2_data()
  G_bad <- matrix(1, 3, 3)
  diag(G_bad) <- 0
  G_bad[1, 2] <- -1

  expect_error(
    fuserplus:::fusedL2DescentGLMNetNew(d$X, d$y, d$groups, lambda = 1e-3, G = G_bad, gamma = 1e-3),
    "non-negative"
  )

  G_asym <- matrix(0, 3, 3)
  G_asym[1, 2] <- 1
  G_asym[2, 1] <- 0.2
  G_asym[2, 3] <- 1
  G_asym[3, 2] <- 1

  expect_warning(
    fit <- fuserplus:::fusedL2DescentGLMNetNew(d$X, d$y, d$groups, lambda = 1e-3, G = G_asym, gamma = 1e-3),
    "not symmetric"
  )
  expect_true(is.matrix(fit))
  expect_equal(dim(fit), c(ncol(d$X), length(unique(d$groups))))
})

test_that("L2-new returns expected coefficient dimensions with and without intercept", {
  d <- make_tiny_l2_data()
  G <- matrix(1, 3, 3)
  diag(G) <- 0

  fit_no_intercept <- fuserplus:::fusedL2DescentGLMNetNew(
    d$X, d$y, d$groups,
    lambda = 1e-3,
    G = G,
    gamma = 1e-3,
    intercept = FALSE
  )
  expect_equal(dim(fit_no_intercept), c(ncol(d$X), length(unique(d$groups))))

  fit_intercept <- fuserplus:::fusedL2DescentGLMNetNew(
    d$X, d$y, d$groups,
    lambda = 1e-3,
    G = G,
    gamma = 1e-3,
    intercept = TRUE
  )
  expect_equal(dim(fit_intercept), c(ncol(d$X) + 1L, length(unique(d$groups))))
})
