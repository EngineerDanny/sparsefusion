library(sparsefusion)

set.seed(1)

k <- 4L
p <- 8L
n_group <- 10L
groups <- rep(seq_len(k), each = n_group)

X <- matrix(rnorm(length(groups) * p), nrow = length(groups), ncol = p)

beta <- matrix(0, nrow = p, ncol = k)
beta[1:2, 1:2] <- 1
beta[1:2, 3:4] <- -1

y <- rowSums(X * t(beta[, groups, drop = FALSE])) +
  rnorm(length(groups), sd = 0.1)

G <- matrix(0, k, k)
for (i in seq_len(k - 1L)) {
  G[i, i + 1L] <- 1
  G[i + 1L, i] <- 1
}

G

fit <- sparse_fusion(
  X, y, groups, G,
  lambda = 1e-3,
  gamma = 1e-2,
  fusion = "l1",
  solver = "active_edge",
  mu = 1e-4,
  tol = 1e-3,
  num.it = 800,
  intercept = FALSE,
  scaling = FALSE
)

dim(fit)

X_new <- X[1:6, , drop = FALSE]
groups_new <- groups[1:6]

y_hat <- numeric(nrow(X_new))
for (g in unique(groups_new)) {
  idx <- groups_new == g
  y_hat[idx] <- X_new[idx, , drop = FALSE] %*% fit[, g]
}

y_hat

# # Exact sparse graph representation, recommended default.
# sparse_fusion(X, y, groups, G, lambda, gamma,
#   fusion = "l1", solver = "active_edge"
# )
# 
# # Full all-pair reference construction.
# sparse_fusion(X, y, groups, G, lambda, gamma,
#   fusion = "l1", solver = "full_pairwise"
# )
# 
# # Approximate L1 chain solver.
# sparse_fusion(X, y, groups, G, lambda, gamma,
#   fusion = "l1", solver = "chain_approx"
# )

fit_l2 <- sparse_fusion(
  X, y, groups, G,
  lambda = 1e-3,
  gamma = 1e-2,
  fusion = "l2",
  solver = "active_edge",
  intercept = FALSE,
  scaling = FALSE
)

dim(fit_l2)
