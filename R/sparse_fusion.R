#' Fit grouped sparse fusion regression.
#'
#' @param X Feature matrix with one row per observation.
#' @param y Response vector.
#' @param groups Group label for each observation.
#' @param G Symmetric fusion graph matrix over groups.
#' @param lambda Sparsity penalty.
#' @param gamma Fusion penalty.
#' @param fusion Fusion penalty family, either `"l1"` or `"l2"`.
#' @param solver Solver representation. Use `"active_edge"` for the sparse
#'   graph solver, `"full_pairwise"` for the reference all-pair construction,
#'   or `"chain_approx"` for the L1 chain approximation.
#' @param intercept Whether to include one intercept per group.
#' @param scaling Whether to scale group losses.
#' @param ... Additional arguments passed to the selected solver.
#'
#' @return Coefficient matrix with one column per group.
#' @export
sparse_fusion <- function(
    X,
    y,
    groups,
    G,
    lambda,
    gamma,
    fusion = c("l1", "l2"),
    solver = NULL,
    intercept = FALSE,
    scaling = FALSE,
    ...) {
  fusion <- match.arg(tolower(fusion), choices = c("l1", "l2"))

  if (is.null(solver)) {
    solver <- "active_edge"
  }
  solver <- tolower(solver)
  solver <- switch(
    solver,
    operator = "active_edge",
    reference = "full_pairwise",
    full_edge = "full_pairwise",
    solver
  )
  solver <- match.arg(solver, choices = c("active_edge", "full_pairwise", "chain_approx"))

  if (fusion == "l2" && solver == "chain_approx") {
    stop("solver = \"chain_approx\" is only available for fusion = \"l1\".")
  }

  if (fusion == "l1") {
    if (solver == "full_pairwise") {
      return(fusedLassoProximal(
        X = X,
        Y = y,
        groups = groups,
        lambda = lambda,
        gamma = gamma,
        G = G,
        intercept = intercept,
        scaling = scaling,
        ...
      ))
    }

    l1_solver <- if (solver == "active_edge") "operator" else "chain_approx"
    return(fusedLassoProximalNew(
      X = X,
      Y = y,
      groups = groups,
      lambda = lambda,
      gamma = gamma,
      G = G,
      intercept = intercept,
      scaling = scaling,
      solver = l1_solver,
      ...
    ))
  }

  if (solver == "full_pairwise") {
    return(fusedL2DescentGLMNet(
      X = X,
      y = y,
      groups = groups,
      lambda = lambda,
      gamma = gamma,
      G = G,
      intercept = intercept,
      scaling = scaling,
      ...
    ))
  }

  fusedL2DescentGLMNetNew(
    X = X,
    y = y,
    groups = groups,
    lambda = lambda,
    gamma = gamma,
    G = G,
    intercept = intercept,
    scaling = scaling,
    ...
  )
}
