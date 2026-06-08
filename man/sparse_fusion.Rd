\name{sparse_fusion}
\alias{sparse_fusion}
\title{Fit grouped sparse fusion regression}
\usage{
sparse_fusion(
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
  ...
)
}
\arguments{
\item{X}{Feature matrix with one row per observation.}
\item{y}{Response vector.}
\item{groups}{Group label for each observation.}
\item{G}{Symmetric fusion graph matrix over groups.}
\item{lambda}{Sparsity penalty.}
\item{gamma}{Fusion penalty.}
\item{fusion}{Fusion penalty family, either \code{"l1"} or \code{"l2"}.}
\item{solver}{Solver representation. Use \code{"active_edge"} for the sparse graph solver, \code{"full_pairwise"} for the reference all-pair construction, or \code{"chain_approx"} for the L1 chain approximation.}
\item{intercept}{Whether to include one intercept per group.}
\item{scaling}{Whether to scale group losses.}
\item{...}{Additional arguments passed to the selected solver.}
}
\value{
Coefficient matrix with one column per group.
}
\description{
Fits grouped fused-regression models through the recommended sparsefusion
interface. The function dispatches to the L1 or L2 solver family and selects
the requested graph representation.
}
