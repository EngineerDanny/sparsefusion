# sparsefusion

`sparsefusion` fits grouped fused-regression models with L1 and L2 fusion
penalties. The package extends the original grouped fused-regression workflow
with sparse graph solver variants that avoid building all-pair fusion objects
when the fusion graph is sparse.

The paper focuses on one computational question:

> If the statistical penalty is defined on a sparse graph, can the solver keep
> that sparse edge representation instead of expanding to all group pairs?

## Install

From a checkout of this repository:

```bash
R CMD INSTALL .
```

Or from GitHub:

```r
install.packages("remotes")
remotes::install_github("EngineerDanny/sparsefusion")
```

Optional packages used by the timing scripts:

```r
install.packages(c("atime", "ggplot2"))
```

## Minimal Reviewer Example

This example creates a small grouped regression problem with a sparse chain
fusion graph. It fits L1 and L2 grouped fused-regression models through the
recommended `sparse_fusion()` interface.

```r
library(sparsefusion)

set.seed(1)

k <- 4L
p <- 8L
n_group <- 8L
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

fit_l1_ref <- sparse_fusion(
  X, y, groups,
  lambda = 1e-3, gamma = 1e-2, G = G,
  fusion = "l1", solver = "full_pairwise",
  mu = 1e-4, tol = 1e-3, num.it = 800,
  intercept = FALSE, scaling = FALSE
)

fit_l1_active <- sparse_fusion(
  X, y, groups,
  lambda = 1e-3, gamma = 1e-2, G = G,
  fusion = "l1", solver = "active_edge",
  mu = 1e-4, tol = 1e-3, num.it = 800,
  intercept = FALSE, scaling = FALSE
)

fit_l1_chain <- sparse_fusion(
  X, y, groups,
  lambda = 1e-3, gamma = 1e-2, G = G,
  fusion = "l1", solver = "chain_approx",
  mu = 1e-4, tol = 1e-3, num.it = 800,
  intercept = FALSE, scaling = FALSE
)

fit_l2_ref <- sparse_fusion(
  X, y, groups,
  lambda = 1e-3, gamma = 1e-2, G = G,
  fusion = "l2", solver = "full_pairwise",
  scaling = FALSE
)

fit_l2_active <- sparse_fusion(
  X, y, groups,
  lambda = 1e-3, gamma = 1e-2, G = G,
  fusion = "l2", solver = "active_edge",
  scaling = FALSE
)

dim(fit_l1_ref)
dim(fit_l1_active)
dim(fit_l1_chain)
dim(fit_l2_ref)
dim(fit_l2_active)
```

Expected dimensions are `p` by `k`, here `8` by `4`.

## Conceptual Timing Panels

The conceptual figures in the paper use timing panels that compare full-pairwise
reference construction with active-edge construction. The labels and colors are
defined in the scripts:

- `Full-Pairwise Ref.` uses `#FA786E`
- `Active-Edge (Ours)` uses `#64A0FF`
- `Chain Approx. (Ours)` uses `#8E63D9`

Run the L1 timing panel from the repository root:

```bash
Rscript scripts/benchmark_atime_l1_variants.R \
  --data_mode synthetic \
  --k_values 10,20,30,40,50,80,100 \
  --p 120 \
  --n_group_train 35 \
  --g_structure sparse_chain \
  --times 3 \
  --out_prefix atime_l1_variants
```

This writes:

```text
build/atime_l1_variants_summary.csv
build/atime_l1_variants_atime_obj.rds
inst/figures/atime_l1_variants_atime.png
```

Run the L2 timing panel:

```bash
Rscript scripts/benchmark_atime_l2_variants.R \
  --data_mode synthetic \
  --k_values 10,20,30,40,50,80,100 \
  --p 120 \
  --n_group_train 35 \
  --g_structure sparse_chain \
  --times 3 \
  --out_prefix atime_l2
```

This writes:

```text
build/atime_l2_l2_summary.csv
build/atime_l2_l2_atime_obj.rds
inst/figures/atime_l2_l2_atime.png
```

The L1 run can be slow because the full-pairwise reference is intentionally
included at larger numbers of groups. For a quick smoke test of the plotting
pipeline, reduce the grid and repeats:

```bash
Rscript scripts/benchmark_atime_l1_variants.R \
  --data_mode synthetic \
  --k_values 10,20,30 \
  --p 120 \
  --n_group_train 35 \
  --g_structure sparse_chain \
  --times 1 \
  --out_prefix atime_l1_smoke
```


## Tests

Run the package tests from the repository root:

```bash
R -q -e "testthat::test_dir('tests/testthat')"
```

## Notes On Solver Variants

- `sparse_fusion(..., fusion = "l1", solver = "active_edge")` runs the exact
  L1 active-edge operator.
- `sparse_fusion(..., fusion = "l1", solver = "full_pairwise")` runs the L1
  full-pairwise reference.
- `sparse_fusion(..., fusion = "l1", solver = "chain_approx")` runs the
  approximate L1 chain solver.
- `sparse_fusion(..., fusion = "l2", solver = "active_edge")` runs the L2
  active-edge augmented-design builder.
- `sparse_fusion(..., fusion = "l2", solver = "full_pairwise")` runs the L2
  full-pairwise augmented-design reference.

For a longer modeling walkthrough, see `vignettes/subgroup_fusion.Rmd`.
