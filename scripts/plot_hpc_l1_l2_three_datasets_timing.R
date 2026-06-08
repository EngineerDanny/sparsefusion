#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    return(default)
  }
  args[[i + 1L]]
}

l1_csv <- get_arg("--l1_csv", "build/hpc/run_20260225_123833_strict_fit_mem_knn/test_results_all.csv")
l2_csv <- get_arg("--l2_csv", "build/hpc/l2_knn_solver_intercept/test_results_all.csv")
meta_csv <- get_arg("--meta_csv", "data/processed/grouped_regression_datasets_summary.csv")
out_png <- get_arg("--out_png", "inst/figures/hpc_l1_l2_three_datasets_timing.png")
out_width <- as.numeric(get_arg("--out_width", "12"))
out_height <- as.numeric(get_arg("--out_height", "4"))

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

datasets_keep <- c("communities_crime", "rdatasets_hsb82", "owid_co2")
dataset_order <- c("communities_crime", "rdatasets_hsb82", "owid_co2")

read_checked <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

clean_dataset_name <- function(x) {
  sub("^rdatasets_", "", x)
}

make_strip_label <- function(dataset, meta) {
  idx <- match(dataset, meta$dataset)
  sprintf(
    "%s\nn=%s, p=%s, k=%s",
    clean_dataset_name(dataset),
    format(meta$n[idx], big.mark = "", scientific = FALSE, trim = TRUE),
    format(meta$p[idx], big.mark = "", scientific = FALSE, trim = TRUE),
    format(meta$k[idx], big.mark = "", scientific = FALSE, trim = TRUE)
  )
}

format_method_axis_labels <- function(x) {
  label_expr <- ifelse(
    grepl("\\(ours\\)", x),
    sprintf("bold('%s')", x),
    sprintf("'%s'", x)
  )
  parse(text = label_expr)
}

format_seconds_label <- function(x) {
  ifelse(
    x < 0.01,
    sprintf("%.3fs", x),
    ifelse(x < 1, sprintf("%.2fs", x), ifelse(x < 10, sprintf("%.2fs", x), sprintf("%.1fs", x)))
  )
}

prep_elapsed <- function(path, panel, method_map) {
  x <- read_checked(path)
  needed <- c("dataset", "method", "status", "elapsed_sec")
  if (!all(needed %in% names(x))) {
    stop("Timing CSV missing columns: ", paste(setdiff(needed, names(x)), collapse = ", "))
  }
  x <- x[
    x$dataset %in% datasets_keep &
      x$method %in% names(method_map) &
      x$status == "ok" &
      is.finite(x$elapsed_sec),
    ,
    drop = FALSE
  ]
  x$elapsed_sec <- as.numeric(x$elapsed_sec)
  x <- x[is.finite(x$elapsed_sec) & x$elapsed_sec > 0, , drop = FALSE]
  if (!nrow(x)) {
    stop("No timing rows left for ", panel)
  }
  x$log10_elapsed_ms <- log10(1000 * x$elapsed_sec)
  x$panel <- panel
  x$method_label <- unname(method_map[x$method])
  x
}

summarize_elapsed <- function(x, strip_map) {
  mean_df <- aggregate(
    log10_elapsed_ms ~ dataset + panel + method_label,
    data = x[, c("dataset", "panel", "method_label", "log10_elapsed_ms"), drop = FALSE],
    FUN = mean
  )
  names(mean_df)[names(mean_df) == "log10_elapsed_ms"] <- "mean_log10_elapsed_ms"
  sd_df <- aggregate(
    log10_elapsed_ms ~ dataset + panel + method_label,
    data = x[, c("dataset", "panel", "method_label", "log10_elapsed_ms"), drop = FALSE],
    FUN = stats::sd
  )
  names(sd_df)[names(sd_df) == "log10_elapsed_ms"] <- "sd_log10_elapsed_ms"
  out <- merge(mean_df, sd_df, by = c("dataset", "panel", "method_label"), all.x = TRUE, sort = FALSE)
  out$sd_log10_elapsed_ms[!is.finite(out$sd_log10_elapsed_ms)] <- 0
  out$xmin <- pmax(0, out$mean_log10_elapsed_ms - out$sd_log10_elapsed_ms)
  out$xmax <- out$mean_log10_elapsed_ms + out$sd_log10_elapsed_ms
  out$dataset_strip <- factor(strip_map[out$dataset], levels = unname(strip_map[dataset_order]))
  out
}

check_required <- function(x, panel, required_methods) {
  required <- expand.grid(
    dataset = datasets_keep,
    method_label = required_methods,
    stringsAsFactors = FALSE
  )
  merged <- merge(required, x, by = c("dataset", "method_label"), all.x = TRUE)
  if (anyNA(merged$mean_log10_elapsed_ms)) {
    missing <- merged[is.na(merged$mean_log10_elapsed_ms), c("dataset", "method_label")]
    print(missing)
    stop("Timing inputs are incomplete for ", panel)
  }
}

meta <- read_checked(meta_csv)
meta <- meta[meta$dataset %in% datasets_keep, c("dataset", "n", "p", "k"), drop = FALSE]
if (nrow(meta) != length(datasets_keep)) {
  stop("Metadata missing for one or more requested datasets.")
}
strip_map <- setNames(vapply(dataset_order, make_strip_label, character(1), meta = meta), dataset_order)

l1_method_map <- c(
  operator = "active_edge (ours)",
  chain_approx = "chain_approx (ours)",
  old_l1 = "full_pairwise"
)
l2_method_map <- c(
  new_l2 = "active_edge (ours)",
  old_l2 = "full_pairwise"
)

l1_raw <- prep_elapsed(l1_csv, "L1", l1_method_map)
l2_raw <- prep_elapsed(l2_csv, "L2", l2_method_map)
l1 <- summarize_elapsed(l1_raw, strip_map)
l2 <- summarize_elapsed(l2_raw, strip_map)

l1_order <- c("active_edge (ours)", "chain_approx (ours)", "full_pairwise")
l2_order <- c("active_edge (ours)", "full_pairwise")
check_required(l1, "L1", l1_order)
check_required(l2, "L2", l2_order)
plot_df <- rbind(l1, l2)
plot_df$panel <- factor(plot_df$panel, levels = c("L1", "L2"))
plot_df$method_label <- factor(
  plot_df$method_label,
  levels = rev(c("active_edge (ours)", "chain_approx (ours)", "full_pairwise"))
)

p <- ggplot(plot_df, aes(x = mean_log10_elapsed_ms, y = method_label)) +
  geom_col(width = 0.62, na.rm = TRUE) +
  geom_errorbar(
    aes(xmin = xmin, xmax = xmax),
    orientation = "y",
    width = 0.25,
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  scale_y_discrete(drop = TRUE, labels = format_method_axis_labels) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
  facet_grid(panel ~ dataset_strip, scales = "free", space = "free_y") +
  labs(x = "log10 elapsed time (milliseconds) over 5 CV Folds (lower is better)", y = NULL) +
  theme_grey(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.7),
    panel.grid.major.x = element_line(color = "white", linewidth = 0.5),
    axis.text.y = element_text(color = "grey30"),
    axis.text.x = element_text(size = 14, color = "grey30"),
    strip.text = element_text(size = 16, color = "grey15", lineheight = 0.95),
    strip.background = element_rect(fill = "#D9D9D9", color = "#D9D9D9"),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    axis.title.x = element_text(size = 19, color = "grey10"),
    plot.margin = margin(12, 12, 12, 12)
  )

ggsave(out_png, plot = p, width = out_width, height = out_height, dpi = 300)

cat(sprintf("Saved timing figure: %s\n", out_png))
