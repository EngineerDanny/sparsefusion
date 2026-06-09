#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    return(default)
  }
  args[[i + 1L]]
}

l1_points_csv <- get_arg(
  "--l1_points_csv",
  "build/hpc/run_20260225_123833_strict_fit_mem_knn/hpc_l1_efficiency_diff_operator_selected10_points.csv"
)
l1_fold_csv <- get_arg(
  "--l1_fold_csv",
  "build/hpc/run_20260225_123833_strict_fit_mem_knn/hpc_l1_efficiency_diff_operator_selected10_fold_points.csv"
)
l2_points_csv <- get_arg(
  "--l2_points_csv",
  "build/hpc/l2_knn_solver_intercept/hpc_l2_efficiency_diff_points_selected10_like_l1op.csv"
)
l2_fold_csv <- get_arg(
  "--l2_fold_csv",
  "build/hpc/l2_knn_solver_intercept/hpc_l2_efficiency_diff_fold_points_selected10_like_l1op.csv"
)
meta_csv <- get_arg("--meta_csv", "data/processed/grouped_regression_datasets_summary.csv")
out_png <- get_arg("--out_png", "inst/figures/hpc_l1_l2_combined_pareto.png")
out_points_csv <- get_arg("--out_points_csv", "build/hpc/hpc_l1_l2_combined_pareto_points.csv")
out_width <- as.numeric(get_arg("--width", "14"))
out_height <- as.numeric(get_arg("--height", "5.6"))

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_points_csv), recursive = TRUE, showWarnings = FALSE)

read_required <- function(path, required) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  x <- read.csv(path, stringsAsFactors = FALSE)
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop("Missing column(s) in ", path, ": ", paste(missing, collapse = ", "))
  }
  x
}

clean_dataset_name <- function(x) {
  sub("^rdatasets_", "", x)
}

format_pct <- function(x) {
  paste0(ifelse(x >= 0, "+", ""), sprintf("%.2f%%", x))
}

make_points <- function(path, panel, method_label) {
  x <- read_required(path, c(
    "dataset", "delta_log10_time", "delta_log10_memory",
    "sd_delta_log10_time", "sd_delta_log10_memory", "pct_rmse_diff", "rmse_dir"
  ))
  x$panel <- panel
  x$method_label <- method_label
  x
}

make_fold_points <- function(path, panel, method_label) {
  x <- read_required(path, c("dataset", "delta_log10_time", "delta_log10_memory", "rmse_dir"))
  x$panel <- panel
  x$method_label <- method_label
  x
}

points <- rbind(
  make_points(l1_points_csv, "L1 Active-Edge", "active_edge (ours)"),
  make_points(l2_points_csv, "L2 Active-Edge", "active_edge (ours)")
)
fold_points <- rbind(
  make_fold_points(l1_fold_csv, "L1 Active-Edge", "active_edge (ours)"),
  make_fold_points(l2_fold_csv, "L2 Active-Edge", "active_edge (ours)")
)

meta <- read_required(meta_csv, c("dataset", "n", "p", "k"))
points <- merge(points, meta[, c("dataset", "n", "p", "k"), drop = FALSE], by = "dataset", all.x = TRUE)
if (anyNA(points$n) || anyNA(points$p) || anyNA(points$k)) {
  bad <- unique(points$dataset[is.na(points$n) | is.na(points$p) | is.na(points$k)])
  stop("Missing metadata for dataset(s): ", paste(bad, collapse = ", "))
}

points$dataset_label <- clean_dataset_name(points$dataset)
points$label <- sprintf(
  "%s (%s, %s, %s, %s)",
  points$dataset_label,
  format(points$n, trim = TRUE, scientific = FALSE),
  format(points$p, trim = TRUE, scientific = FALSE),
  format(points$k, trim = TRUE, scientific = FALSE),
  format_pct(points$pct_rmse_diff)
)
points$rmse_dir <- factor(points$rmse_dir, levels = c("better_or_equal_rmse", "worse_rmse"))
fold_points$rmse_dir <- factor(fold_points$rmse_dir, levels = c("better_or_equal_rmse", "worse_rmse"))
points$panel <- factor(points$panel, levels = c("L1 Active-Edge", "L2 Active-Edge"))
fold_points$panel <- factor(fold_points$panel, levels = levels(points$panel))
points$dataset_label <- factor(points$dataset_label, levels = sort(unique(points$dataset_label)))
fold_points$dataset_label <- clean_dataset_name(fold_points$dataset)
fold_points$dataset_label <- factor(fold_points$dataset_label, levels = levels(points$dataset_label))

write.csv(points, out_points_csv, row.names = FALSE)

p <- ggplot(points, aes(x = delta_log10_time, y = delta_log10_memory)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
  geom_point(
    data = fold_points,
    aes(x = delta_log10_time, y = delta_log10_memory, color = dataset_label),
    shape = 16,
    size = 1.15,
    alpha = 0.18,
    inherit.aes = FALSE
  ) +
  geom_segment(
    aes(
      x = delta_log10_time - sd_delta_log10_time,
      xend = delta_log10_time + sd_delta_log10_time,
      y = delta_log10_memory,
      yend = delta_log10_memory,
      color = dataset_label
    ),
    linewidth = 0.55,
    alpha = 0.95
  ) +
  geom_segment(
    aes(
      x = delta_log10_time,
      xend = delta_log10_time,
      y = delta_log10_memory - sd_delta_log10_memory,
      yend = delta_log10_memory + sd_delta_log10_memory,
      color = dataset_label
    ),
    linewidth = 0.4,
    alpha = 0.95
  ) +
  geom_point(
    aes(color = dataset_label),
    shape = 21,
    fill = NA,
    stroke = 0.85,
    size = 2.6,
    alpha = 0.98,
    show.legend = FALSE
  ) +
  ggrepel::geom_label_repel(
    aes(label = label, color = dataset_label),
    size = 2.8,
    label.size = 0.18,
    fill = "white",
    alpha = 0.88,
    box.padding = 0.25,
    point.padding = 0.14,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.alpha = 0.82,
    segment.size = 0.32,
    show.legend = FALSE
  ) +
  facet_wrap(~panel, nrow = 1, scales = "free_x", strip.position = "top") +
  labs(
    x = "log(Time / Reference Time)\n(left is faster)",
    y = "log(Memory / Reference Memory)\n(down is lower memory)",
    caption = "Labels: dataset (n, p, k, ΔRMSE%). ΔRMSE% is relative to the reference; negative is better."
  ) +
  theme_grey(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_line(color = "white", linewidth = 0.35),
    panel.grid.major = element_line(color = "white", linewidth = 0.65),
    strip.text = element_text(size = 15, color = "grey10"),
    strip.background = element_rect(fill = "#D9D9D9", color = "#D9D9D9"),
    axis.text = element_text(size = 12, color = "grey30"),
    axis.title.x = element_text(size = 16, color = "black", lineheight = 0.9),
    axis.title.y = element_text(size = 16, color = "black", lineheight = 0.9),
    plot.caption = element_text(size = 11, color = "black", hjust = 0.5, lineheight = 0.95),
    panel.spacing.y = unit(0.9, "lines"),
    plot.margin = margin(8, 16, 12, 16)
  )

ggsave(out_png, plot = p, width = out_width, height = out_height, dpi = 500)

cat("Saved combined Pareto figure: ", out_png, "\n", sep = "")
cat("Saved combined Pareto points: ", out_points_csv, "\n", sep = "")
