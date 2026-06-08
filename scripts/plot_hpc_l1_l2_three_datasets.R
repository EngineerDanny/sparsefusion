#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    return(default)
  }
  args[[i + 1L]]
}

l1_csv <- get_arg("--l1_csv", "build/hpc/run_20260225_123833_strict_fit_mem_knn/final_test_summary.csv")
l2_csv <- get_arg("--l2_csv", "build/hpc/l2_knn_solver_intercept/final_test_summary.csv")
meta_csv <- get_arg("--meta_csv", "data/processed/grouped_regression_datasets_summary.csv")
out_png <- get_arg("--out_png", "inst/figures/hpc_l1_l2_three_datasets.png")
out_width <- as.numeric(get_arg("--out_width", "12"))
out_height <- as.numeric(get_arg("--out_height", "4"))

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

datasets_keep <- c("communities_crime", "rdatasets_hsb82", "owid_co2")

read_checked <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

clean_dataset_name <- function(x) {
  sub("^rdatasets_", "", x)
}

format_method_axis_labels <- function(x) {
  label_expr <- ifelse(
    x %in% c("active_edge (ours)", "chain_approx (ours)"),
    sprintf("bold('%s')", x),
    sprintf("'%s'", x)
  )
  parse(text = label_expr)
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

prep_l1 <- function(path) {
  x <- read_checked(path)
  x <- x[x$dataset %in% datasets_keep & x$method %in% c("operator", "old_l1", "chain_approx", "featureless"), , drop = FALSE]
  method_map <- c(
    operator = "active_edge (ours)",
    old_l1 = "full_pairwise",
    chain_approx = "chain_approx (ours)",
    featureless = "featureless"
  )
  x$panel <- "L1"
  x$method_label <- unname(method_map[x$method])
  x
}

prep_l2 <- function(path) {
  x <- read_checked(path)
  x <- x[x$dataset %in% datasets_keep & x$method %in% c("new_l2", "old_l2", "featureless"), , drop = FALSE]
  method_map <- c(
    new_l2 = "active_edge (ours)",
    old_l2 = "full_pairwise",
    featureless = "featureless"
  )
  x$panel <- "L2"
  x$method_label <- unname(method_map[x$method])
  x
}

meta <- read_checked(meta_csv)
meta <- meta[meta$dataset %in% datasets_keep, c("dataset", "n", "p", "k"), drop = FALSE]
if (nrow(meta) != length(datasets_keep)) {
  stop("Metadata missing for one or more requested datasets.")
}

l1 <- prep_l1(l1_csv)
l2 <- prep_l2(l2_csv)
plot_df <- rbind(
  l1[, c("dataset", "panel", "method_label", "mean_test_rmse", "sd_test_rmse"), drop = FALSE],
  l2[, c("dataset", "panel", "method_label", "mean_test_rmse", "sd_test_rmse"), drop = FALSE]
)

required_rows <- expand.grid(
  dataset = datasets_keep,
  panel = c("L1", "L2"),
  method_label = c("active_edge (ours)", "full_pairwise", "chain_approx (ours)", "featureless"),
  stringsAsFactors = FALSE
)
required_rows <- required_rows[
  !(required_rows$panel == "L2" & required_rows$method_label == "chain_approx (ours)"),
  ,
  drop = FALSE
]
merged_check <- merge(required_rows, plot_df, by = c("dataset", "panel", "method_label"), all.x = TRUE)
if (anyNA(merged_check$mean_test_rmse) || anyNA(merged_check$sd_test_rmse)) {
  missing_rows <- merged_check[is.na(merged_check$mean_test_rmse) | is.na(merged_check$sd_test_rmse), c("dataset", "panel", "method_label")]
  print(missing_rows)
  stop("Plot inputs are incomplete.")
}

dataset_order <- c("communities_crime", "rdatasets_hsb82", "owid_co2")
panel_order <- c("L1", "L2")
method_order <- c("active_edge (ours)", "chain_approx (ours)", "full_pairwise", "featureless")

strip_map <- setNames(vapply(dataset_order, make_strip_label, character(1), meta = meta), dataset_order)
plot_df$dataset_strip <- factor(strip_map[plot_df$dataset], levels = unname(strip_map[dataset_order]))
plot_df$panel <- factor(plot_df$panel, levels = panel_order)
plot_df$method_label <- factor(plot_df$method_label, levels = rev(method_order))
plot_df$xmin <- plot_df$mean_test_rmse - plot_df$sd_test_rmse
plot_df$xmax <- plot_df$mean_test_rmse + plot_df$sd_test_rmse
plot_df$label <- sprintf("%.3f\u00b1%.3f", plot_df$mean_test_rmse, plot_df$sd_test_rmse)

panel_ranges <- aggregate(cbind(xmin, xmax) ~ panel + dataset_strip, data = plot_df, FUN = range)
panel_ranges$xmin <- panel_ranges$xmin[, 1L]
panel_ranges$xmax <- panel_ranges$xmax[, 2L]
panel_ranges$pad <- pmax((panel_ranges$xmax - panel_ranges$xmin) * 0.18, 0.02)
plot_df <- merge(plot_df, panel_ranges, by = c("panel", "dataset_strip"), suffixes = c("", "_range"))
plot_df$label_x <- ifelse(
  plot_df$method_label == "featureless",
  plot_df$xmin - plot_df$pad * 0.04,
  plot_df$xmax + plot_df$pad * 0.06
)
plot_df$label_hjust <- ifelse(plot_df$method_label == "featureless", 1, 0)
plot_df$blank_x_left <- plot_df$xmin_range - plot_df$pad
plot_df$blank_x <- plot_df$xmax_range + plot_df$pad

p <- ggplot(plot_df, aes(x = mean_test_rmse, y = method_label)) +
  geom_blank(aes(x = blank_x_left)) +
  geom_blank(aes(x = blank_x)) +
  geom_vline(
    data = subset(plot_df, method_label == "active_edge (ours)"),
    aes(xintercept = mean_test_rmse),
    inherit.aes = FALSE,
    color = "red",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_segment(aes(x = xmin, xend = xmax, yend = method_label), linewidth = 0.5, color = "grey25") +
  geom_point(shape = 21, size = 2.3, stroke = 0.5, fill = "white", color = "grey20") +
  geom_text(
    aes(x = label_x, label = label, hjust = label_hjust),
    size = 4.0,
    color = "grey20"
  ) +
  facet_grid(panel ~ dataset_strip, scales = "free", space = "free_y") +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  scale_y_discrete(labels = format_method_axis_labels) +
  labs(x = "Test RMSE (mean \u00b1 SD) over 5 CV Folds", y = NULL) +
  theme_grey(base_size = 16) +
  theme(
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
    plot.margin = margin(12, 52, 12, 12)
  )

ggsave(out_png, plot = p, width = out_width, height = out_height, dpi = 300)
cat(sprintf("Saved figure: %s\n", out_png))
