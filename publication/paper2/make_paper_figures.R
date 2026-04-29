#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
fig_dir <- file.path(root, "publication", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

short_dataset <- function(x) {
  map <- c(
    communities_crime = "Communities",
    school_ilea = "School",
    king_county_house_sales = "King County",
    beijing_air_quality = "Beijing Air",
    owid_co2 = "OWID CO2",
    radon_minnesota = "Radon",
    rdatasets_grunfeld = "Grunfeld",
    rdatasets_dietox = "Dietox",
    world_bank_wdi = "WDI"
  )
  out <- unname(map[x])
  ifelse(is.na(out), x, out)
}

theme_paper <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

make_conceptual_figure <- function() {
  png(file.path(fig_dir, "fusetune_concept.png"), width = 2400, height = 1200, res = 220)
  grid.newpage()

  draw_box <- function(x, y, w, h, label = NULL, fill = "#f8f9fa", col = "#495057",
                       fontsize = 13, lwd = 2, fontface = "plain", just = "centre") {
    grid.roundrect(
      x = unit(x, "npc"), y = unit(y, "npc"),
      width = unit(w, "npc"), height = unit(h, "npc"),
      r = unit(0.02, "snpc"),
      gp = gpar(fill = fill, col = col, lwd = lwd)
    )
    if (!is.null(label)) {
      grid.text(label, x = unit(x, "npc"), y = unit(y, "npc"), just = just,
                gp = gpar(col = "#212529", fontsize = fontsize, fontface = fontface))
    }
  }

  draw_arrow <- function(x0, y0, x1, y1, col = "#495057", lwd = 2) {
    grid.segments(
      x0 = unit(x0, "npc"), y0 = unit(y0, "npc"),
      x1 = unit(x1, "npc"), y1 = unit(y1, "npc"),
      arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
      gp = gpar(col = col, lwd = lwd)
    )
  }

  palette <- list(
    ink = "#1f2933",
    muted = "#52606d",
    line = "#9aa5b1",
    panel = "#f8fafc",
    blue = "#1d4ed8",
    blue_fill = "#eff6ff",
    red = "#b42318",
    red_fill = "#fff1f2",
    gold = "#a16207",
    gold_fill = "#fffbeb"
  )

  grid.text("FuseTune: solver-aware low-budget tuning for fused regression",
            x = unit(0.5, "npc"), y = unit(0.92, "npc"),
            gp = gpar(fontsize = 15.5, fontface = "bold", col = palette$ink))
  grid.text("Shared objective, branch-specific search policy",
            x = unit(0.5, "npc"), y = unit(0.875, "npc"),
            gp = gpar(fontsize = 10.8, fontface = "italic", col = palette$muted))

  draw_box(0.15, 0.60, 0.19, 0.17,
           fill = "#ffffff", col = palette$line, lwd = 1.5)
  grid.text("Inputs", x = unit(0.15, "npc"), y = unit(0.66, "npc"),
            gp = gpar(fontsize = 12.2, fontface = "bold", col = palette$ink))
  grid.text("grouped data\nfixed folds\nsearch bounds",
            x = unit(0.15, "npc"), y = unit(0.56, "npc"),
            gp = gpar(fontsize = 10.2, col = palette$muted))

  draw_box(0.41, 0.60, 0.24, 0.19,
           fill = palette$blue_fill, col = palette$blue, lwd = 1.8)
  grid.text("Cached CV objective", x = unit(0.41, "npc"), y = unit(0.675, "npc"),
            gp = gpar(fontsize = 12.2, fontface = "bold", col = palette$blue))
  grid.text("optimize in log(lambda), log(gamma)\nvalidation RMSE with pruning and reuse",
            x = unit(0.41, "npc"), y = unit(0.57, "npc"),
            gp = gpar(fontsize = 9.9, col = palette$ink))

  draw_box(0.71, 0.58, 0.32, 0.43,
           fill = "#ffffff", col = palette$line, lwd = 1.6)
  grid.text("Solver-specific inner policy", x = unit(0.71, "npc"), y = unit(0.75, "npc"),
            gp = gpar(fontsize = 12.6, fontface = "bold", col = palette$ink))

  draw_box(0.63, 0.56, 0.14, 0.25,
           fill = palette$red_fill, col = palette$red, lwd = 1.5)
  grid.text("l1 branch", x = unit(0.63, "npc"), y = unit(0.645, "npc"),
            gp = gpar(fontsize = 11.4, fontface = "bold", col = palette$red))
  grid.text("center start\naxis-aligned direct search\nnonmonotone acceptance\nadaptive step sizes",
            x = unit(0.63, "npc"), y = unit(0.53, "npc"),
            gp = gpar(fontsize = 8.8, col = palette$ink))

  draw_box(0.79, 0.56, 0.14, 0.25,
           fill = palette$gold_fill, col = palette$gold, lwd = 1.5)
  grid.text("l2 branch", x = unit(0.79, "npc"), y = unit(0.645, "npc"),
            gp = gpar(fontsize = 11.4, fontface = "bold", col = palette$gold))
  grid.text("glmnet lambda seed\nsmall gamma screen\ngamma-first refinement\nconditional lambda touch-up",
            x = unit(0.79, "npc"), y = unit(0.53, "npc"),
            gp = gpar(fontsize = 8.8, col = palette$ink))

  draw_box(0.92, 0.60, 0.07, 0.17,
           fill = "#ffffff", col = palette$line, lwd = 1.5)
  grid.text("Output", x = unit(0.92, "npc"), y = unit(0.66, "npc"),
            gp = gpar(fontsize = 12.2, fontface = "bold", col = palette$ink))
  grid.text("selected\nhyperparameters\nfinal fit",
            x = unit(0.92, "npc"), y = unit(0.56, "npc"),
            gp = gpar(fontsize = 9.1, col = palette$muted))

  draw_arrow(0.245, 0.60, 0.29, 0.60, col = palette$muted, lwd = 1.8)
  draw_arrow(0.53, 0.60, 0.56, 0.60, col = palette$muted, lwd = 1.8)
  draw_arrow(0.87, 0.60, 0.885, 0.60, col = palette$muted, lwd = 1.8)

  grid.text("Branch choice is determined by the fitted solver, not by a second optimization layer.",
            x = unit(0.5, "npc"), y = unit(0.26, "npc"),
            gp = gpar(fontsize = 10.2, col = palette$muted))
  dev.off()
}

# Figure 1: Real-data FuseTune vs glmnet_seeded.
make_real_efficiency <- function() {
  files <- list(
    l1 = file.path(root, "build", "lit100_round3", "experiments", "real_l1_tier1_retained", "summary_long.csv"),
    l2 = file.path(root, "build", "lit100_round3", "experiments", "real_l2_tier1_retained", "summary_long.csv")
  )
  parts <- lapply(names(files), function(sv) {
    dt <- fread(files[[sv]])
    dt <- dt[method %in% c("glmnet_seeded", "fusetune")]
    wide <- dcast(
      dt,
      dataset + solver ~ method,
      value.var = c("Elapsed time (seconds)", "Test RMSE")
    )
    wide[, speedup := `Elapsed time (seconds)_glmnet_seeded` / `Elapsed time (seconds)_fusetune`]
    wide[, rmse_ratio := `Test RMSE_fusetune` / `Test RMSE_glmnet_seeded`]
    wide[, dataset_label := short_dataset(dataset)]
    wide
  })
  cmp <- rbindlist(parts, fill = TRUE)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    label_geom <- ggrepel::geom_text_repel(
      size = 2.6, max.overlaps = 20, seed = 42,
      min.segment.length = 0.3, box.padding = 0.35,
      point.padding = 0.2, segment.color = "grey70"
    )
  } else {
    label_geom <- geom_text(nudge_y = 0.008, size = 2.6, check_overlap = FALSE, hjust = 0.5)
  }
  p <- ggplot(cmp, aes(x = speedup, y = rmse_ratio, label = dataset_label)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.3, color = "#0b7285") +
    label_geom +
    facet_wrap(~ solver, scales = "free_x") +
    labs(
      title = "FuseTune versus glmnet-seeded on real datasets",
      x = "Speedup over glmnet-seeded (>1 favors FuseTune)",
      y = "Test RMSE ratio (FuseTune / glmnet-seeded)"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    theme_paper()
  ggsave(file.path(fig_dir, "real_efficiency_compare.png"), p, width = 10, height = 4.6, dpi = 300)
}

# Figure 2: Synthetic family summary, retained baselines only.
make_synth_family <- function() {
  l1 <- fread(file.path(root, "build", "synth_regime_design", "benchmark_v1", "gain_efficiency_long.csv"))[solver == "l1"]
  l2 <- fread(file.path(root, "build", "synth_regime_design", "benchmark_v1_l2_hybrid_20260401", "gain_efficiency_long.csv"))[solver == "l2"]
  dt <- rbindlist(list(l1, l2), fill = TRUE)
  keep_methods <- c("grid", "random", "hooke_jeeves", "lbfgsb_multistart", "glmnet_seeded", "fusetune")
  dt <- dt[method %in% keep_methods]
  summ <- dt[, .(mean_gain_per_second = mean(gain_per_second, na.rm = TRUE)), by = .(solver, family, method)]
  summ[, method := factor(method, levels = rev(keep_methods))]
  p <- ggplot(summ, aes(x = method, y = mean_gain_per_second, fill = method == "fusetune")) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#c92a2a", "FALSE" = "#adb5bd")) +
    facet_grid(solver ~ family, scales = "free_x") +
    labs(
      title = "Synthetic family summary by normalized gain per second",
      x = NULL,
      y = "Mean normalized gain per second"
    ) +
    theme_paper()
  ggsave(file.path(fig_dir, "synthetic_family_gain.png"), p, width = 10, height = 5.6, dpi = 300)
}

make_real_efficiency()
make_synth_family()
make_conceptual_figure()
