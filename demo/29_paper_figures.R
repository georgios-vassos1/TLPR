## ============================================================
## 29_paper_figures.R — Publication figures for the MSTP paper
##
## Loads pre-computed RDS results from demo/results/ and
## generates all figures for the computational study:
##
##   F1 — Corrmat heatmap (2-block 6×6 Gaussian copula)
##   F2 — Sensitivity 3-panel (UB vs λ, ρ_cross, spot_mult)
##   F3 — Regret density/violin across topologies (demo/21)
##   F4 — Regret% vs SDDP iterations (demo/19)
##   F5 — Capacity optimisation convergence (demo/28 Part B)
##
## Run all demos first (they save RDS to demo/results/).
## Outputs go to demo/figures/ as PDF (7×5 inches).
## ============================================================

suppressPackageStartupMessages({
  library(MSTP)
  library(ggplot2)
  library(data.table)
})

## ── Output directory ──────────────────────────────────────────────────────────

FIGURES_DIR <- "demo/figures"
RESULTS_DIR <- "demo/results"
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

## ── Publication theme ─────────────────────────────────────────────────────────

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor    = element_blank(),
      strip.background    = element_rect(fill = "grey92", colour = "grey70"),
      legend.position     = "bottom",
      legend.key.width    = unit(1.5, "cm"),
      plot.title          = element_text(size = base_size, face = "bold"),
      axis.title          = element_text(size = base_size - 1),
      axis.text           = element_text(size = base_size - 2)
    )
}

save_fig <- function(p, name, width = 7, height = 5) {
  path <- file.path(FIGURES_DIR, paste0(name, ".pdf"))
  ggsave(path, p, width = width, height = height)
  cat(sprintf("  Saved: %s\n", path))
  invisible(path)
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("29_paper_figures.R — Generating publication figures\n")
cat(strrep("=", 70), "\n\n")

## ══════════════════════════════════════════════════════════════════════════════
## F1 — Corrmat heatmap
## ══════════════════════════════════════════════════════════════════════════════

cat("F1: Corrmat heatmap...\n")

set.seed(3L)  # seed 3 -> rho_OO=0.6, rho_DD=0.7 for reproducible caption
cm <- MSTP::gen_corrmat(n_blocks = 2L, block_size = 6L, cross_corr = 0.0)
n  <- nrow(cm)

cm_dt <- data.table(
  row = rep(seq_len(n), each = n),
  col = rep(seq_len(n), times = n),
  val = as.vector(cm)
)
cm_dt[, row_rev := (n + 1L) - row]

f1 <- ggplot(cm_dt, aes(x = col, y = row_rev, fill = val)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white", high = "#b2182b",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    name = "Correlation"
  ) +
  scale_x_continuous(breaks = c(1, 6, 7, 12),
                     labels = c("O1", "O6", "D1", "D6"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = c(n + 1L - c(1L, 6L, 7L, 12L)),
                     labels = c("O1", "O6", "D1", "D6"),
                     expand = c(0, 0)) +
  geom_vline(xintercept = 6.5, colour = "grey40", linewidth = 0.6, linetype = "dashed") +
  geom_hline(yintercept = 6.5, colour = "grey40", linewidth = 0.6, linetype = "dashed") +
  labs(
    title = expression("Gaussian copula corrmat: " * rho[OO] * " = 0.6, " * rho[DD] * " = 0.7, " * rho[OD] * " = 0"),
    x = NULL, y = NULL
  ) +
  coord_fixed() +
  theme_pub()

save_fig(f1, "F1_corrmat", width = 5.5, height = 5)

## ══════════════════════════════════════════════════════════════════════════════
## F2 — Sensitivity 3-panel
## ══════════════════════════════════════════════════════════════════════════════

cat("F2: Sensitivity 3-panel...\n")

rds_25 <- file.path(RESULTS_DIR, "25_results.rds")
if (!file.exists(rds_25)) {
  cat("  SKIP: demo/results/25_results.rds not found (run demo/25 first)\n")
} else {
  s25 <- readRDS(rds_25)

  ## Panel A: UB vs lambda
  dt_lam <- rbindlist(lapply(s25$lambda, function(r)
    data.table(x = r$lambda, ub = r$ub, ub_sd = r$ub_sd, lb = r$lb,
               gap = r$gap_pct)))

  pa <- ggplot(dt_lam, aes(x = factor(x), y = ub)) +
    geom_errorbar(aes(ymin = ub - ub_sd, ymax = ub + ub_sd),
                  width = 0.2, colour = "steelblue") +
    geom_point(size = 3, colour = "steelblue4") +
    labs(x = expression(lambda ~ "(demand rate)"),
         y = "UB (mean ± sd)", title = "(a) Demand rate") +
    theme_pub()

  ## Panel B: UB vs rho_cross
  dt_cor <- rbindlist(lapply(s25$corr, function(r)
    data.table(x = r$cross_corr, ub = r$ub, ub_sd = r$ub_sd, lb = r$lb,
               gap = r$gap_pct)))

  pb <- ggplot(dt_cor, aes(x = factor(x), y = ub)) +
    geom_errorbar(aes(ymin = ub - ub_sd, ymax = ub + ub_sd),
                  width = 0.2, colour = "darkorange") +
    geom_point(size = 3, colour = "darkorange4") +
    labs(x = expression(rho[cross] ~ "(copula correlation)"),
         y = "UB (mean ± sd)", title = "(b) Cross-block correlation") +
    theme_pub()

  ## Panel C: UB vs spot multiplier
  dt_spt <- rbindlist(lapply(s25$spot, function(r)
    data.table(x = r$spot_mult, ub = r$ub, ub_sd = r$ub_sd, lb = r$lb,
               gap = r$gap_pct)))

  pc <- ggplot(dt_spt, aes(x = factor(x), y = ub)) +
    geom_errorbar(aes(ymin = ub - ub_sd, ymax = ub + ub_sd),
                  width = 0.2, colour = "forestgreen") +
    geom_point(size = 3, colour = "darkgreen") +
    labs(x = expression(italic(m) ~ "(spot rate multiplier)"),
         y = "UB (mean ± sd)", title = "(c) Spot rate multiplier") +
    theme_pub()

  ## Combine with patchwork if available, otherwise save separately
  if (requireNamespace("patchwork", quietly = TRUE)) {
    f2 <- patchwork::wrap_plots(pa, pb, pc, nrow = 1L) +
      patchwork::plot_annotation(
        title = "Sensitivity analysis: SDDP upper bound vs key parameters",
        subtitle = "6x6x20 instance, tau=12, 500 SDDP iterations, 200 OOB trials"
      )
    save_fig(f2, "F2_sensitivity", width = 12, height = 4.5)
  } else {
    save_fig(pa, "F2a_sensitivity_lambda",   width = 4, height = 4)
    save_fig(pb, "F2b_sensitivity_corr",     width = 4, height = 4)
    save_fig(pc, "F2c_sensitivity_spot",     width = 4, height = 4)
    cat("  Tip: install patchwork for a combined 3-panel figure\n")
  }
}

## ══════════════════════════════════════════════════════════════════════════════
## F3 — Regret density/violin across topologies (demo/21)
## ══════════════════════════════════════════════════════════════════════════════

cat("F3: Regret density across topologies...\n")

rds_21 <- file.path(RESULTS_DIR, "21_results.rds")
if (!file.exists(rds_21)) {
  cat("  SKIP: demo/results/21_results.rds not found (run demo/21 first)\n")
} else {
  r21 <- readRDS(rds_21)

  dt_reg <- rbindlist(lapply(r21, function(res)
    data.table(topology = res$label,
               regret   = res$regret * 100)))

  ## Order topologies by increasing state space
  topo_levels <- sapply(r21, `[[`, "label")
  dt_reg[, topology := factor(topology, levels = topo_levels)]

  f3 <- ggplot(dt_reg, aes(x = topology, y = regret, fill = topology)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7, colour = "grey30") +
    geom_boxplot(width = 0.12, outlier.size = 0.5,
                 fill = "white", colour = "grey20") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title    = "SDDP regret distribution across network topologies",
      subtitle = "Regret = (SDDP cost − clairvoyant LP) / |clairvoyant LP|",
      x        = "Network topology",
      y        = "Regret (%)"
    ) +
    theme_pub()

  save_fig(f3, "F3_regret_topology", width = 7, height = 5)
}

## ══════════════════════════════════════════════════════════════════════════════
## F4 — Regret% vs SDDP iterations (demo/19)
## ══════════════════════════════════════════════════════════════════════════════

cat("F4: Regret vs SDDP iterations...\n")

rds_19 <- file.path(RESULTS_DIR, "19_results.rds")
if (!file.exists(rds_19)) {
  cat("  SKIP: demo/results/19_results.rds not found (run demo/19 first)\n")
} else {
  r19 <- readRDS(rds_19)

  ## Summary stats per iteration count
  dt_it <- rbindlist(lapply(r19, function(res) {
    rp <- res$regret * 100
    data.table(
      iters      = res$iters,
      gap_pct    = res$gap_pct,
      reg_mean   = mean(rp),
      reg_median = median(rp),
      reg_lo     = quantile(rp, 0.25),
      reg_hi     = quantile(rp, 0.75)
    )
  }))

  f4 <- ggplot(dt_it, aes(x = iters)) +
    geom_ribbon(aes(ymin = reg_lo, ymax = reg_hi),
                fill = "steelblue", alpha = 0.25) +
    geom_line(aes(y = reg_mean,   colour = "Mean"),   linewidth = 1) +
    geom_line(aes(y = reg_median, colour = "Median"), linewidth = 1, linetype = "dashed") +
    geom_point(aes(y = reg_mean,   colour = "Mean"),   size = 3) +
    geom_point(aes(y = reg_median, colour = "Median"), size = 3) +
    scale_colour_manual(values = c(Mean = "steelblue4", Median = "tomato3"),
                        name = NULL) +
    scale_x_continuous(breaks = dt_it$iters) +
    labs(
      title    = "SDDP regret vs training iterations",
      subtitle = "2x2 instance, R=10, tau=4 | Ribbon = IQR | 500 OOB trials",
      x        = "SDDP iterations",
      y        = "Regret (%)"
    ) +
    theme_pub()

  save_fig(f4, "F4_regret_iterations", width = 7, height = 5)
}

## (F5 — capacity optimisation convergence removed: demo/31 archived)

## ══════════════════════════════════════════════════════════════════════════════
## Done
## ══════════════════════════════════════════════════════════════════════════════

cat(sprintf("\nAll figures written to %s/\n", FIGURES_DIR))
cat("Figure inventory:\n")
figs <- list.files(FIGURES_DIR, pattern = "\\.pdf$", full.names = FALSE)
for (f in figs) cat(sprintf("  %s\n", f))

invisible(NULL)
