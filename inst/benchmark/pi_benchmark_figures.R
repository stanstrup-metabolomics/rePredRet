#!/usr/bin/env Rscript
# Example figures of the "Prediction interval benchmark" article
# (vignettes/benchmark.qmd): the bands of the interval methods on pairs with
# unequal scatter and on pairs of different sizes, and the width along the
# model.
#
# Usage (from the package root, with the package installed):
#   Rscript inst/benchmark/pi_benchmark_figures.R [report_path] [out_dir]
#
#   report_path  Extracted RepoRT repository (default RepoRT_data/RepoRT-master)
#   out_dir      Output directory (default vignettes/benchmark)

args <- commandArgs(TRUE)
report_path <- if (length(args) >= 1) args[1] else "RepoRT_data/RepoRT-master"
out_dir <- if (length(args) >= 2) args[2] else "vignettes/benchmark"

suppressMessages(library(rePredRet))
library(ggplot2)
alpha <- 0.05
z <- stats::qnorm(1 - alpha / 2)

rd <- suppressMessages(load_report_data(report_path, method_types = "RP"))
lookup <- rePredRet:::prepare_rt_lookup(rd$datasets)
get_xy <- function(id1, id2) {
  m <- rePredRet:::common_rt_matrix(lookup[[id1]], lookup[[id2]])
  unclass(m)[stats::complete.cases(m), 1:2, drop = FALSE]
}

# Pairs shown in the article
heteroscedastic <- data.frame(
  id1 = c("0345", "0311", "0125", "0068", "0243", "0125"),
  id2 = c("0433", "0433", "0195", "0189", "0433", "0155"),
  class = "heteroscedastic"
)
sizes <- data.frame(
  id1 = c("0286", "0050", "0198", "0095", "0322", "0249"),
  id2 = c("0320", "0259", "0213", "0101", "0336", "0254"),
  class = rep(c("small", "mid", "large"), each = 2)
)

methods <- list(
  "fast_ci (old)" = function(x, y, grid) {
    fit <- rePredRet:::fit_gam_mono(x, y)
    mu <- rePredRet:::predict_gam_mono(fit, grid)
    se <- as.numeric(stats::predict(fit$fit_unc, data.frame(x = grid), se.fit = TRUE)$se.fit)
    r <- fit$y - rePredRet:::predict_gam_mono(fit, fit$x)
    h <- z * sqrt(se^2 + stats::var(r))
    cbind(mu, cummax(mu - h), cummax(mu + h))
  },
  fast_ci = function(x, y, grid) gam_fast_ci(cbind(x, y), grid, alpha),
  adaptive_pi = function(x, y, grid) gam_adaptive_pi(cbind(x, y), grid, alpha),
  bootstrap = function(x, y, grid) {
    set.seed(1)
    bt <- boot::boot(cbind(x, y), rePredRet:::gam_mono_con_fun, R = 200, newdata = c(grid, x))
    rePredRet:::boot2ci_PI(bt, rePredRet:::fit_gam_mono(x, y), grid, alpha)
  }
)
colours <- c("fast_ci (old)" = "grey60", fast_ci = "#2c7bdb", adaptive_pi = "#f0662f", bootstrap = "#1aa37a")
line_types <- c("fast_ci (old)" = "dashed", fast_ci = "solid", adaptive_pi = "solid", bootstrap = "solid")

band_data <- function(pairs) {
  purrr::pmap(pairs, function(id1, id2, class) {
    m <- get_xy(id1, id2)
    x <- m[, 1]; y <- m[, 2]
    grid <- seq(min(x), max(x), length.out = 300)
    panel <- sprintf("%s -> %s (n = %d, %s)", id1, id2, length(x), class)
    bands <- purrr::imap(methods, function(fn, meth) {
      o <- fn(x, y, grid)
      data.frame(panel = panel, method = meth, x = grid, fit = o[, 1], lower = o[, 2], upper = o[, 3])
    })
    list(points = data.frame(panel = panel, x = x, y = y), bands = purrr::list_rbind(bands))
  })
}

order_panels <- function(d, panels) transform(d, panel = factor(panel, levels = panels),
                                              method = factor(method, levels = names(methods)))

plot_bands <- function(dat) {
  panels <- unique(purrr::map_chr(dat, function(p) p$points$panel[1]))
  points <- transform(purrr::list_rbind(purrr::map(dat, "points")), panel = factor(panel, levels = panels))
  bands <- order_panels(purrr::list_rbind(purrr::map(dat, "bands")), panels)
  fit <- bands[bands$method == "adaptive_pi", ]
  long <- rbind(transform(bands, bound = "lower", value = lower),
                transform(bands, bound = "upper", value = upper))
  ggplot() +
    geom_point(data = points, aes(x, y), colour = "grey40", alpha = 0.7, size = 1.5) +
    geom_line(data = fit, aes(x, fit), linetype = "dashed", colour = "black") +
    geom_line(data = long, aes(x, value, colour = method, linetype = method,
                               group = interaction(method, bound)), linewidth = 0.8) +
    facet_wrap(~panel, ncol = 2, scales = "free") +
    scale_colour_manual(values = colours, name = "95% PI bounds") +
    scale_linetype_manual(values = line_types, name = "95% PI bounds") +
    labs(x = "RT in source system", y = "RT in target system",
         caption = "Points: calibration compounds. Dashed: monotone GAM fit (shared by all methods).") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
}

het <- band_data(heteroscedastic)
ggsave(file.path(out_dir, "fig_examples_heteroscedastic.png"), plot_bands(het),
       width = 11, height = 11, dpi = 150, bg = "white")
ggsave(file.path(out_dir, "fig_examples_sizes.png"), plot_bands(band_data(sizes)),
       width = 11, height = 11, dpi = 150, bg = "white")

panels <- unique(purrr::map_chr(het, function(p) p$points$panel[1]))
widths <- order_panels(purrr::list_rbind(purrr::map(het, "bands")), panels)
p_width <- ggplot(widths, aes(x, upper - lower, colour = method, linetype = method)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~panel, ncol = 2, scales = "free") +
  scale_colour_manual(values = colours, name = NULL) +
  scale_linetype_manual(values = line_types, name = NULL) +
  labs(x = "RT in source system", y = "PI width", title = "Interval width along the model") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "fig_width_along_x.png"), p_width,
       width = 11, height = 10, dpi = 150, bg = "white")
message("Wrote figures to ", out_dir)
