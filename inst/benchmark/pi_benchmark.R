#!/usr/bin/env Rscript
# Benchmark of prediction interval methods on RepoRT system pairs.
#
# Produces the results shown in the "Prediction interval benchmark" article
# (vignettes/benchmark.qmd): cross-validated coverage, width along the model,
# robustness and runtime, for the package methods and for the alternatives
# that were considered (suggested in B. Hamidova's MSc thesis evaluation).
#
# Usage (from the package root, with the package installed or loaded):
#   Rscript inst/benchmark/pi_benchmark.R [report_path] [out_dir] [n_workers]
#
#   report_path  Extracted RepoRT repository (default RepoRT_data/RepoRT-master)
#   out_dir      Output directory (default output/pi_benchmark, gitignored).
#                Each pair's result is saved in out_dir/pairs as soon as it is
#                done, so an interrupted run resumes where it stopped.
#   n_workers    Parallel workers (default: the allocated cores minus one, from
#                parallelly::availableCores(); pass fewer to keep a desktop usable)
#
# The alternatives need mgcv (gaulss, scat) and qgam. Running everything
# takes about one hour on 7 cores and 3-4 hours on 2, mostly qgam and
# bootstrap. Set the environment variable PI_BENCHMARK_QUICK=1 for a
# 2-pairs-per-group smoke test.

args <- commandArgs(TRUE)
report_path <- if (length(args) >= 1) args[1] else "RepoRT_data/RepoRT-master"
out_dir <- if (length(args) >= 2) args[2] else "output/pi_benchmark"
n_workers <- if (length(args) >= 3) as.integer(args[3]) else parallelly::availableCores(omit = 1)
pair_dir <- file.path(out_dir, "pairs")
dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)
pair_dir <- normalizePath(pair_dir)

suppressMessages(library(rePredRet))
alpha <- 0.05


# ---- Data: RP pairs with >= 10 shared compounds ------------------------------

rd <- suppressMessages(load_report_data(report_path, method_types = "RP"))
lookup <- rePredRet:::prepare_rt_lookup(rd$datasets)
pairs <- rePredRet:::find_model_pairs(lookup, 10)
pairs <- data.frame(id1 = purrr::map_chr(pairs, "id1"), id2 = purrr::map_chr(pairs, "id2"))
pairs <- pairs[pairs$id1 < pairs$id2, ]  # one direction per pair

get_xy <- function(id1, id2) {
  m <- rePredRet:::common_rt_matrix(lookup[[id1]], lookup[[id2]])
  unclass(m)[stats::complete.cases(m), 1:2, drop = FALSE]
}
pairs$n <- purrr::map2_int(pairs$id1, pairs$id2, function(id1, id2) nrow(get_xy(id1, id2)))


# ---- Candidate methods ---------------------------------------------------------
# Each takes the data and a grid and returns list(lower, upper) at the grid.
# The alternatives use the package's monotone fit and only change the interval.

mono_fit_parts <- function(x, y) {
  fit <- rePredRet:::fit_gam_mono(x, y)
  list(
    fit = fit,
    mu = function(xn) rePredRet:::predict_gam_mono(fit, xn),
    se = function(xn) as.numeric(stats::predict(fit$fit_unc, data.frame(x = xn), se.fit = TRUE)$se.fit),
    r = fit$y - rePredRet:::predict_gam_mono(fit, fit$x)
  )
}
from_matrix <- function(o) list(lower = o[, 2], upper = o[, 3])
z <- stats::qnorm(1 - alpha / 2)

# Interval with a global scale s2 (and SE), bounds made monotone like fast_ci
global_scale_pi <- function(b, grid, s2, q = z) {
  from_matrix(rePredRet:::monotone_interval(grid, b$mu(grid), q * sqrt(b$se(grid)^2 + s2)))
}

methods <- list(
  # Package methods
  fast_ci = function(x, y, grid) from_matrix(gam_fast_ci(cbind(x, y), grid, alpha)),
  adaptive_pi = function(x, y, grid) from_matrix(gam_adaptive_pi(cbind(x, y), grid, alpha)),
  bootstrap = function(x, y, grid) {
    bt <- boot::boot(cbind(x, y), rePredRet:::gam_mono_con_fun, R = 200, newdata = c(grid, x))
    from_matrix(rePredRet:::boot2ci_PI(bt, rePredRet:::fit_gam_mono(x, y), grid, alpha))
  },

  # fast_ci before version 2: in-sample residual variance, running maximum
  fast_ci_v1 = function(x, y, grid) {
    b <- mono_fit_parts(x, y)
    mu <- b$mu(grid)
    h <- z * sqrt(b$se(grid)^2 + stats::var(b$r))
    list(lower = cummax(mu - h), upper = cummax(mu + h))
  },

  # bootstrap before version 2: spread of the bootstrap fits at the
  # compounds and along the grid, i.e. a confidence interval of the curve
  bootstrap_v1 = function(x, y, grid) {
    bt <- boot::boot(cbind(x, y), rePredRet:::gam_mono_con_fun, R = 200, newdata = grid)
    sd_grid <- apply(bt$t, 2, stats::sd)
    res_sd <- apply(apply(bt$t, 1, function(t) stats::approx(grid, t, x)$y), 1, stats::sd)
    fit_sd <- ifelse(x %in% grid, 0, stats::approx(grid, sd_grid, x)$y)
    s <- stats::approx(x, sqrt(fit_sd^2 + res_sd^2), grid, ties = function(v) v[1])$y
    list(lower = bt$t0 - z * s, upper = bt$t0 + z * s)
  },

  # Alternatives suggested in the thesis evaluation
  weighted_sigma = function(x, y, grid) {
    b <- mono_fit_parts(x, y)
    global_scale_pi(b, grid, sum(b$fit$w * b$r^2) / sum(b$fit$w))
  },
  mad_sigma = function(x, y, grid) {
    b <- mono_fit_parts(x, y)
    global_scale_pi(b, grid, stats::mad(b$r)^2)
  },
  scat = function(x, y, grid) {
    b <- mono_fit_parts(x, y)
    f <- mgcv::gam(y ~ s(x, k = b$fit$k, bs = "tp"), data = data.frame(x = b$fit$x, y = y),
                   family = mgcv::scat())
    th <- f$family$getTheta(TRUE)
    global_scale_pi(b, grid, th[2]^2, q = stats::qt(1 - alpha / 2, df = th[1]))
  },
  qgam = function(x, y, grid) {
    b <- mono_fit_parts(x, y)
    f <- suppressMessages(qgam::mqgam(y ~ s(x, k = b$fit$k, bs = "tp"),
                                      data = data.frame(x = b$fit$x, y = y),
                                      qu = c(alpha / 2, 1 - alpha / 2),
                                      control = list(progress = FALSE)))
    lo <- as.numeric(qgam::qdo(f, alpha / 2, stats::predict, newdata = data.frame(x = grid)))
    up <- as.numeric(qgam::qdo(f, 1 - alpha / 2, stats::predict, newdata = data.frame(x = grid)))
    mu <- b$mu(grid)
    list(lower = cummax(pmin(lo, mu)), upper = cummax(pmax(up, mu)))
  },
  gaulss = function(x, y, grid) {
    b <- mono_fit_parts(x, y)
    f <- mgcv::gam(list(y ~ s(x, k = b$fit$k, bs = "tp"), ~ s(x, k = max(3, min(b$fit$k, 5)), bs = "tp")),
                   data = data.frame(x = b$fit$x, y = y), family = mgcv::gaulss())
    pr <- stats::predict(f, data.frame(x = grid), type = "response", se.fit = TRUE)
    from_matrix(rePredRet:::monotone_interval(grid, b$mu(grid),
                                              z * sqrt(pr$se.fit[, 1]^2 + 1 / pr$fit[, 2]^2)))
  }
)


# ---- Pair selection ------------------------------------------------------------

set.seed(2026)
quick <- nzchar(Sys.getenv("PI_BENCHMARK_QUICK"))
pick <- function(lo, hi, k) {
  if (quick) k <- min(k, 2)
  cand <- pairs[pairs$n >= lo & pairs$n <= hi, ]
  cand[sample(nrow(cand), min(k, nrow(cand))), ]
}
sel <- rbind(
  transform(pick(10, 15, 40), class = "small (10-15)"),
  transform(pick(30, 60, 40), class = "mid (30-60)"),
  transform(pick(200, Inf, 30), class = "large (200+)")
)

# Pairs with clearly non-constant scatter: largest ratio of the median
# |residual| between x-tertiles, among good fits (R2 >= 0.8)
cand <- pick(30, 400, if (quick) 20 else 600)
cand <- cand[!paste(cand$id1, cand$id2) %in% paste(sel$id1, sel$id2), ]
cand$het <- purrr::map2_dbl(cand$id1, cand$id2, function(id1, id2) {
  m <- get_xy(id1, id2)
  b <- tryCatch(mono_fit_parts(m[, 1], m[, 2]), error = function(e) NULL)
  if (is.null(b) || 1 - sum(b$r^2) / sum((m[, 2] - mean(m[, 2]))^2) < 0.8) return(NA_real_)
  s <- dplyr::tibble(r = abs(b$r), third = cut(rank(b$fit$x, ties.method = "first"), 3, labels = FALSE)) |>
    dplyr::group_by(third) |>
    dplyr::summarize(med = stats::median(r)) |>
    dplyr::pull(med)
  max(s) / max(min(s), 1e-3)
})
het <- cand[order(-cand$het), ][seq_len(if (quick) 2 else 25), ]
het$het <- NULL
sel <- rbind(sel, transform(het, class = "heteroscedastic"))


# ---- Evaluation ----------------------------------------------------------------

# One pair: its index i (seeds the folds), group and RT matrix m
eval_pair <- function(i, class, m) {
  x <- m[, 1]; y <- m[, 2]; n <- length(x)
  set.seed(i)

  # Well-annotated label: sigmoid outlier weight of the full fit >= 0.5.
  # Outliers are taken to be wrong annotations that the PI need not cover.
  b0 <- mono_fit_parts(x, y)
  wa <- b0$fit$w >= 0.5

  # Quiet/noisy label: local median |residual| of the full fit, split at the median
  r0 <- b0$r
  kk <- max(5, round(n / 5))
  loc <- purrr::map_dbl(seq_len(n), function(j) stats::median(abs(r0[order(abs(x - x[j]))[2:(kk + 1)]])))
  noisy <- loc > stats::median(loc)

  K <- if (n <= 30) n else if (n <= 100) 10 else 5
  folds <- if (K == n) seq_len(n) else sample(rep_len(seq_len(K), n))

  purrr::list_rbind(purrr::imap(methods, function(fn, meth) {
    grid <- seq(min(x), max(x), length.out = 1000)
    t0 <- proc.time()[["elapsed"]]
    # Keep the first error message, so that failures are not silent
    error <- NA_character_
    catch <- function(e) {
      if (is.na(error)) error <<- conditionMessage(e)
      NULL
    }
    full <- tryCatch(fn(x, y, grid), error = catch)
    runtime <- proc.time()[["elapsed"]] - t0
    full_ok <- !is.null(full) && all(is.finite(c(full$lower, full$upper)))

    cv_fail <- 0
    cv <- purrr::list_rbind(purrr::map(seq_len(K), function(f) {
      tr <- folds != f
      g <- seq(min(x[tr]), max(x[tr]), length.out = 1000)
      o <- tryCatch(fn(x[tr], y[tr], g), error = catch)
      if (is.null(o) || any(!is.finite(c(o$lower, o$upper)))) {
        cv_fail <<- cv_fail + 1
        return(NULL)
      }
      te <- which(!tr & x >= min(g) & x <= max(g))
      if (length(te) == 0) return(NULL)
      lo <- stats::approx(g, o$lower, x[te])$y
      up <- stats::approx(g, o$upper, x[te])$y
      miss <- pmax(lo - y[te], 0) + pmax(y[te] - up, 0)
      data.frame(covered = y[te] >= lo & y[te] <= up, noisy = noisy[te], wa = wa[te],
                 iscore = (up - lo) + 2 / alpha * miss)
    }))

    if (ncol(cv) == 0) {  # every fold failed
      cv <- data.frame(covered = logical(), noisy = logical(), wa = logical(), iscore = numeric())
    }

    yr <- diff(range(y))
    w <- if (full_ok) full$upper - full$lower else NA
    data.frame(
      method = meth, class = class, n = n, runtime = runtime,
      fail = (!full_ok) + cv_fail,
      n_test = NROW(cv), n_cov = sum(cv$covered),
      n_wa = sum(cv$wa), cov_wa = sum(cv$covered & cv$wa),
      # Quiet/noisy coverage among the well-annotated compounds
      n_quiet = sum(!cv$noisy & cv$wa), cov_quiet = sum(cv$covered & !cv$noisy & cv$wa),
      n_noisy = sum(cv$noisy & cv$wa), cov_noisy = sum(cv$covered & cv$noisy & cv$wa),
      iscore = if (any(cv$wa)) mean(cv$iscore[cv$wa]) / yr else NA,
      med_width = stats::median(w) / yr,
      width_ratio = if (full_ok) stats::quantile(w, 0.95)[[1]] / max(stats::quantile(w, 0.05)[[1]], 1e-9) else NA,
      absurd = full_ok && stats::median(w) > yr,
      error = error
    )
  }))
}

# Saved results of earlier, interrupted runs are reused
run_pair <- function(i, class, m) {
  f <- file.path(pair_dir, sprintf("pair_%03d.rds", i))
  if (file.exists(f)) return(readRDS(f))
  res <- eval_pair(i, class, m)
  saveRDS(res, f)
  res
}

# Only each pair's RT matrix is sent to the workers; one future per pair
# balances the load between small and large pairs
future::plan(future::multisession, workers = n_workers)
per_pair <- purrr::list_rbind(furrr::future_pmap(
  list(i = seq_len(nrow(sel)), class = sel$class, m = purrr::map2(sel$id1, sel$id2, get_xy)),
  run_pair,
  # The method functions are in a list, whose functions furrr does not scan
  # for the objects they use: these are sent to the workers explicitly
  .options = furrr::furrr_options(
    seed = TRUE, scheduling = Inf,
    globals = c("eval_pair", "run_pair", "pair_dir", "methods", "mono_fit_parts",
                "from_matrix", "global_scale_pi", "alpha", "z"),
    packages = "rePredRet"
  )
))
future::plan(future::sequential)
utils::write.csv(per_pair, file.path(out_dir, "per_pair.csv"), row.names = FALSE)

# A method that fails on every fit points to a broken setup, not a weak method
broken <- per_pair |>
  dplyr::group_by(method) |>
  dplyr::summarize(all_failed = all(n_test == 0), error = dplyr::first(stats::na.omit(error)))
if (any(broken$all_failed)) {
  stop("Every fit failed for: ",
       paste0(broken$method[broken$all_failed], " (", broken$error[broken$all_failed], ")",
              collapse = "; "))
}

labels <- c(fast_ci_v1 = "fast_ci (before fix)", bootstrap_v1 = "bootstrap (before fix)",
            weighted_sigma = "weighted sigma", mad_sigma = "MAD sigma")
per_pair$method <- ifelse(per_pair$method %in% names(labels), labels[per_pair$method], per_pair$method)

per_pair_all <- rbind(per_pair, transform(per_pair, class = "all"))
summ <- per_pair_all |>
  dplyr::group_by(class, method) |>
  dplyr::summarize(
    pairs = dplyr::n(),
    coverage = sum(cov_wa) / sum(n_wa),
    cov_all = sum(n_cov) / sum(n_test),
    cov_outliers = (sum(n_cov) - sum(cov_wa)) / (sum(n_test) - sum(n_wa)),
    outlier_share = 1 - sum(n_wa) / sum(n_test),
    cov_quiet = sum(cov_quiet) / sum(n_quiet),
    cov_noisy = sum(cov_noisy) / sum(n_noisy),
    med_iscore = stats::median(iscore, na.rm = TRUE),
    med_width_rel = stats::median(med_width, na.rm = TRUE),
    width_ratio_q95q05 = stats::median(width_ratio, na.rm = TRUE),
    fail_fits = sum(fail), absurd = sum(absurd),
    .groups = "drop"
  ) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.double), \(v) round(v, 3)))
utils::write.csv(summ, file.path(out_dir, "benchmark_summary.csv"), row.names = FALSE)
message("Wrote ", file.path(out_dir, "benchmark_summary.csv"))
