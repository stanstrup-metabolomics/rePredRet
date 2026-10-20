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
#   out_dir      Output directory (default benchmark_output)
#   n_workers    Parallel workers (default: available cores - 1)
#
# The alternatives need mgcv (gaulss, scat) and qgam. Running everything
# takes about one hour on 7 cores, mostly qgam and bootstrap. Set the
# environment variable PI_BENCHMARK_QUICK=1 for a 2-pairs-per-group smoke test.

args <- commandArgs(TRUE)
report_path <- if (length(args) >= 1) args[1] else "RepoRT_data/RepoRT-master"
out_dir <- if (length(args) >= 2) args[2] else "benchmark_output"
n_workers <- if (length(args) >= 3) as.integer(args[3]) else max(1, parallelly::availableCores() - 1)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(library(rePredRet))
alpha <- 0.05


# ---- Data: RP pairs with >= 10 shared compounds ------------------------------

rd <- suppressMessages(load_report_data(report_path, method_types = "RP"))
lookup <- rePredRet:::prepare_rt_lookup(rd$datasets)
pairs <- rePredRet:::find_model_pairs(lookup, 10)
pairs <- data.frame(id1 = vapply(pairs, `[[`, "", "id1"), id2 = vapply(pairs, `[[`, "", "id2"))
pairs <- pairs[pairs$id1 < pairs$id2, ]  # one direction per pair

get_xy <- function(id1, id2) {
  m <- rePredRet:::common_rt_matrix(lookup[[id1]], lookup[[id2]])
  unclass(m)[stats::complete.cases(m), 1:2, drop = FALSE]
}
pairs$n <- vapply(seq_len(nrow(pairs)), function(i) nrow(get_xy(pairs$id1[i], pairs$id2[i])), numeric(1))


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
cand$het <- vapply(seq_len(nrow(cand)), function(i) {
  m <- get_xy(cand$id1[i], cand$id2[i])
  b <- tryCatch(mono_fit_parts(m[, 1], m[, 2]), error = function(e) NULL)
  if (is.null(b) || 1 - sum(b$r^2) / sum((m[, 2] - mean(m[, 2]))^2) < 0.8) return(NA_real_)
  s <- tapply(abs(b$r), cut(rank(b$fit$x, ties.method = "first"), 3, labels = FALSE), stats::median)
  max(s) / max(min(s), 1e-3)
}, numeric(1))
het <- cand[order(-cand$het), ][seq_len(if (quick) 2 else 25), ]
het$het <- NULL
sel <- rbind(sel, transform(het, class = "heteroscedastic"))


# ---- Evaluation ----------------------------------------------------------------

eval_pair <- function(i) {
  m <- get_xy(sel$id1[i], sel$id2[i])
  x <- m[, 1]; y <- m[, 2]; n <- length(x)
  set.seed(i)

  # Quiet/noisy label: local median |residual| of the full fit, split at the median
  r0 <- mono_fit_parts(x, y)$r
  kk <- max(5, round(n / 5))
  loc <- vapply(seq_len(n), function(j) stats::median(abs(r0[order(abs(x - x[j]))[2:(kk + 1)]])), numeric(1))
  noisy <- loc > stats::median(loc)

  K <- if (n <= 30) n else if (n <= 100) 10 else 5
  folds <- if (K == n) seq_len(n) else sample(rep_len(seq_len(K), n))

  do.call(rbind, lapply(names(methods), function(meth) {
    fn <- methods[[meth]]
    grid <- seq(min(x), max(x), length.out = 1000)
    t0 <- proc.time()[["elapsed"]]
    full <- tryCatch(fn(x, y, grid), error = function(e) NULL)
    runtime <- proc.time()[["elapsed"]] - t0
    full_ok <- !is.null(full) && all(is.finite(c(full$lower, full$upper)))

    cv_fail <- 0
    cv <- do.call(rbind, lapply(seq_len(K), function(f) {
      tr <- folds != f
      g <- seq(min(x[tr]), max(x[tr]), length.out = 1000)
      o <- tryCatch(fn(x[tr], y[tr], g), error = function(e) NULL)
      if (is.null(o) || any(!is.finite(c(o$lower, o$upper)))) {
        cv_fail <<- cv_fail + 1
        return(NULL)
      }
      te <- which(!tr & x >= min(g) & x <= max(g))
      if (length(te) == 0) return(NULL)
      lo <- stats::approx(g, o$lower, x[te])$y
      up <- stats::approx(g, o$upper, x[te])$y
      miss <- pmax(lo - y[te], 0) + pmax(y[te] - up, 0)
      data.frame(covered = y[te] >= lo & y[te] <= up, noisy = noisy[te],
                 iscore = (up - lo) + 2 / alpha * miss)
    }))

    yr <- diff(range(y))
    w <- if (full_ok) full$upper - full$lower else NA
    data.frame(
      method = meth, class = sel$class[i], n = n, runtime = runtime,
      fail = (!full_ok) + cv_fail,
      n_test = NROW(cv), n_cov = sum(cv$covered),
      n_quiet = sum(!cv$noisy), cov_quiet = sum(cv$covered & !cv$noisy),
      n_noisy = sum(cv$noisy), cov_noisy = sum(cv$covered & cv$noisy),
      iscore = if (is.null(cv)) NA else mean(cv$iscore) / yr,
      med_width = stats::median(w) / yr,
      width_ratio = if (full_ok) stats::quantile(w, 0.95)[[1]] / max(stats::quantile(w, 0.05)[[1]], 1e-9) else NA,
      absurd = full_ok && stats::median(w) > yr
    )
  }))
}

cl <- parallel::makeCluster(n_workers)
parallel::clusterEvalQ(cl, suppressMessages(library(rePredRet)))
parallel::clusterExport(cl, setdiff(ls(), "cl"))
per_pair <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(sel)), eval_pair))
parallel::stopCluster(cl)
utils::write.csv(per_pair, file.path(out_dir, "per_pair.csv"), row.names = FALSE)

per_pair_all <- rbind(per_pair, transform(per_pair, class = "all"))
summ <- do.call(rbind, lapply(split(per_pair_all, list(per_pair_all$method, per_pair_all$class), drop = TRUE), function(d) {
  data.frame(
    method = d$method[1], class = d$class[1], pairs = nrow(d),
    coverage = sum(d$n_cov) / sum(d$n_test),
    cov_quiet = sum(d$cov_quiet) / sum(d$n_quiet),
    cov_noisy = sum(d$cov_noisy) / sum(d$n_noisy),
    med_iscore = stats::median(d$iscore, na.rm = TRUE),
    med_width_rel = stats::median(d$med_width, na.rm = TRUE),
    width_ratio = stats::median(d$width_ratio, na.rm = TRUE),
    fail_fits = sum(d$fail), absurd = sum(d$absurd)
  )
}))
utils::write.csv(summ, file.path(out_dir, "benchmark_summary.csv"), row.names = FALSE)
message("Wrote ", file.path(out_dir, "benchmark_summary.csv"))
