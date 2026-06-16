#!/usr/bin/env Rscript
# Tests for the R MAPA implementation.
# Run with: Rscript test_mapa.R  (from the reference/r/ directory, or anywhere)

library(testthat)

# Source mapa.R relative to this script's location, whether run via Rscript
# or sourced interactively.
this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(0)$ofile)),  # when source()d interactively
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      getwd()
    }
  }
)

source(file.path(this_dir, "mapa.R"))

FIXTURES_DIR   <- normalizePath(file.path(this_dir, "..", "fixtures"))
BAYESIAN_K     <- 10
MIN_OBS        <- 50
MIN_BADS       <- 10
MIN_CONFIDENCE <- 0.95

# ---------------------------------------------------------------------------
# Fixture loaders
# ---------------------------------------------------------------------------

load_observations <- function() {
  read.csv(file.path(FIXTURES_DIR, "raw_observations.csv"), stringsAsFactors = FALSE)
}

load_bins <- function(filename) {
  df <- read.csv(file.path(FIXTURES_DIR, filename), stringsAsFactors = FALSE)
  # Ensure integer columns match R's numeric representation
  df$n_obs  <- as.integer(df$n_obs)
  df$n_bads <- as.integer(df$n_bads)
  df
}

load_calibrated_bins <- function(filename) {
  df <- load_bins(filename)
  df
}

bins_equal <- function(a, b) {
  if (nrow(a) != nrow(b)) return(FALSE)
  all(a$score_min == b$score_min) &&
    all(a$score_max == b$score_max) &&
    all(a$n_obs    == b$n_obs) &&
    all(a$n_bads   == b$n_bads)
}

calibrated_bins_equal <- function(a, b, tol = 1e-9) {
  if (!bins_equal(a, b)) return(FALSE)
  all(abs(a$pd - b$pd) < tol | (abs(a$pd - b$pd) / pmax(abs(b$pd), 1e-15)) < tol)
}

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

test_that("bins_from_observations matches expected", {
  obs      <- load_observations()
  result   <- bins_from_observations(obs)
  expected <- load_bins("expected_initial_bins.csv")

  expect_true(bins_equal(result, expected))
})

test_that("calibrate matches expected pooled bins", {
  obs      <- load_observations()
  result   <- calibrate(obs)
  expected <- load_bins("expected_pooled_bins.csv")

  expect_true(bins_equal(result, expected))
})

test_that("result is monotone non-increasing", {
  obs    <- load_observations()
  result <- calibrate(obs)
  rates  <- result$n_bads / result$n_obs

  expect_equal(rates, sort(rates, decreasing = TRUE))
})

test_that("pooling preserves totals", {
  obs     <- load_observations()
  initial <- bins_from_observations(obs)
  result  <- mapa(initial)

  expect_equal(sum(result$n_obs),  sum(initial$n_obs))
  expect_equal(sum(result$n_bads), sum(initial$n_bads))
  expect_equal(sum(result$n_obs),  nrow(obs))
})

test_that("enforce_minimum_size matches expected", {
  obs    <- load_observations()
  pooled <- calibrate(obs)
  result <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)

  expected <- load_bins("expected_min_size_bins.csv")
  expect_true(bins_equal(result, expected))
})

test_that("enforce_minimum_size satisfies thresholds", {
  obs    <- load_observations()
  pooled <- calibrate(obs)
  result <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)

  if (nrow(result) > 1) {
    expect_true(all(result$n_obs  >= MIN_OBS))
    expect_true(all(result$n_bads >= MIN_BADS))
  }
})

test_that("enforce_minimum_size preserves totals and monotonicity", {
  obs    <- load_observations()
  pooled <- calibrate(obs)
  result <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)

  expect_equal(sum(result$n_obs),  sum(pooled$n_obs))
  expect_equal(sum(result$n_bads), sum(pooled$n_bads))

  rates <- result$n_bads / result$n_obs
  expect_equal(rates, sort(rates, decreasing = TRUE))
})

test_that("enforce_minimum_size is noop with default thresholds", {
  obs    <- load_observations()
  pooled <- calibrate(obs)
  result <- enforce_minimum_size(pooled)

  expect_true(bins_equal(result, pooled))
})

test_that("Bayesian adjustment matches expected", {
  obs      <- load_observations()
  pooled   <- calibrate(obs)
  sized    <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)
  result   <- apply_bayesian_adjustment(sized, k = BAYESIAN_K)
  expected <- load_calibrated_bins("expected_calibrated_bins.csv")

  expect_equal(nrow(result), nrow(expected))
  expect_equal(result$score_min, expected$score_min)
  expect_equal(result$score_max, expected$score_max)
  expect_equal(result$n_obs,     expected$n_obs)
  expect_equal(result$n_bads,    expected$n_bads)
  expect_equal(result$pd, expected$pd, tolerance = 1e-9)
})

test_that("Bayesian adjustment shrinks toward prior", {
  obs    <- load_observations()
  pooled <- calibrate(obs)
  prior  <- sum(pooled$n_bads) / sum(pooled$n_obs)
  result <- apply_bayesian_adjustment(pooled, k = BAYESIAN_K, prior = prior)

  for (i in seq_len(nrow(pooled))) {
    orig_rate <- pooled$n_bads[i] / pooled$n_obs[i]
    adj_pd    <- result$pd[i]
    lo <- min(orig_rate, prior)
    hi <- max(orig_rate, prior)
    expect_true(adj_pd >= lo - .Machine$double.eps * 100)
    expect_true(adj_pd <= hi + .Machine$double.eps * 100)
  }
})

test_that("repool_calibrated_bins matches expected", {
  obs        <- load_observations()
  pooled     <- calibrate(obs)
  sized      <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)
  calibrated <- apply_bayesian_adjustment(sized, k = BAYESIAN_K)
  result     <- repool_calibrated_bins(calibrated)
  expected   <- load_calibrated_bins("expected_repooled_calibrated_bins.csv")

  expect_equal(nrow(result), nrow(expected))
  expect_equal(result$score_min, expected$score_min)
  expect_equal(result$score_max, expected$score_max)
  expect_equal(result$n_obs,     expected$n_obs)
  expect_equal(result$n_bads,    expected$n_bads)
  expect_equal(result$pd, expected$pd, tolerance = 1e-9)
})

test_that("repool_calibrated_bins restores monotonicity", {
  obs        <- load_observations()
  pooled     <- calibrate(obs)
  sized      <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)
  calibrated <- apply_bayesian_adjustment(sized, k = BAYESIAN_K)

  # Fixture deliberately has a non-monotone pd sequence after Bayesian adjustment
  pds_before <- calibrated$pd
  expect_false(identical(pds_before, sort(pds_before, decreasing = TRUE)))

  result  <- repool_calibrated_bins(calibrated)
  pds_after <- result$pd
  expect_equal(pds_after, sort(pds_after, decreasing = TRUE))
})

test_that("repool_calibrated_bins preserves totals", {
  obs        <- load_observations()
  pooled     <- calibrate(obs)
  sized      <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)
  calibrated <- apply_bayesian_adjustment(sized, k = BAYESIAN_K)
  result     <- repool_calibrated_bins(calibrated)

  expect_equal(sum(result$n_obs),  sum(calibrated$n_obs))
  expect_equal(sum(result$n_bads), sum(calibrated$n_bads))
})

test_that("interpolate_pd matches expected smoothed PDs", {
  obs        <- load_observations()
  pooled     <- calibrate(obs)
  sized      <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)
  calibrated <- apply_bayesian_adjustment(sized, k = BAYESIAN_K)
  repooled   <- repool_calibrated_bins(calibrated)

  expected <- read.csv(file.path(FIXTURES_DIR, "expected_smoothed_pds.csv"),
                       stringsAsFactors = FALSE)

  for (i in seq_len(nrow(expected))) {
    score    <- expected$score[i]
    exp_pd   <- expected$pd[i]
    result   <- interpolate_pd(repooled, score)
    # relative tolerance 1e-9
    expect_equal(result, exp_pd, tolerance = 1e-9)
  }
})

test_that("interpolate_pd is monotone non-increasing", {
  obs        <- load_observations()
  pooled     <- calibrate(obs)
  sized      <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)
  calibrated <- apply_bayesian_adjustment(sized, k = BAYESIAN_K)
  repooled   <- repool_calibrated_bins(calibrated)

  scores <- sort(unique(c(repooled$score_min, repooled$score_max)))
  pds    <- sapply(scores, function(s) interpolate_pd(repooled, s))
  expect_equal(pds, sort(pds, decreasing = TRUE))
})

test_that("interpolate_pd at pool midpoint matches pool pd", {
  obs        <- load_observations()
  pooled     <- calibrate(obs)
  sized      <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)
  calibrated <- apply_bayesian_adjustment(sized, k = BAYESIAN_K)
  repooled   <- repool_calibrated_bins(calibrated)

  for (i in seq_len(nrow(repooled))) {
    midpoint <- (repooled$score_min[i] + repooled$score_max[i]) / 2
    result   <- interpolate_pd(repooled, midpoint)
    expect_equal(result, repooled$pd[i], tolerance = 1e-9)
  }
})

test_that("run_pipeline bands match repool_calibrated_bins", {
  obs        <- load_observations()
  pooled     <- calibrate(obs)
  sized      <- enforce_minimum_size(pooled, min_obs = MIN_OBS, min_bads = MIN_BADS)
  calibrated <- apply_bayesian_adjustment(sized, k = BAYESIAN_K)
  repooled   <- repool_calibrated_bins(calibrated)

  pipeline <- run_pipeline(obs, k = BAYESIAN_K, min_obs = MIN_OBS, min_bads = MIN_BADS)

  expect_true(calibrated_bins_equal(pipeline$bands, repooled))
})

test_that("run_pipeline pd_for_score matches interpolate_pd", {
  obs      <- load_observations()
  pipeline <- run_pipeline(obs, k = BAYESIAN_K, min_obs = MIN_OBS, min_bads = MIN_BADS)

  expected <- read.csv(file.path(FIXTURES_DIR, "expected_smoothed_pds.csv"),
                       stringsAsFactors = FALSE)

  for (i in seq_len(nrow(expected))) {
    score  <- expected$score[i]
    exp_pd <- expected$pd[i]
    result <- pipeline$pd_for_score(score)
    expect_equal(result, interpolate_pd(pipeline$bands, score), tolerance = 1e-9)
    expect_equal(result, exp_pd, tolerance = 1e-9)
  }
})

test_that("mapa min_confidence matches expected", {
  obs     <- load_observations()
  initial <- bins_from_observations(obs)
  result  <- mapa(initial, min_confidence = MIN_CONFIDENCE)

  expected <- load_bins("expected_pooled_bins_confidence.csv")
  expect_true(bins_equal(result, expected))
})

test_that("mapa min_confidence preserves totals and monotonicity", {
  obs     <- load_observations()
  initial <- bins_from_observations(obs)
  result  <- mapa(initial, min_confidence = MIN_CONFIDENCE)

  expect_equal(sum(result$n_obs),  sum(initial$n_obs))
  expect_equal(sum(result$n_bads), sum(initial$n_bads))

  rates <- result$n_bads / result$n_obs
  expect_equal(rates, sort(rates, decreasing = TRUE))
})

test_that("mapa min_confidence merges at least as much as plain mapa", {
  obs     <- load_observations()
  initial <- bins_from_observations(obs)
  plain     <- mapa(initial)
  confident <- mapa(initial, min_confidence = MIN_CONFIDENCE)

  expect_lte(nrow(confident), nrow(plain))
})

test_that("increasing direction works", {
  obs            <- load_observations()
  flipped        <- data.frame(score = -obs$score, bad = obs$bad)
  result         <- calibrate(flipped, increasing = TRUE)
  rates          <- result$n_bads / result$n_obs

  expect_equal(rates, sort(rates))
  expect_equal(sum(result$n_obs), nrow(obs))
})
