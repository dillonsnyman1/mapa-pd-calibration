# Monotone Adjacent Pooling Algorithm (MAPA)
#
# MAPA turns raw (score, bad) observations into a score-to-PD calibration
# curve that is guaranteed to be monotone: as the score improves, the
# calibrated PD never gets worse.
#
# Bins are represented as data.frames with columns:
#   score_min, score_max, n_obs, n_bads, count, count_bads
# For unweighted data, count == n_obs and count_bads == n_bads.
# For value-weighted data, n_obs and n_bads are weighted sums while
# count and count_bads are raw observation/default counts.
# Calibrated bins also have a `pd` column.
#
# Run the full pipeline with run_pipeline(), or call the steps individually.
# See the README for usage examples.

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

.bad_rate <- function(bin) {
  bin$n_bads / bin$n_obs
}

# Whether `upper` (the higher-scoring bin) violates monotonicity relative
# to `lower`.
.violates <- function(lower, upper, increasing) {
  if (increasing) {
    .bad_rate(upper) < .bad_rate(lower)
  } else {
    .bad_rate(upper) > .bad_rate(lower)
  }
}

# Two-proportion z-test: TRUE means the bad rates of a and b are NOT
# significantly different at the given confidence level (so merge them).
# Uses raw counts (count, count_bads) for sample sizes, matching the
# Python reference. bad_rate (n_bads/n_obs) is still used for the rate
# difference in the numerator.
.not_significant <- function(a, b, confidence) {
  pooled_rate <- (a$count_bads + b$count_bads) / (a$count + b$count)
  if (pooled_rate <= 0 || pooled_rate >= 1) {
    return(TRUE)
  }
  se <- sqrt(pooled_rate * (1 - pooled_rate) * (1 / a$count + 1 / b$count))
  z <- abs(.bad_rate(a) - .bad_rate(b)) / se
  z_critical <- qnorm((1 + confidence) / 2)
  z < z_critical
}

# Merge two adjacent bins into one spanning their combined score range.
.merge_bins <- function(a, b) {
  data.frame(
    score_min  = a$score_min,
    score_max  = b$score_max,
    n_obs      = a$n_obs + b$n_obs,
    n_bads     = a$n_bads + b$n_bads,
    count      = a$count + b$count,
    count_bads = a$count_bads + b$count_bads,
    stringsAsFactors = FALSE
  )
}

# Whether `upper` violates monotonicity of pd relative to `lower`.
.violates_pd <- function(lower, upper, increasing) {
  if (increasing) {
    upper$pd < lower$pd
  } else {
    upper$pd > lower$pd
  }
}

# Merge two adjacent calibrated bins (n_obs-weighted pd).
.merge_calibrated <- function(a, b) {
  n_obs  <- a$n_obs + b$n_obs
  n_bads <- a$n_bads + b$n_bads
  data.frame(
    score_min  = a$score_min,
    score_max  = b$score_max,
    n_obs      = n_obs,
    n_bads     = n_bads,
    pd         = (a$pd * a$n_obs + b$pd * b$n_obs) / n_obs,
    count      = a$count + b$count,
    count_bads = a$count_bads + b$count_bads,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Group raw (score, bad) observations into one bin per unique score,
#' ordered by score ascending.
#'
#' @param observations A data.frame or matrix with columns `score` and `bad`,
#'   and an optional third column `weight`. Or a list of two-/three-element
#'   vectors c(score, bad [, weight]). The `bad` column should be 1 for a
#'   default and 0 otherwise. When `weight` is provided, `n_obs` and `n_bads`
#'   become weighted sums; `count` and `count_bads` track the raw (unweighted)
#'   observation and default counts. When no weight column is present,
#'   weight = 1 for every observation (so count == n_obs, count_bads == n_bads).
#' @return A data.frame with columns score_min, score_max, n_obs, n_bads,
#'   count, count_bads — one row per unique score value.
bins_from_observations <- function(observations) {
  if (is.data.frame(observations)) {
    scores  <- observations$score
    bads    <- observations$bad
    weights <- if (!is.null(observations$weight)) observations$weight else NULL
  } else {
    # Assume matrix or list of pairs/triples
    observations <- as.data.frame(observations)
    scores  <- observations[[1]]
    bads    <- observations[[2]]
    weights <- if (ncol(observations) >= 3) observations[[3]] else NULL
  }

  if (is.null(weights)) {
    weights <- rep(1, length(scores))
  }

  unique_scores <- sort(unique(scores))
  rows <- lapply(unique_scores, function(s) {
    idx <- scores == s
    data.frame(
      score_min  = s,
      score_max  = s,
      n_obs      = sum(weights[idx]),
      n_bads     = sum(as.integer(bads[idx]) * weights[idx]),
      count      = sum(idx),
      count_bads = sum(as.integer(bads[idx])),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Run the Monotone Adjacent Pooling Algorithm (PAVA-style).
#'
#' @param bins A data.frame of bins (score_min, score_max, n_obs, n_bads),
#'   ordered by score ascending (e.g. from bins_from_observations()).
#' @param increasing If FALSE (default), bad rate must be non-increasing.
#'   If TRUE, bad rate must be non-decreasing.
#' @param min_confidence Optional numeric in (0,1). Adjacent bins whose bad
#'   rates are not distinguishable at this confidence level (two-proportion
#'   z-test) are merged even if they don't violate monotonicity.
#' @return A data.frame of pooled bins satisfying monotonicity.
mapa <- function(bins, increasing = FALSE, min_confidence = NULL) {
  # Stack implemented as a list of single-row data.frames
  stack <- list()

  for (i in seq_len(nrow(bins))) {
    b <- bins[i, , drop = FALSE]
    stack <- c(stack, list(b))

    repeat {
      n <- length(stack)
      if (n < 2) break

      lower <- stack[[n - 1]]
      upper <- stack[[n]]

      should_merge <- .violates(lower, upper, increasing) ||
        (!is.null(min_confidence) && .not_significant(lower, upper, min_confidence))

      if (!should_merge) break

      merged <- .merge_bins(lower, upper)
      stack[[n - 1]] <- merged
      stack[[n]] <- NULL
    }
  }

  do.call(rbind, stack)
}

#' Convenience wrapper: bins_from_observations + mapa.
#'
#' @param observations Raw observations; see bins_from_observations().
#' @param increasing See mapa().
#' @param min_confidence See mapa().
#' @return A data.frame of pooled bins.
calibrate <- function(observations, increasing = FALSE, min_confidence = NULL) {
  mapa(bins_from_observations(observations), increasing, min_confidence)
}

#' Pool bins below minimum size thresholds, then re-run mapa.
#'
#' A bin "violates" if its observation or bad measure falls below the
#' threshold. Each violating bin is merged into whichever adjacent bin
#' has the closer bad rate.
#'
#' @param bins A data.frame of bins, typically from mapa().
#' @param min_obs Minimum observations required per bin (default 0).
#' @param min_bads Minimum bads required per bin (default 0).
#' @param increasing See mapa().
#' @param min_confidence See mapa().
#' @param use_counts If TRUE (default), check thresholds against raw
#'   observation/bad counts (count, count_bads). If FALSE, check against
#'   weighted sums (n_obs, n_bads). For number-weighted data both are
#'   identical.
#' @return A data.frame of pooled bins, monotone.
enforce_minimum_size <- function(bins, min_obs = 0, min_bads = 0,
                                  increasing = FALSE, min_confidence = NULL,
                                  use_counts = TRUE) {
  # Work on a list of single-row data.frames for easy splicing
  bin_list <- lapply(seq_len(nrow(bins)), function(i) bins[i, , drop = FALSE])

  repeat {
    n <- length(bin_list)
    if (n <= 1) break

    # Find first violating bin
    violator <- NULL
    for (i in seq_len(n)) {
      b <- bin_list[[i]]
      obs_val  <- if (use_counts) b$count else b$n_obs
      bads_val <- if (use_counts) b$count_bads else b$n_bads
      if (obs_val < min_obs || bads_val < min_bads) {
        violator <- i
        break
      }
    }

    if (is.null(violator)) break

    if (violator == 1) {
      neighbour <- 2
    } else if (violator == n) {
      neighbour <- n - 1
    } else {
      rate      <- .bad_rate(bin_list[[violator]])
      left_diff  <- abs(rate - .bad_rate(bin_list[[violator - 1]]))
      right_diff <- abs(rate - .bad_rate(bin_list[[violator + 1]]))
      neighbour <- if (left_diff <= right_diff) violator - 1 else violator + 1
    }

    lo <- min(violator, neighbour)
    hi <- max(violator, neighbour)
    merged <- .merge_bins(bin_list[[lo]], bin_list[[hi]])

    # Rebuild list: everything before lo, merged bin, everything after hi
    before <- if (lo > 1)  bin_list[seq_len(lo - 1)]   else list()
    after  <- if (hi < n)  bin_list[seq(hi + 1, n)]    else list()
    bin_list <- c(before, list(merged), after)
  }

  result_bins <- do.call(rbind, bin_list)
  mapa(result_bins, increasing, min_confidence)
}

#' Shrink each bin's empirical bad rate toward a prior (Bayesian credibility).
#'
#' Each bin's adjusted PD:
#'   pd = (n_bads + k * prior) / (n_obs + k)
#'
#' @param bins A data.frame of bins, typically from mapa().
#' @param k Credibility weight (equivalent observations of the prior).
#' @param prior PD to shrink toward. Defaults to overall bad rate.
#' @return A data.frame with an additional `pd` column.
apply_bayesian_adjustment <- function(bins, k, prior = NULL) {
  if (is.null(prior)) {
    prior <- sum(bins$n_bads) / sum(bins$n_obs)
  }
  pd <- (bins$n_bads + k * prior) / (bins$n_obs + k)
  cbind(bins, pd = pd)
}

#' Re-pool calibrated bins to restore monotonicity of pd.
#'
#' Runs the same adjacent-pooling algorithm as mapa(), but on `pd` instead
#' of bad_rate, merging by n_obs-weighted average pd.
#'
#' @param bins A data.frame with pd column, typically from apply_bayesian_adjustment().
#' @param increasing See mapa().
#' @return A data.frame of pooled calibrated bins with monotone pd.
repool_calibrated_bins <- function(bins, increasing = FALSE) {
  stack <- list()

  for (i in seq_len(nrow(bins))) {
    b <- bins[i, , drop = FALSE]
    stack <- c(stack, list(b))

    repeat {
      n <- length(stack)
      if (n < 2) break
      if (!.violates_pd(stack[[n - 1]], stack[[n]], increasing)) break

      merged <- .merge_calibrated(stack[[n - 1]], stack[[n]])
      stack[[n - 1]] <- merged
      stack[[n]] <- NULL
    }
  }

  do.call(rbind, stack)
}

#' Interpolate a smoothed PD for an individual score via log-odds interpolation.
#'
#' Each pool is anchored at its midpoint score. log-odds is linearly
#' interpolated between the two bracketing midpoints, then converted back to
#' a probability. Flat extrapolation beyond the first/last midpoint.
#'
#' @param bins A data.frame with pd column (from repool_calibrated_bins()).
#' @param score The score to compute a PD for.
#' @return A single numeric PD value.
interpolate_pd <- function(bins, score) {
  mids     <- (bins$score_min + bins$score_max) / 2
  log_odds <- log((1 - bins$pd) / bins$pd)
  n        <- nrow(bins)

  if (score <= mids[1])  return(bins$pd[1])
  if (score >= mids[n])  return(bins$pd[n])

  for (i in seq_len(n - 1)) {
    if (mids[i] <= score && score <= mids[i + 1]) {
      t            <- (score - mids[i]) / (mids[i + 1] - mids[i])
      interpolated <- log_odds[i] + t * (log_odds[i + 1] - log_odds[i])
      return(1 / (1 + exp(interpolated)))
    }
  }

  stop("unreachable: score must lie between mids[1] and mids[n] here")
}

#' Run the full MAPA pipeline.
#'
#' Chains: bins_from_observations -> mapa -> enforce_minimum_size ->
#'         apply_bayesian_adjustment -> repool_calibrated_bins
#'
#' @param observations Raw observations; see bins_from_observations().
#' @param k Bayesian credibility weight; see apply_bayesian_adjustment().
#' @param min_obs Minimum observations per bin; see enforce_minimum_size().
#' @param min_bads Minimum bads per bin; see enforce_minimum_size().
#' @param prior PD to shrink toward; see apply_bayesian_adjustment().
#' @param increasing Direction of monotonicity; see mapa().
#' @param min_confidence Confidence-based pooling threshold; see mapa().
#' @param use_counts See enforce_minimum_size().
#' @return A list with:
#'   $bands       — data.frame of final calibrated bins
#'   $pd_for_score — function(score) returning a smoothed PD
run_pipeline <- function(observations, k, min_obs = 0, min_bads = 0,
                          prior = NULL, increasing = FALSE,
                          min_confidence = NULL, use_counts = TRUE) {
  pooled     <- calibrate(observations, increasing, min_confidence)
  sized      <- enforce_minimum_size(pooled, min_obs, min_bads, increasing,
                                     min_confidence, use_counts)
  calibrated <- apply_bayesian_adjustment(sized, k, prior)
  bands      <- repool_calibrated_bins(calibrated, increasing)

  list(
    bands         = bands,
    pd_for_score  = function(score) interpolate_pd(bands, score)
  )
}
