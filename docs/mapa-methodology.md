# MAPA: Monotone Adjacent Pooling Algorithm

## Background

When a credit score is calibrated to a probability of default (PD), the
resulting score-to-PD mapping is expected to be monotone: a higher
score should never imply a higher PD than a lower score. This is both a
business expectation (the score is supposed to rank-order risk) and, in
most jurisdictions, a regulatory one.

In practice, empirical bad rates observed within score bands are noisy,
especially in the tails where observation counts are small. A direct
mapping from "observed bad rate per score band" to PD can therefore be
non-monotone even when the underlying score itself rank-orders risk
well overall.

The Monotone Adjacent Pooling Algorithm (MAPA) resolves this by merging
("pooling") adjacent score bands wherever the observed bad rates violate
monotonicity, until the resulting step function is monotone. Each pooled
band's bad rate then becomes the calibrated PD for every score within it.

## Algorithm

**Input**: raw observations, each a `(score, bad)` pair (or
`(score, bad, weight)` - see "Number vs. value weighting" below), where
`bad` is 1 for a default and 0 otherwise, and `weight` (default 1) is an
optional observation weight (e.g. exposure at default for value-weighted
calibration).

**Output**: a set of pooled score bands, in score order, partitioning the
full observed score range, whose bad rates are monotone (by default,
non-increasing as score increases).

**Procedure**:

1. **Bin**: group the raw observations into one band per unique score,
   ordered by score ascending. This is the finest possible starting point
   - every distinct score value is initially its own band. Each bin tracks
   both weighted sums (`n_obs` = sum of weights, `n_bads` = sum of
   weights for defaulted observations) and raw counts (`count` = number
   of observations, `count_bads` = number of defaulted observations).
   For number-weighted data (all weights = 1), these are identical.
2. **Pool**: process bands left to right (lowest score first), maintaining
   a stack of pooled bands built so far.
   - Push the next band onto the stack.
   - While the top two bands on the stack violate the required
     monotonicity (i.e. the higher-scoring band has a *worse* bad rate than
     the lower-scoring band it's compared against), merge them into a
     single band by summing their weighted sums and raw counts, and push
     the merged band back onto the stack.
   - Repeat until all bands have been processed.

Because merging recomputes the bad rate (`n_bads / n_obs`) of the combined
band, a single merge can resolve (or trigger) multiple violations, so the
pooling step repeats until the top of the stack is consistent with the band
below it. When value-weighted, this bad rate is an exposure-weighted default
rate.

If you already have pre-aggregated score bands (e.g. from an existing
binning scheme) rather than raw observations, you can skip the binning step
and run the pooling step directly on those bands.

## Relationship to PAVA / isotonic regression

MAPA is an application of the Pool Adjacent Violators Algorithm (PAVA)
(the classical algorithm for isotonic regression) to grouped, binary
outcome data ordered by score. PAVA finds the monotone step function that
best fits the data in a least-squares sense; pooling adjacent violators
until none remain is exactly how that optimum is reached. MAPA applies the
same mechanics with credit-scoring terminology (score bands, bad rates,
PD) and framing (producing a usable score-to-PD lookup table).

## Using the output

The pooled bands form a step function: for any score, find the pooled
band whose original score range contains it, and use that band's bad rate
as the calibrated PD. Because pooling only merges *adjacent* bands and
preserves total observation and bad counts, the pooled bands still
partition the full score range and the full population exactly.

## Minimum bin size

Pooling stops as soon as monotonicity is satisfied, it does not guarantee
that every resulting band has *enough* observations to support a stable
rate estimate. A band can end up monotone but tiny (e.g. a handful of
observations at the extreme low- or high-risk end of the score range).

`enforce_minimum_size` addresses this with a further pooling pass applied
after `mapa`: any band with fewer than `min_obs` observations or fewer than
`min_bads` bads is merged into whichever adjacent band has the *closer bad
rate* - minimizing the distortion this introduces - and the process
repeats until every band meets both thresholds (or only one band remains).

The `use_counts` parameter (default `true`) controls what the thresholds
are checked against. When `true`, thresholds compare against raw
observation counts (`count`, `count_bads`) - this is the natural choice
for statistical reliability, since each raw observation is an independent
trial. When `false`, thresholds compare against weighted sums (`n_obs`,
`n_bads`) - useful when thresholds represent economic significance (e.g.
minimum exposure volume). For number-weighted data the two are identical.

Because merging toward the closer-rate neighbour can introduce a new
monotonicity violation (the merged band's rate moves toward, and can cross,
its other neighbour), the result is passed back through `mapa` to restore
monotonicity before being returned.

With the default thresholds (`min_obs = 0, min_bads = 0`), this step is a
no-op - every band trivially satisfies the thresholds regardless of
whether counts or weighted sums are used.

## Confidence-based pooling

Monotonicity alone does not guarantee that two adjacent bands' bad rates
are *meaningfully different* as a band can be merely-monotone-by-luck rather
than genuinely distinguishable from its neighbour, especially when either
band is small.

`mapa` (and therefore `calibrate`, `enforce_minimum_size`, and
`run_pipeline`) accepts an optional `min_confidence` parameter (e.g. `0.95`
for 95%). When given, pooling also merges adjacent bands whose bad rates are
not statistically significantly different at that confidence level, even if
they don't violate monotonicity.

The test used is the standard two-proportion z-test. The bad rates being
compared are `p_a = n_bads_a / n_obs_a` and `p_b = n_bads_b / n_obs_b`
(which may be value-weighted), but the sample sizes in the test use raw
observation counts (`count_a`, `count_b`), not weighted sums. This is
because the z-test assumes independent Bernoulli trials - weighted sums
don't represent independent trials and would make the test statistically
meaningless. The pooled rate for the test is
`p = (n_bads_a + n_bads_b) / (n_obs_a + n_obs_b)`, giving a standard
error

```
se = sqrt(p * (1 - p) * (1 / count_a + 1 / count_b))
```

and a z-statistic `z = |p_a - p_b| / se`. If `z` is below the critical value
for the requested confidence level (e.g. `1.96` for 95%), the difference is
not significant and the bands are merged - using the same merge rule as
ordinary pooling (summing weighted sums and raw counts).

This is applied during the same left-to-right pooling pass as the
monotonicity check, so the result remains monotone: merging two adjacent
bands always produces a rate between the two, which cannot introduce a new
violation against either band's other neighbour.

With the default (`min_confidence` not given), only monotonicity violations
are merged, exactly as described above. See
`reference/fixtures/expected_pooled_bins_confidence.csv` for the result of
applying `min_confidence=0.95` to the bundled example - the 41 per-score
bands collapse to 6 bands whose bad rates are pairwise distinguishable at
the 95% confidence level.

## Bayesian adjustment

Pooling fixes *monotonicity*, but not *noise*: a pooled band can still
have very few observations, especially near the tails of the score
distribution, making its empirical bad rate an unreliable PD estimate on
its own.

`apply_bayesian_adjustment` addresses this with a standard credibility
(Bayesian shrinkage) estimate: each band's empirical bad rate is pulled
toward a prior - typically the overall population bad rate - by an amount
that depends on the band's size:

```
pd = (n_bads + k * prior) / (n_obs + k)
```

`k` is the credibility weight, expressed as a number of "equivalent
observations" of the prior. A band with `n_obs == k` is shrunk halfway
between its own bad rate and the prior; bands much larger than `k` are
barely adjusted, and bands much smaller than `k` end up close to the
prior. When value-weighted, `n_obs` and `n_bads` in the formula are
weighted sums, so the shrinkage strength is determined by total weighted
volume rather than raw count. The `k` parameter is still denominated in
the same units as `n_obs` (i.e. total weight, not number of observations).

### Note on monotonicity after shrinkage

Shrinking each band independently toward a single global prior is not
guaranteed to preserve the monotonicity established by the pooling step. A
small band whose empirical rate happens to be close to a much larger
neighbour's can be pulled past it toward the prior, re-introducing a small
local violation.

The bundled reference fixtures intentionally include such a case (two
adjacent bands near the low-risk end of the score range, with very
similar bad rates but very different sizes) to make this concrete -
running `apply_bayesian_adjustment` on the size-enforced bands produces a
pair of adjacent PDs that are no longer non-increasing.

In practice this is usually addressed by one or more of:

- choosing `k` conservatively relative to typical band sizes,
- re-running the pooling step on the *adjusted* rates (treating each band's
  shrunk PD and observation count as a new input bin) to restore strict
  monotonicity, or
- accepting small, economically immaterial violations in low-risk bands
  where the practical impact on capital/provisions is negligible.

This reference implementation takes the second approach via
`repool_calibrated_bins`, described next.

## Re-pooling after shrinkage

`repool_calibrated_bins` re-applies the pooling step to the
Bayesian-adjusted bands, this time using each band's shrunk `pd` (rather
than its raw bad rate) as the quantity that must be monotone. It is
otherwise the same algorithm as `mapa`: adjacent bands whose `pd` violates
the required monotonicity are merged, repeatedly, until the sequence of
`pd` values is monotone.

The only difference is how merged bands are combined. `mapa` sums
observation and bad counts and lets the bad rate fall out of those sums.
Here there is no single underlying count to sum - `pd` is already a shrunk
estimate - so merging instead takes the `n_obs`-weighted average of the two
bands' `pd` values:

```
pd = (pd_a * n_obs_a + pd_b * n_obs_b) / (n_obs_a + n_obs_b)
```

Weighted sums and raw counts are still summed as before, so totals
continue to be preserved.

Running `repool_calibrated_bins` on the bundled fixture's calibrated bins
merges the two crossing bands from the example above into one, restoring a
fully non-increasing `pd` sequence (see
`reference/fixtures/expected_repooled_calibrated_bins.csv`).

## Smoothing: log-odds interpolation

Even after re-pooling, the calibrated PD curve is a step function:
every score within a pool gets that pool's single PD, with a discontinuous
jump at each pool boundary. Two scores one point apart, either side of a
pool boundary, can get noticeably different PDs, while two scores far apart
within the same wide pool get identical PDs.

`interpolate_pd` smooths this into a continuous curve using log-odds
interpolation, a standard scorecard scaling technique. Each pool is
reduced to a single anchor point:

- its midpoint score, `(score_min + score_max) / 2`, and
- the log-odds of its `pd`: `log_odds = ln((1 - pd) / pd)`.

For an individual score, log-odds is linearly interpolated between the
anchor points of the two pools whose midpoints bracket it, then converted
back to a PD via `pd = 1 / (1 + exp(log_odds))`. Scores at or beyond the
first or last pool's midpoint get that pool's PD unchanged (flat
extrapolation).

Because log-odds is a monotonic transform of `pd`, interpolating in
log-odds space preserves the monotonicity of the underlying pool PDs - the
smoothed curve is monotone wherever the pool PDs are.

This step is purely cosmetic: it doesn't change which pool a score
"belongs to" for reporting purposes, only the PD assigned to individual
scores within a pool. It is typically applied last, after
`repool_calibrated_bins`.

## Number vs. value weighting

By default, every observation carries equal weight (`weight = 1`), and the
pipeline produces a standard number-weighted PD calibration: bad rates are
simple proportions, and all size thresholds refer to observation counts.
This is the traditional approach for scorecard calibration.

For some applications - notably IFRS 9 expected credit loss (ECL)
calculations - PD must be calibrated on an exposure-weighted basis, where
each observation's contribution is proportional to its exposure at default
(EAD). Setting `weight` to the observation's EAD produces a
**value-weighted** calibration: bad rates become exposure-weighted default
rates, and weighted sums (`n_obs`, `n_bads`) reflect total exposure volume
rather than counts.

How weighting flows through the pipeline:

- **Binning and pooling** use `n_bads / n_obs` as the bad rate, which is
  exposure-weighted when weights are used. Monotonicity is enforced on
  this weighted bad rate.
- **Confidence-based pooling** uses raw counts (`count`, `count_bads`) for
  sample sizes in the z-test, since the statistical test requires
  independent trials. The bad rates being compared are still
  `n_bads / n_obs`.
- **Minimum bin size** (`enforce_minimum_size`) defaults to checking raw
  counts (`use_counts = true`), but accepts `use_counts = false` to check
  weighted sums instead. `run_pipeline` passes `use_counts` through.
- **Bayesian adjustment** uses `n_obs` and `n_bads` (weighted sums) in the
  credibility formula. The `k` parameter should be set relative to typical
  `n_obs` values (total weight per band, not observation count).
- **Re-pooling** and **smoothing** operate on the adjusted `pd` values and
  are unaffected by the weighting mode.

## Attribution

MAPA, as implemented here, is based on the method described by
**Raymond Anderson** in *The Credit Scoring Toolkit* (Oxford University
Press, 2007) and related published work on score calibration. Anderson
developed and applied this approach during his work in retail credit risk,
including at institutions where the author of this repository has worked.
All credit for the underlying methodology - pooling, minimum bin size,
Bayesian/credibility adjustment, and log-odds smoothing - belongs to
Anderson's published work.

The confidence-based pooling option (see above) is the one piece of this
implementation the author isn't aware of from Anderson's published work
specifically - it's a standard statistical technique (a two-proportion
z-test) applied to the same pooling step, included here as an optional
extension rather than as an attributed part of the original method.

This repository is an independent, from-scratch implementation using
synthetic data, intended as an educational reference. It is not affiliated
with, endorsed by, or derived from any proprietary source code.
