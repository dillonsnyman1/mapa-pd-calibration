/******************************************************************************
 Monotone Adjacent Pooling Algorithm (MAPA) - SAS macro reference
 implementation.

 Mirrors the Python (reference/python/mapa.py) and C++
 (reference/cpp/mapa.hpp / mapa.cpp) reference implementations:

   %mapa_bins_from_observations  - group raw (score, bad) observations into
                                    one bin per unique score
   %mapa_pool                    - pool adjacent bins until the bad rate is
                                    monotone (the core PAVA-style algorithm),
                                    optionally also pooling bins whose bad
                                    rates aren't statistically distinguishable
   %mapa_calibrate                - bins_from_observations + pool
   %mapa_enforce_minimum_size     - further pool bins below given size
                                    thresholds, then restore monotonicity
   %mapa_bayesian_adjustment       - shrink each bin's bad rate toward a
                                    prior using credibility weighting
   %mapa_repool_calibrated         - re-pool Bayesian-adjusted bins to
                                    restore monotonicity of pd
   %mapa_interpolate_pd             - smooth the pooled PD step function via
                                    log-odds interpolation between pools
   %mapa_run_pipeline               - run the full pipeline above in one
                                    call, producing both a band table and
                                    (optionally) smoothed per-score PDs

 All datasets use the columns: score_min, score_max, n_obs, n_bads, count,
 count_bads (and, for %mapa_bayesian_adjustment's and
 %mapa_repool_calibrated's output, pd). For unweighted data n_obs = count
 and n_bads = count_bads; for value-weighted data n_obs = sum(weight) and
 n_bads = sum(bad * weight) while count and count_bads track raw
 observation numbers.

 Implementation style: rather than simulating a stack with _temporary_
 arrays, %mapa_pool and %mapa_enforce_minimum_size are implemented as
 macro-driven "iterate until converged" loops - each pass is a plain DATA
 step + PROC SQL, and the loop stops once a pass makes no further changes.
 This keeps every individual step close to ordinary SAS data manipulation,
 at the cost of needing several passes over (typically very small) bin
 datasets.

 Note: these macros create working datasets in WORK prefixed with
 _mapa_pool_ and _mapa_sized_, which are deleted again before each macro
 returns. Avoid using those prefixes for your own dataset names.

 See ../../docs/mapa-methodology.md for background and attribution.
******************************************************************************/

%macro mapa_bins_from_observations(in=, out=);
    /* Group raw (score, bad) observations into one bin per unique score,
       ordered by score ascending. `in` must have columns: score, bad
       (bad = 1 for a default, 0 otherwise). If `in` also has a column
       named weight, n_obs = sum(weight) and n_bads = sum(bad * weight);
       otherwise weight defaults to 1 (so n_obs = count(*) and
       n_bads = sum(bad)). count and count_bads always track the raw
       number of observations regardless of weighting. */
    %local has_weight;
    proc sql noprint;
        select count(*) into :has_weight
        from dictionary.columns
        where libname='WORK' and upcase(memname)=upcase("&in") and upcase(name)='WEIGHT';
    quit;

    proc sql;
        create table &out as
        select score as score_min,
               score as score_max,
               %if &has_weight > 0 %then %do;
                   sum(weight) as n_obs,
                   sum(bad * weight) as n_bads,
               %end;
               %else %do;
                   count(*) as n_obs,
                   sum(bad) as n_bads,
               %end;
               count(*) as count,
               sum(bad) as count_bads
        from &in
        group by score_min
        order by score_min;
    quit;
%mend mapa_bins_from_observations;


%macro mapa_pool(in=, out=, increasing=0, min_confidence=0);
    /* Pool adjacent bins until the bad rate is monotone.

       `in` must be sorted by score_min ascending, with columns:
       score_min, score_max, n_obs, n_bads.

       If increasing=0 (default), the bad rate must be non-increasing as
       score increases (higher score = lower risk). If increasing=1, the
       bad rate must be non-decreasing instead.

       If min_confidence > 0 (e.g. 0.95 for 95%), adjacent bins whose bad
       rates do not differ at this confidence level (two-proportion z-test)
       are merged as well, even if they don't violate monotonicity. This
       produces fewer, larger bins whose bad rates are more reliably
       distinguishable from their neighbours. With the default (0), only
       monotonicity violations are merged.

       Each pass assigns every row to a group: a row starts a new group
       unless its bad rate violates monotonicity relative to the running
       (possibly already-merged) rate of the current group, or (if
       min_confidence > 0) is statistically indistinguishable from it, in
       which case it is merged into that group. Groups are then
       re-aggregated. This repeats until a pass produces no merges (row
       count unchanged), which happens exactly when no violations or
       insignificant differences remain. */
    %local iter converged n_before n_after src z_crit;
    %let iter = 0;
    %let converged = 0;
    %let src = &in;

    %if &min_confidence > 0 %then
        %let z_crit = %sysfunc(probit(%sysevalf((1 + &min_confidence) / 2)));
    %else
        %let z_crit = 0;

    %do %until (&converged);
        %let iter = %eval(&iter + 1);

        data _mapa_pool_grp;
            set &src;
            retain _grp 1 _run_nobs _run_nbads _run_count _run_count_bads;
            if _n_ = 1 then do;
                _run_nobs  = n_obs;
                _run_nbads = n_bads;
                _run_count = count;
                _run_count_bads = count_bads;
            end;
            else do;
                _violates = (&increasing = 0 and (n_bads / n_obs) > (_run_nbads / _run_nobs))
                            or (&increasing = 1 and (n_bads / n_obs) < (_run_nbads / _run_nobs));

                if &min_confidence > 0 then do;
                    _p_pool = (_run_count_bads + count_bads) / (_run_count + count);
                    if _p_pool <= 0 or _p_pool >= 1 then _not_sig = 1;
                    else do;
                        _se = sqrt(_p_pool * (1 - _p_pool) * (1 / _run_count + 1 / count));
                        _z = abs((n_bads / n_obs) - (_run_nbads / _run_nobs));
                        _not_sig = (_z / _se < &z_crit);
                    end;
                end;
                else _not_sig = 0;

                if _violates or _not_sig then do;
                    _run_nobs  + n_obs;
                    _run_nbads + n_bads;
                    _run_count + count;
                    _run_count_bads + count_bads;
                end;
                else do;
                    _grp + 1;
                    _run_nobs  = n_obs;
                    _run_nbads = n_bads;
                    _run_count = count;
                    _run_count_bads = count_bads;
                end;
            end;
        run;

        proc sql noprint;
            select count(*) into :n_before from &src;
        quit;

        proc sql;
            create table _mapa_pool_&iter as
            select min(score_min) as score_min,
                   max(score_max) as score_max,
                   sum(n_obs) as n_obs,
                   sum(n_bads) as n_bads,
                   sum(count) as count,
                   sum(count_bads) as count_bads
            from _mapa_pool_grp
            group by _grp
            order by score_min;
        quit;

        proc sql noprint;
            select count(*) into :n_after from _mapa_pool_&iter;
        quit;

        %if &n_after = &n_before %then %let converged = 1;
        %let src = _mapa_pool_&iter;
    %end;

    data &out;
        set &src;
        keep score_min score_max n_obs n_bads count count_bads;
    run;

    proc datasets lib=work nolist;
        delete _mapa_pool_: ;
    quit;
%mend mapa_pool;


%macro mapa_calibrate(in=, out=, increasing=0, min_confidence=0);
    /* Convenience: group raw (score, bad) observations into per-score
       bins, then pool them. */
    %mapa_bins_from_observations(in=&in, out=_mapa_initial_bins)
    %mapa_pool(in=_mapa_initial_bins, out=&out, increasing=&increasing, min_confidence=&min_confidence)

    proc datasets lib=work nolist;
        delete _mapa_initial_bins;
    quit;
%mend mapa_calibrate;


%macro mapa_enforce_minimum_size(in=, out=, min_obs=0, min_bads=0, increasing=0, min_confidence=0, use_counts=1);
    /* Further pool bins that don't meet minimum size thresholds, even if
       they don't violate monotonicity.

       `in` must be sorted by score_min ascending and already pooled (e.g.
       the output of %mapa_pool), with columns: score_min, score_max,
       n_obs, n_bads, count, count_bads.

       When use_counts=1 (the default), a bin "violates" if count < min_obs
       or count_bads < min_bads (raw observation counts). When use_counts=0,
       n_obs and n_bads (weighted sums) are checked instead. Each pass finds
       the first violating bin, decides which adjacent bin has the closer bad
       rate, and merges the two. This repeats until no bin violates the
       thresholds or only one bin remains. Because merging toward the
       closer-rate neighbour can re-introduce a monotonicity violation, the
       result is passed back through %mapa_pool. min_confidence is forwarded
       to that final pass; see %mapa_pool. */
    %local iter src n violator neighbour lo hi done;
    %let iter = 0;
    %let src = &in;
    %let done = 0;

    %do %until (&done);
        data _mapa_sized_seq;
            set &src;
            _seq = _n_;
        run;

        proc sql noprint;
            select count(*) into :n from _mapa_sized_seq;
        quit;

        %if &n <= 1 %then %do;
            %let src = _mapa_sized_seq;
            %let done = 1;
        %end;
        %else %do;
            /* Annotate each bin with its immediate neighbours' counts. */
            proc sql;
                create table _mapa_sized_ann as
                select a.*,
                       p.n_obs as prev_nobs, p.n_bads as prev_nbads,
                       x.n_obs as next_nobs, x.n_bads as next_nbads
                from _mapa_sized_seq a
                left join _mapa_sized_seq p on p._seq = a._seq - 1
                left join _mapa_sized_seq x on x._seq = a._seq + 1;
            quit;

            %let violator = ;
            proc sql noprint;
                select min(_seq) into :violator
                from _mapa_sized_ann
                %if &use_counts = 1 %then %do;
                    where count < &min_obs or count_bads < &min_bads;
                %end;
                %else %do;
                    where n_obs < &min_obs or n_bads < &min_bads;
                %end;
            quit;

            %if %length(&violator) = 0 %then %do;
                %let src = _mapa_sized_seq;
                %let done = 1;
            %end;
            %else %do;
                /* Pick whichever neighbour has the closer bad rate. */
                data _null_;
                    set _mapa_sized_ann;
                    where _seq = &violator;
                    rate = n_bads / n_obs;
                    if _seq = 1 then neighbour = 2;
                    else if _seq = &n then neighbour = _seq - 1;
                    else do;
                        prate = prev_nbads / prev_nobs;
                        xrate = next_nbads / next_nobs;
                        if abs(rate - prate) <= abs(rate - xrate) then neighbour = _seq - 1;
                        else neighbour = _seq + 1;
                    end;
                    call symputx('neighbour', neighbour);
                run;

                %let lo = %sysfunc(min(&violator, &neighbour));
                %let hi = %sysfunc(max(&violator, &neighbour));
                %let iter = %eval(&iter + 1);

                /* Merge bins &lo and &hi into a single bin. */
                data _mapa_sized_merged;
                    set _mapa_sized_seq(where=(_seq = &lo) keep=score_min n_obs n_bads count count_bads);
                    set _mapa_sized_seq(where=(_seq = &hi) keep=score_max n_obs n_bads count count_bads
                                         rename=(n_obs=n_obs_hi n_bads=n_bads_hi count=count_hi count_bads=count_bads_hi));
                    n_obs  = n_obs + n_obs_hi;
                    n_bads = n_bads + n_bads_hi;
                    count  = count + count_hi;
                    count_bads = count_bads + count_bads_hi;
                    drop n_obs_hi n_bads_hi count_hi count_bads_hi;
                run;

                proc sql;
                    create table _mapa_sized_&iter as
                    select score_min, score_max, n_obs, n_bads, count, count_bads
                    from _mapa_sized_seq
                    where _seq < &lo
                    union all
                    select score_min, score_max, n_obs, n_bads, count, count_bads
                    from _mapa_sized_merged
                    union all
                    select score_min, score_max, n_obs, n_bads, count, count_bads
                    from _mapa_sized_seq
                    where _seq > &hi
                    order by score_min;
                quit;

                %let src = _mapa_sized_&iter;
            %end;
        %end;
    %end;

    /* Restore monotonicity. */
    %mapa_pool(in=&src, out=&out, increasing=&increasing, min_confidence=&min_confidence)

    proc datasets lib=work nolist;
        delete _mapa_sized_: ;
    quit;
%mend mapa_enforce_minimum_size;


%macro mapa_bayesian_adjustment(in=, out=, k=, prior=);
    /* Shrink each bin's empirical bad rate toward `prior` using Bayesian
       (credibility) weighting:

           pd = (n_bads + k * prior) / (n_obs + k)

       `k` is the credibility weight, in equivalent observations of the
       prior. If `prior` is not given, it defaults to the overall bad rate
       across all bins (sum(n_bads) / sum(n_obs)).

       Note: shrinking each bin independently toward a single global prior
       is not guaranteed to preserve the monotonicity established by
       %mapa_pool - see ../../docs/mapa-methodology.md. */
    %local prior_value;

    %if %length(&prior) = 0 %then %do;
        proc sql noprint;
            select sum(n_bads) / sum(n_obs) into :prior_value
            from &in;
        quit;
    %end;
    %else %do;
        %let prior_value = &prior;
    %end;

    data &out;
        set &in;
        pd = (n_bads + &k * &prior_value) / (n_obs + &k);
    run;
%mend mapa_bayesian_adjustment;


%macro mapa_repool_calibrated(in=, out=, increasing=0);
    /* Re-apply pooling to Bayesian-adjusted bins, restoring monotonicity
       of pd.

       %mapa_bayesian_adjustment shrinks each bin's bad rate toward a
       shared prior independently, which is not guaranteed to preserve the
       monotonicity established by %mapa_pool - see
       ../../docs/mapa-methodology.md. This runs the same adjacent-pooling
       algorithm again, but on pd instead of the bad rate, merging
       violating bins by taking the n_obs-weighted average of their pd
       values.

       `in` must be sorted by score_min ascending, with columns:
       score_min, score_max, n_obs, n_bads, pd (e.g. the output of
       %mapa_bayesian_adjustment). `increasing` has the same meaning as in
       %mapa_pool, applied to pd instead of the bad rate. */
    %local iter converged n_before n_after src;
    %let iter = 0;
    %let converged = 0;
    %let src = &in;

    %do %until (&converged);
        %let iter = %eval(&iter + 1);

        data _mapa_repool_grp;
            set &src;
            retain _grp 1 _run_nobs _run_pd _run_count _run_count_bads;
            if _n_ = 1 then do;
                _run_nobs = n_obs;
                _run_pd   = pd;
                _run_count = count;
                _run_count_bads = count_bads;
            end;
            else do;
                if (&increasing = 0 and pd > _run_pd)
                   or (&increasing = 1 and pd < _run_pd) then do;
                    _run_pd = (_run_pd * _run_nobs + pd * n_obs) / (_run_nobs + n_obs);
                    _run_nobs + n_obs;
                    _run_count + count;
                    _run_count_bads + count_bads;
                end;
                else do;
                    _grp + 1;
                    _run_nobs = n_obs;
                    _run_pd   = pd;
                    _run_count = count;
                    _run_count_bads = count_bads;
                end;
            end;
        run;

        proc sql noprint;
            select count(*) into :n_before from &src;
        quit;

        proc sql;
            create table _mapa_repool_&iter as
            select min(score_min) as score_min,
                   max(score_max) as score_max,
                   sum(n_obs) as n_obs,
                   sum(n_bads) as n_bads,
                   sum(pd * n_obs) / sum(n_obs) as pd,
                   sum(count) as count,
                   sum(count_bads) as count_bads
            from _mapa_repool_grp
            group by _grp
            order by score_min;
        quit;

        proc sql noprint;
            select count(*) into :n_after from _mapa_repool_&iter;
        quit;

        %if &n_after = &n_before %then %let converged = 1;
        %let src = _mapa_repool_&iter;
    %end;

    data &out;
        set &src;
        keep score_min score_max n_obs n_bads pd count count_bads;
    run;

    proc datasets lib=work nolist;
        delete _mapa_repool_: ;
    quit;
%mend mapa_repool_calibrated;


%macro mapa_interpolate_pd(bins=, scores=, out=);
    /* Smooth the pooled PD step function via log-odds interpolation
       between pools.

       `bins` must be sorted by score_min ascending, with columns:
       score_min, score_max, n_obs, n_bads, pd (e.g. the output of
       %mapa_repool_calibrated). `scores` must have a single column,
       score, with one row per individual score to compute a smoothed PD
       for. `out` gets columns: score, pd.

       Each pool is represented by a single anchor point: its midpoint
       score, (score_min + score_max) / 2, and the log-odds of its pd,
       log_odds = log((1 - pd) / pd). For each value in `scores`:

         - if it is at or before the first pool's midpoint, or at or after
           the last pool's midpoint, the nearest pool's pd is returned
           unchanged (flat extrapolation beyond the anchor points).
         - otherwise, log_odds is linearly interpolated between the
           midpoints of the two pools bracketing the score, and converted
           back to a pd via pd = 1 / (1 + exp(log_odds)).

       Because log-odds is a monotonic transform of pd, this preserves the
       monotonicity of a monotone input sequence of pool PDs. */
    %local min_mid max_mid first_log_odds last_log_odds;

    data _mapa_interp_bins;
        set &bins;
        _mid = (score_min + score_max) / 2;
        _log_odds = log((1 - pd) / pd);
        keep _mid _log_odds;
    run;

    proc sql noprint;
        select min(_mid), max(_mid) into :min_mid, :max_mid from _mapa_interp_bins;
        select _log_odds into :first_log_odds from _mapa_interp_bins where _mid = &min_mid;
        select _log_odds into :last_log_odds from _mapa_interp_bins where _mid = &max_mid;
    quit;

    /* For each score, find the nearest pool midpoint at or below it (lo)
       and at or above it (hi). */
    proc sql;
        create table _mapa_interp_join as
        select s.score,
               lo._mid as _lo_mid, lo._log_odds as _lo_log_odds,
               hi._mid as _hi_mid, hi._log_odds as _hi_log_odds
        from &scores s
        left join _mapa_interp_bins lo
            on lo._mid = (select max(_mid) from _mapa_interp_bins where _mid <= s.score)
        left join _mapa_interp_bins hi
            on hi._mid = (select min(_mid) from _mapa_interp_bins where _mid >= s.score);
    quit;

    data &out;
        set _mapa_interp_join;
        if score <= &min_mid then _log_odds = &first_log_odds;
        else if score >= &max_mid then _log_odds = &last_log_odds;
        else if _lo_mid = _hi_mid then _log_odds = _lo_log_odds;
        else _log_odds = _lo_log_odds
                         + (score - _lo_mid) / (_hi_mid - _lo_mid) * (_hi_log_odds - _lo_log_odds);
        pd = 1 / (1 + exp(_log_odds));
        keep score pd;
    run;

    proc datasets lib=work nolist;
        delete _mapa_interp_bins _mapa_interp_join;
    quit;
%mend mapa_interpolate_pd;


%macro mapa_run_pipeline(in=, out_bands=, k=, min_obs=0, min_bads=0, prior=, increasing=0,
                          min_confidence=0, scores=, out_smoothed=, use_counts=1);
    /* Run the full MAPA pipeline in one call: bin, pool, enforce minimum
       size, apply Bayesian adjustment, and re-pool.

       This chains %mapa_calibrate, %mapa_enforce_minimum_size,
       %mapa_bayesian_adjustment and %mapa_repool_calibrated. The result is
       two independent outputs - use whichever suits the consumer:

         - `out_bands`: the band table (score_min, score_max, n_obs,
           n_bads, pd) - the typical deliverable for reporting and
           governance.
         - `out_smoothed` (optional): a smoothed, continuous PD per score,
           via %mapa_interpolate_pd. Only produced if both `scores` (a
           dataset with a single column, score) and `out_smoothed` are
           given.

       min_confidence is forwarded to %mapa_calibrate and
       %mapa_enforce_minimum_size; see %mapa_pool. */
    %mapa_calibrate(in=&in, out=_mapa_pipeline_pooled, increasing=&increasing, min_confidence=&min_confidence)
    %mapa_enforce_minimum_size(in=_mapa_pipeline_pooled, out=_mapa_pipeline_sized,
                                min_obs=&min_obs, min_bads=&min_bads, increasing=&increasing,
                                min_confidence=&min_confidence, use_counts=&use_counts)
    %mapa_bayesian_adjustment(in=_mapa_pipeline_sized, out=_mapa_pipeline_calibrated, k=&k, prior=&prior)
    %mapa_repool_calibrated(in=_mapa_pipeline_calibrated, out=&out_bands, increasing=&increasing)

    %if %length(&scores) > 0 and %length(&out_smoothed) > 0 %then %do;
        %mapa_interpolate_pd(bins=&out_bands, scores=&scores, out=&out_smoothed)
    %end;

    proc datasets lib=work nolist;
        delete _mapa_pipeline_pooled _mapa_pipeline_sized _mapa_pipeline_calibrated;
    quit;
%mend mapa_run_pipeline;
