function tests = test_timeseries_rev1
%TEST_TIMESERIES_REV1 Tests for plan step 6: time series in both models.
%   Covers the shared scan clock (scan_clock_rev1, scan_clock_state_rev1),
%   the event windows and library, the pair time series with the harm
%   metric (mag_pair_timeseries_rev1), and the wrapper's ts_method =
%   "events", parfor trials and random-walk drift. Plan section 5.4: the
%   event method reproduces fixed-grid maxima exactly above the floor and
%   reproduces the Rev 5 120-hour persistence tables.
%
%   Run: results = runtests("test_timeseries_rev1");
tests = functiontests(localfunctions);
end

% =========================================================================
function setupOnce(tc)
mode = struct("vcp", 212, "sails", 3, "avset_top_deg", 6.4);
tc.TestData.mode = mode;
o = mode;
o.seed = 3;
o.n_samples = 1e3;
o.duration_s = 600;
o.n_trials = 2;
tc.TestData.o = o;
tc.TestData.fixed = nexrad_asr11_mag_rev1(o);
end

% =========================================================================
function testEventsEqualFixedInWrapper(tc)
% Same clock, same grid: the event subset gives the same per-scan maxima.
d = tc.TestData;
o = d.o;
o.ts_method = "events";
e = nexrad_asr11_mag_rev1(o);
f = d.fixed.scan_max;
s = e.scan_max;
verifyLessThan(tc, s.share, 0.15, "event windows should keep a small share of the grid");
fl = max(s.floor21_dBi, s.floor12_dBi);
verifyLessThan(tc, fl, min(f.asr_rev.max_g21_dBi(f.asr_rev.complete)));
% Same rows; maxima equal above the floor, censored (NaN) at or below it.
verifyEqual(tc, s.asr_rev(:, ["trial", "rev", "complete"]), f.asr_rev(:, ["trial", "rev", "complete"]));
check_max(tc, s.asr_rev.max_g21_dBi, f.asr_rev.max_g21_dBi, fl, "revolutions");
check_max(tc, s.asr_rev.max_g12_dBi, f.asr_rev.max_g12_dBi, fl, "revolutions g12");
hi = f.asr_rev.max_g21_dBi > fl;
verifyEqual(tc, s.asr_rev.t_at_max_s(hi), f.asr_rev.t_at_max_s(hi));
verifyEqual(tc, s.nexrad_cut(:, ["trial", "volume", "cut", "el_deg", "complete"]), ...
    f.nexrad_cut(:, ["trial", "volume", "cut", "el_deg", "complete"]));
check_max(tc, s.nexrad_cut.max_g21_dBi, f.nexrad_cut.max_g21_dBi, fl, "cuts");
verifyEqual(tc, s.persistence, f.persistence);
verifyEqual(tc, s.persistence_volume, f.persistence_volume);
verifyTrue(tc, all(isnan(e.ts_pct.gain_dBi)));
end

function testPairTimeseriesMatchesWrapper(tc)
% The two-part time series (library + gains per pair) against the
% wrapper's fixed grid, at a horizontal path and at asymmetric angles.
d = tc.TestData;
w = d.fixed;
lib = scan_event_library_rev1(w.scan_max.clock_spec, struct("seed", d.o.seed, ...
    "n_trials", d.o.n_trials, "duration_s", d.o.duration_s));
R = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, 0, 0);
check_against(tc, R, w.scan_max, "0/0");
o = d.o;
o.el_nexrad_deg = 1.7;
o.el_atc_deg = 2.0;
w2 = nexrad_asr11_mag_rev1(o);
R2 = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, 1.7, 2.0);
check_against(tc, R2, w2.scan_max, "1.7/2.0");
end

function testReproducesRev5Persistence(tc)
% Plan 5.4 test 2: the Rev 5 120-hour run (worst mode, seed 20261105,
% 120 trials of 1 h; report_figures_rev5) from the event library.
d = tc.TestData;
lib = scan_event_library_rev1(d.fixed.scan_max.clock_spec, struct("seed", 20261105, ...
    "n_trials", 120, "duration_s", 3600));
R = mag_pair_timeseries_rev1(lib, d.fixed.sys_nexrad.tx, d.fixed.sys_asr11.tx, 0, 0);
here = fileparts(which("test_timeseries_rev1"));
P = readtable(fullfile(here, "docs", "figures", "persistence_120h_rev5.csv"));
verifyEqual(tc, R.persistence.rev_pct, P.rev_pct, "RelTol", 1e-12);
verifyEqual(tc, R.persistence.n_runs, P.n_runs);
verifyEqual(tc, R.persistence.max_run_revs, P.max_run_revs);
V = readtable(fullfile(here, "docs", "figures", "persistence_volume_120h_rev5.csv"));
verifyEqual(tc, R.persistence_volume.vol_pct, V.vol_pct, "RelTol", 1e-12);
verifyEqual(tc, R.persistence_volume.n_runs, V.n_runs);
verifyEqual(tc, R.persistence_volume.max_run_vols, V.max_run_vols);
end

function testHarmMetricMatchesBruteForce(tc)
% Corrupted radials at 0.5 deg, against every sample of the fixed grid.
d = tc.TestData;
w = d.fixed;
spec = w.scan_max.clock_spec;
lib = scan_event_library_rev1(spec, struct("seed", d.o.seed, "n_trials", 1, "duration_s", d.o.duration_s));
thr = [40 50];
R = mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, 0, 0, struct("radial_thr_dBi", thr));
C = scan_clock_rev1(spec, d.o.seed, 1, d.o.duration_s);
t = (0:0.002:d.o.duration_s)';
[st_n, st_a, ck] = scan_clock_state_rev1(C, t);
[~, ~, ck0] = scan_clock_state_rev1(C, 0);
pn = antenna_prepare_rev1(w.sys_nexrad.tx);
pa = antenna_prepare_rev1(w.sys_asr11.tx);
g21 = antenna_gain_rev1(pa, 180, 0, st_a) + antenna_gain_rev1(pn, 0, 0, st_n);
n_cut = numel(spec.rate_n);
key = (ck.vol - ck0.vol)*n_cut + ck.cut - ck0.cut + 1;
low = abs(spec.el_n(ck.cut) - 0.5) < 1e-6;
radial = floor(mod(st_n.mech_az, 360)) + 1;
for i = 1:numel(thr)
    hit = low & g21 >= thr(i);
    u = unique([key(hit), radial(hit)], "rows");
    n_by_key = accumarray(u(:, 1), 1, [max(key) 1]);
    keys_low = unique(key(low));
    expect = n_by_key(keys_low);
    got = R.radial_cuts.("n_radials_ge_" + thr(i));
    verifyEqual(tc, got, expect, sprintf("threshold %g dBi", thr(i)));
end
verifyGreaterThan(tc, R.radial_summary.cuts_with_hit_pct(1), 0);
end

function testRandomWalkDrift(tc)
% Correlated drift stays in the tolerance band and is autocorrelated.
d = tc.TestData;
spec = d.fixed.scan_max.clock_spec;
spec.jitter_mode = "random_walk";
spec.drift_step_a = 0.02;
C = scan_clock_rev1(spec, 1, 1, 3600);
r = C.r_rev;
verifyLessThanOrEqual(tc, max(abs(r - spec.rate_a)), spec.jit_a + 1e-12);
c = corrcoef(r(1:end - 1), r(2:end));
verifyGreaterThan(tc, c(1, 2), 0.9);
o = d.o;
o.asr_jitter_mode = "random_walk";
o.n_trials = 1;
o.duration_s = 120;
w = nexrad_asr11_mag_rev1(o);
verifyEqual(tc, height(w.scan_max.asr_rev), 26, "AbsTol", 1);   % 120 s at about 4.8 s per revolution
end

function testParallelEqualsSerial(tc)
assumeTrue(tc, ~isempty(ver("parallel")) && canUseParallelPool, "No parallel pool available.");
if isempty(gcp("nocreate"))
    parpool("Processes", 2);
end
o = tc.TestData.o;
o.parallel_flag = 1;
p = nexrad_asr11_mag_rev1(o);
verifyEqual(tc, p.scan_max.asr_rev, tc.TestData.fixed.scan_max.asr_rev);
verifyEqual(tc, p.scan_max.lobe_counts, tc.TestData.fixed.scan_max.lobe_counts);
verifyEqual(tc, p.ts_pct, tc.TestData.fixed.ts_pct);
end

function testInputErrors(tc)
d = tc.TestData;
w = d.fixed;
lib = scan_event_library_rev1(w.scan_max.clock_spec, struct("seed", 1, "n_trials", 1, "duration_s", 30));
verifyError(tc, @() mag_pair_timeseries_rev1(lib, w.sys_nexrad.tx, w.sys_asr11.tx, 0, 0, ...
    struct("persist_thr_dBi", 5)), "mag_pair_ts:thresholdBelowFloor");
verifyError(tc, @() nexrad_asr11_mag_rev1(struct("seed", 1, "ts_method", "fast")), "radar_coupling:badTsMethod");
verifyError(tc, @() scan_event_windows_rev1(scan_clock_rev1(w.scan_max.clock_spec, 1, 1, 10), 10, 0.002, ...
    struct("dt_coarse_s", 0.2)), "scan_events:coarseStep");
end

% =========================================================================
function check_against(tc, R, S, label)
% Every ATC revolution and NEXRAD volume of the fixed grid: same rows, same
% maxima above the floor, censored (NaN) at or below it.
fl = max(R.floor21_dBi, R.floor12_dBi);
verifyEqual(tc, R.asr_rev.complete, S.asr_rev.complete);
check_max(tc, R.asr_rev.max_g21_dBi, S.asr_rev.max_g21_dBi, fl, "revolution maxima " + label);
check_max(tc, R.asr_rev.max_g12_dBi, S.asr_rev.max_g12_dBi, fl, "revolution maxima g12 " + label);
check_max(tc, R.nexrad_volume.max_g21_dBi, S.nexrad_volume.max_g21_dBi, fl, "volume maxima " + label);
check_max(tc, R.nexrad_cut.max_g21_dBi, S.nexrad_cut.max_g21_dBi, fl, "cut maxima " + label);
verifyEqual(tc, R.persistence.rev_pct, S.persistence.rev_pct);
verifyEqual(tc, R.persistence.max_run_revs, S.persistence.max_run_revs);
end

function check_max(tc, got, ref, fl, label)
verifyEqual(tc, numel(got), numel(ref), label);
hi = ref > fl;
verifyEqual(tc, got(hi), ref(hi), label + " above the floor");
verifyTrue(tc, all(isnan(got(~hi)) | got(~hi) == ref(~hi)), label + " at or below the floor");
end
