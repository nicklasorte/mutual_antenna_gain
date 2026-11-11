function tests = test_sectors_rev1
%TEST_SECTORS_REV1 Tests for plan step 5: sector scans in both models.
%   Covers sector_sweep_rev1 (shared sweep law), the engine's
%   bidirectional sweep with turnaround (random and time series), sector
%   tables and bearings in the two-part model, and the wrapper options
%   nexrad_sector_center_deg, nexrad_sweep_mode, nexrad_turnaround_s,
%   atc_sector_deg and atc_sector_center_deg. Plan section 6.4.
%
%   Run: results = runtests("test_sectors_rev1");
tests = functiontests(localfunctions);
end

% =========================================================================
function setupOnce(tc)
tc.TestData.s1 = struct("tx", struct("g_max_dBi", 40), "scan", struct("type", "mech_rotating", ...
    "rate_deg_s", 30, "sector_deg", 60, "sector_center_deg", 0, "sweep_mode", "bidirectional", ...
    "turnaround_s", 1));
tc.TestData.s2 = struct("tx", struct("g_max_dBi", 34), "scan", struct("type", "mech_rotating", "rate_deg_s", 47.3));
tc.TestData.pct = [90 99 99.9 99.99];
tc.TestData.n_mc = 4e6;
end

% =========================================================================
function testRestartMatchesRotationFormula(tc)
% Restart law = the rev 1 engine and wrapper formula, bit for bit.
t = (0:0.013:50)';
r = 12.524;
u = sector_sweep_rev1(t, r, 90, 33.3);
verifyEqual(tc, u, mod(33.3 + r.*t, 90));
end

function testBidirectionalStaysInSectorAndDwells(tc)
% Bidirectional law: stays in [0, sector], moves at the rate, and spends
% turnaround/(sector/rate + turnaround) of the time at the edges.
dt = 1e-4;
t = (0:dt:198 - dt)';                           % 33 whole 6 s cycles (2 x (2 s sweep + 1 s dwell))
sec = 60; r = 30; turn = 1;
u = sector_sweep_rev1(t, r, sec, 0, "bidirectional", turn);
verifyGreaterThanOrEqual(tc, min(u), 0);
verifyLessThanOrEqual(tc, max(u), sec);
moving = u > 0 & u < sec;
verifyEqual(tc, max(abs(diff(u(moving & [moving(2:end); false])))), r*dt, "RelTol", 1e-6);
f = turn/(sec/r + turn);
verifyEqual(tc, mean(~moving), f, "AbsTol", 2e-3);
verifyEqual(tc, mean(u == sec), f/2, "AbsTol", 2e-3);
end

function testEngineRandomMatchesExact(tc)
% Bidirectional sector with turnaround: engine random mode vs the exact
% CDF, path at the sector center and at the sector edge (where the
% turnaround dwell parks the main beam on the path).
d = tc.TestData;
check_vs_mc(tc, d.s1, d.s2, "bidirectional, path at center");
s1 = d.s1;
s1.scan.sector_center_deg = 30;                 % edge on the path
check_vs_mc(tc, s1, d.s2, "bidirectional, path at edge");
s1.scan.sweep_mode = "restart";
s1.scan.turnaround_s = 0;
check_vs_mc(tc, s1, d.s2, "restart, path at edge");
end

function testEdgeDwellRaisesTailAtEdge(tc)
% With the path at the edge, the dwell adds main-beam time: the 99%
% value must rise against the same sector without turnaround.
d = tc.TestData;
s1 = d.s1;
s1.scan.sector_center_deg = 30;
a = mag_exact_cdf_rev1(s1, d.s2, struct("pct", 99));
s1.scan.turnaround_s = 0;
b = mag_exact_cdf_rev1(s1, d.s2, struct("pct", 99));
verifyGreaterThan(tc, a.pct12.gain_dBi - b.pct12.gain_dBi, 3);
end

function testEngineTimeSeriesEdgeShare(tc)
% Time series with the shared sweep law: the beam sits at the two edges
% for the expected share of the time, and the long-run time-percentile
% CDF agrees with the exact CDF at 90% and 99%.
d = tc.TestData;
o = struct("seed", 5, "mode", "time_series", "duration_s", 600, "n_trials", 4, "dt_s", 0.002);
r = radar_mutual_gain_cdf_rev1(d.s1, d.s2, o);
x = mag_exact_cdf_rev1(d.s1, d.s2, struct("pct", [90 99]));
q = prctile(r.g12_dBi, [90 99]);
verifyLessThanOrEqual(tc, abs(q(:) - x.pct12.gain_dBi), 0.5);
end

function testBidirectionalClockSweepsWholeSector(tc)
% Review R1: the wrapper's scan clock restarted the back-and-forth sweep at
% a fixed point every cut, so part of the sector was never visited. Now
% the beam reaches both edges and dwells at each for the random-mode share.
o = struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, "cdf_method", "exact", ...
    "nexrad_sector_deg", 30, "nexrad_sweep_mode", "bidirectional", "nexrad_turnaround_s", 0.5);
w = nexrad_asr11_mag_rev1(o);
spec = w.scan_max.clock_spec;
pass = spec.sec_n./abs(spec.rate_n) + spec.turn_n;
share = spec.turn_n*numel(pass)/2/sum(pass);     % time at one edge, per edge
t = (0:0.005:600)';
for j = 1:4
    C = scan_clock_rev1(spec, 1, j, 600);
    st = scan_clock_state_rev1(C, t);
    off = wrapTo180(st.mech_az - spec.center_n);
    verifyLessThanOrEqual(tc, min(off), -15 + 1e-6, sprintf("trial %d near edge", j));
    verifyGreaterThanOrEqual(tc, max(off), 15 - 1e-6, sprintf("trial %d far edge", j));
    verifyEqual(tc, mean(off <= -15 + 1e-6), share, "RelTol", 0.15);
    verifyEqual(tc, mean(off >= 15 - 1e-6), share, "RelTol", 0.15);
    h = histcounts(off(abs(off) < 14), -14:7:14);
    verifyEqual(tc, h/mean(h), ones(1, 4), "AbsTol", 0.1);
end
end

function testEngineBidirectionalOddPassesMatchExact(tc)
% Review R1 in the engine: with one pass per step (n_rev = 1) the time
% series left part of the sector unvisited. With the path 10 deg off the
% sector center the old code missed the exact CDF by 3.4 dB at the median
% and 10 dB at 99%; it must now agree.
o = struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, "cdf_method", "exact", ...
    "nexrad_sector_deg", 30, "nexrad_sector_center_deg", 10, "nexrad_sweep_mode", "bidirectional", ...
    "nexrad_turnaround_s", 0.5);
w = nexrad_asr11_mag_rev1(o);
r = radar_mutual_gain_cdf_rev1(w.sys_nexrad, w.sys_asr11, struct("seed", 3, "mode", "time_series", ...
    "duration_s", 900, "n_trials", 4, "pct", [50 90 99]));
x = mag_exact_cdf_rev1(w.sys_nexrad, w.sys_asr11, struct("pct", [50 90 99]));
verifyLessThanOrEqual(tc, abs(r.pct12.gain_dBi - x.pct12.gain_dBi), 0.5);
end

function testTwoPartSectorMatchesExact(tc)
% A sector table stores no slices; mag_pair_cdf_rev1 makes them from the
% bearing. Receiver = WSR-88D sweeping a 90 deg bidirectional sector with
% 0.5 s turnaround; brg_rx_deg = 20 means the transmitter sits 20 deg
% clockwise from the sector center, so the exact engine's absolute center
% is 180 - 20 (the transmitter is at azimuth 180 from the receiver).
o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, "nexrad_pattern", "nssl_measured"));
scan = nexrad_vcp_rev1(212, sails=3, avset_top_deg=6.4);
scan.sector_deg = 90;
scan.sweep_mode = "bidirectional";
scan.turnaround_s = 0.5;
Tn = gain_dist_table_rev1(o.sys_nexrad.tx, scan, struct("name", "WSR-88D sector", "el_grid_deg", [-1 3]));
verifyEmpty(tc, Tn.P_tx);
verifyGreaterThan(tc, sum(Tn.W_edge(:)), 0);
% The transmitter's takeoff angle is a grid node, so its slice is exact too.
Ta = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, struct("name", "ASR-11", "el_grid_deg", [1.7 2.0]));
pct = [90 99 99.9 99.99];
R = mag_pair_cdf_rev1(Ta, Tn, [2.0; 2.0], [1.7; 1.7], struct("pct", pct, "brg_rx_deg", [20; 0]));
for k = 1:2
    b = [20 0];
    sn = struct("tx", o.sys_nexrad.tx, "scan", scan);
    sn.scan.sector_center_deg = 180 - b(k);
    x = mag_exact_cdf_rev1(o.sys_asr11, sn, struct("el12_deg", 2.0, "el21_deg", 1.7, "pct", pct));
    verifyEqual(tc, R.pct{k, 3:2 + numel(pct)}, x.pct12.gain_dBi', "AbsTol", 1e-9, sprintf("bearing %g", b(k)));
end
% Path outside the 90 deg sector (70 deg from its center): the NEXRAD
% main beam never points at the transmitter, so the tail drops sharply.
R2 = mag_pair_cdf_rev1(Ta, Tn, 2.0, 1.7, struct("pct", 99.99, "brg_rx_deg", 70));
verifyGreaterThan(tc, R.pct.p99_99(2) - R2.pct.p99_99, 20);
end

function testWrapperSectorReproducesRev5(tc)
% Plan 6.4 test 3: the Rev 5 sensitivity row "NEXRAD sector 10 deg"
% (docs/figures/sensitivity_sweep_rev5.csv) with the new center option at
% its default, bit for bit.
o = struct("seed", 20261105, "vcp", 212, "sails", 3, "avset_top_deg", 6.4, ...
    "n_samples", 5e6, "duration_s", 3600, "n_trials", 4, "nexrad_sector_deg", 10, "nexrad_pattern", "nssl_measured", "cdf_method", "monte_carlo");
r = nexrad_asr11_mag_rev1(o);
q = r.cdf.pct21;
got = arrayfun(@(p) q.gain_dBi(abs(q.pct - p) < 1e-9), [99 99.9 99.99 99.999]);
ref = [45.4154465467263, 57.8777634115202, 73.939273560345, 75.8138178986038];
verifyEqual(tc, got, ref, "RelTol", 1e-10);
end

function testWrapperSectorCenterAndExact(tc)
% A sector aimed away from the ATC radar cannot reach the main-beam
% plateau; the exact CDF agrees with Monte Carlo off center too.
base = struct("seed", 7, "vcp", 212, "sails", 3, "avset_top_deg", 6.4, "n_samples", 4e6, ...
    "duration_s", 120, "n_trials", 1, "nexrad_sector_deg", 30);
on = nexrad_asr11_mag_rev1(base);
o = base;
o.nexrad_sector_center_deg = 120;
off = nexrad_asr11_mag_rev1(o);
verifyGreaterThan(tc, pick(on, 99.99) - pick(off, 99.99), 20);
o.nexrad_sector_center_deg = 12;
o.nexrad_sweep_mode = "bidirectional";
o.nexrad_turnaround_s = 0.3;
m = nexrad_asr11_mag_rev1(o);
o.cdf_method = "exact";
x = nexrad_asr11_mag_rev1(o);
for p = [90 99 99.9]
    j = m.cdf.pct21.pct == p;
    lo = m.cdf.pct21.ci95_lo_dBi(j) - 0.02;
    hi = m.cdf.pct21.ci95_hi_dBi(j) + x.cdf.bias_dB + 0.02;
    verifyTrue(tc, pick(x, p) >= lo && pick(x, p) <= hi, sprintf("%g%%: exact %.3f, MC [%.3f %.3f]", ...
        p, pick(x, p), lo, hi));
end
end

function testWrapperAtcSector(tc)
o = struct("seed", 4, "vcp", 212, "n_samples", 4e6, "duration_s", 120, "n_trials", 1, ...
    "atc_sector_deg", 90, "atc_sector_center_deg", 10);
m = nexrad_asr11_mag_rev1(o);
o.cdf_method = "exact";
x = nexrad_asr11_mag_rev1(o);
for p = [90 99 99.9]
    j = m.cdf.pct21.pct == p;
    lo = m.cdf.pct21.ci95_lo_dBi(j) - 0.02;
    hi = m.cdf.pct21.ci95_hi_dBi(j) + x.cdf.bias_dB + 0.02;
    verifyTrue(tc, pick(x, p) >= lo && pick(x, p) <= hi, sprintf("%g%%: exact %.3f, MC [%.3f %.3f]", ...
        p, pick(x, p), lo, hi));
end
% One row per 90 deg sweep: 1.2 s at 75 deg/s (jitter 0.5 deg/s).
n_rows = height(m.scan_max.asr_rev);
verifyEqual(tc, n_rows, 120/1.2, "RelTol", 0.02);
% Lobe-pair shares from the patterns agree with the time series (120 s,
% so a few percentage points of sampling noise).
verifyLessThanOrEqual(tc, max(abs(m.lobe_pairs.time_pct - m.lobe_pairs.time_pct_ts)), 3);
end

function testInputErrors(tc)
s1 = tc.TestData.s1;
s1.scan.sector_deg = 360;
verifyError(tc, @() radar_mutual_gain_cdf_rev1(s1, tc.TestData.s2, struct("seed", 1, "mode", "random", ...
    "n_samples", 10)), "radar_coupling:badSweepMode");
verifyError(tc, @() sector_sweep_rev1(1, 30, 360, 0, "bidirectional"), "sector_sweep:bidirectional360");
verifyError(tc, @() nexrad_asr11_mag_rev1(struct("seed", 1, "atc_sector_deg", 90, ...
    "asr_jitter_mode", "per_rev")), "radar_coupling:sectorJitter");
verifyError(tc, @() nexrad_asr11_mag_rev1(struct("seed", 1, "nexrad_turnaround_s", 1)), ...
    "radar_coupling:badTurnaround");
end

% =========================================================================
function v = pick(r, p)
v = r.cdf.pct21.gain_dBi(abs(r.cdf.pct21.pct - p) < 1e-9);
end

function check_vs_mc(tc, s1, s2, label)
d = tc.TestData;
w = warning("off", "radar_coupling:commensurateRates");
cleanup = onCleanup(@() warning(w));
x = mag_exact_cdf_rev1(s1, s2, struct("pct", d.pct));
mc = radar_mutual_gain_cdf_rev1(s1, s2, struct("seed", 9, "mode", "random", "n_samples", d.n_mc));
g = sort(mc.g12_dBi);
n = numel(g);
for j = 1:numel(d.pct)
    p = d.pct(j)/100;
    h = 3.5*sqrt(n*p*(1 - p));
    lo = g(max(1, floor(n*p - h)));
    hi = g(min(n, ceil(n*p + h)));
    got = x.pct12.gain_dBi(j);
    verifyTrue(tc, got >= lo - 1e-9 && got <= hi + x.bias_dB + 1e-9, ...
        sprintf("%s at %.2f%%: exact %.3f, MC [%.3f %.3f]", label, d.pct(j), got, lo, hi));
end
end
