function tests = test_mag_lookup_rev1
%TEST_MAG_LOOKUP_REV1 Tests for gain_dist_table_rev1 and mag_pair_cdf_rev1.
%   The reference is an independent code path: the random-mode Monte Carlo
%   of radar_mutual_gain_cdf_rev1 (random draws, prctile). An exact
%   percentile must lie inside the Monte Carlo order-statistic interval,
%   widened upward by the conservative grid bias R.bias_dB.
%
%   Run: results = runtests("test_mag_lookup_rev1");
tests = functiontests(localfunctions);
end

% =========================================================================
function setupOnce(tc)
o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, "nexrad_pattern", "nssl_measured"));
scan_n = nexrad_vcp_rev1(212, sails=3, avset_top_deg=6.4);   % Rev 5 worst mode
tc.TestData.sys_a = o.sys_asr11;
tc.TestData.sys_n = struct("name", "WSR-88D", "tx", o.sys_nexrad.tx, "rx", o.sys_nexrad.tx, "scan", scan_n);
grid = -0.5:0.05:2.5;
tc.TestData.grid = grid;
tc.TestData.Ta = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, ...
    struct("name", "ASR-11", "el_grid_deg", grid));
tc.TestData.Tn = gain_dist_table_rev1(o.sys_nexrad.tx, scan_n, ...
    struct("name", "WSR-88D", "mode_names", "VCP 212 SAILSx3 AVSET 6.4", "el_grid_deg", grid));
[scans, names] = nexrad_mode_scans_rev1("lowest");   % the 74 modes of the Rev 5 study
tc.TestData.Tn_all = gain_dist_table_rev1(o.sys_nexrad.tx, scans, ...
    struct("name", "WSR-88D", "mode_names", names, "el_grid_deg", -0.1:0.05:0.1));
tc.TestData.Ta_small = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, ...
    struct("name", "ASR-11", "el_grid_deg", -0.1:0.05:0.1));
tc.TestData.pct = [90 99 99.9 99.99 99.999];
tc.TestData.n_mc = 4e6;
end

% =========================================================================
function testMatchesMonteCarloOnGrid(tc)
% Takeoff angles on grid nodes: 0/0 (the Rev 5 baseline) and 2.0/1.7.
for e = [0 0; 2.0 1.7]'
    check_vs_mc(tc, e(1), e(2));
end
end

function testMatchesMonteCarloOffGrid(tc)
% Takeoff angles between grid nodes exercise the quantile interpolation.
check_vs_mc(tc, 1.73, 2.02);
end

function testEnvelopeBoundsEveryMode(tc)
d = tc.TestData;
R = mag_pair_cdf_rev1(d.Ta_small, d.Tn_all, 0, 0, struct("per_mode", true, "pct", d.pct));
env = R.pct{1, 3:2 + numel(d.pct)};
per_mode = reshape(R.mode_pct(1, :, :), [], numel(d.pct));
verifyEqual(tc, env, max(per_mode, [], 1), "Envelope must equal the largest mode percentile.");
verifyEqual(tc, numel(R.mode_pairs), 74);
end

function testSingleBinSliceInterpolates(tc)
% Regression: at high tilts (19.5 deg) the whole WSR-88D slice sits at the
% far-sidelobe floor, one bin, so its CDF has one distinct value. Off-node
% elevations must still interpolate.
d = tc.TestData;
R = mag_pair_cdf_rev1(d.Ta_small, d.Tn_all, 0.025, 0.075, struct("per_mode", true, "pct", d.pct));
verifyTrue(tc, all(isfinite(R.pct{1, 3:end})));
end

function testEnvelopeMatchesRev5ModeEnvelope(tc)
% Rev 5 Monte Carlo envelope over the same 74 modes (5e6 samples per mode,
% docs/figures/mode_envelope_rev5.csv): 20.156, 43.949, 45.492, 68.857 dBi.
% Allowed differences: the Rev 5 order-statistic half-widths at the
% baseline mode (convergence_rev5.csv: 0.15, 0.09, 0.002, 2.1 dB), plus the
% upward grid bias.
d = tc.TestData;
R = mag_pair_cdf_rev1(d.Ta_small, d.Tn_all, 0, 0, struct("pct", [99 99.9 99.99 99.999]));
got = R.pct{1, 3:6};
ref = [20.156 43.949 45.492 68.857];
tol = [0.15 0.09 0.002 2.1];
verifyLessThanOrEqual(tc, got, ref + tol + R.bias_dB);
verifyGreaterThanOrEqual(tc, got, ref - tol);
end

function testAvsetTerminationsCloseEnvelopeGap(tc)
% Review R3: an AVSET termination above 6.4 deg gives the cuts just above
% 6.4 deg more time than any of the 74 Rev 5 modes; near 7.25 deg the
% 74-mode envelope was 2.6 dB low at 99.9%. The default mode set now holds
% every termination.
[s_all, n_all] = nexrad_mode_scans_rev1();
[~, n_74] = nexrad_mode_scans_rev1("lowest");
verifyEqual(tc, numel(n_74), 74);
verifyTrue(tc, all(ismember(n_74, n_all)));
verifyTrue(tc, ismember("VCP 12 AVSET 8.0", n_all));
d = tc.TestData;
g = [7.2 7.25 7.3];
T_all = gain_dist_table_rev1(d.sys_n.tx, s_all, struct("mode_names", n_all, "el_grid_deg", g));
Ta = gain_dist_table_rev1(d.sys_a.tx, d.sys_a.scan, struct("name", "ASR-11", "el_grid_deg", -7.3:0.05:-7.2));
R_all = mag_pair_cdf_rev1(Ta, T_all, -7.25, 7.25, struct("pct", 99.9));
R_74 = mag_pair_cdf_rev1(Ta, T_all, -7.25, 7.25, struct("pct", 99.9, "rx_modes", n_74));
verifyGreaterThan(tc, R_all.pct.p99_9, R_74.pct.p99_9 + 1);
verifyFalse(tc, contains(R_all.worst_mode, "AVSET 6.4") || ~contains(R_all.worst_mode, "AVSET"));
end

function testNodeIsExactTable(tc)
% On a node the pair CDF uses the stored distributions directly, so
% "quantile" and "upper" agree.
d = tc.TestData;
a = mag_pair_cdf_rev1(d.Ta, d.Tn, 1.0, 0.5);
b = mag_pair_cdf_rev1(d.Ta, d.Tn, 1.0, 0.5, struct("interp", "upper"));
verifyEqual(tc, a.pct, b.pct);
end

function testUpperAndRoundingNeverBelowQuantile(tc)
% Review R4: "nearest" and plain rounding came out low. "upper" and the
% angle cells of el_round_deg take an upper envelope, so every percentile
% is at least the unrounded quantile-interpolated value, and every pair in
% a cell gets the same result.
d = tc.TestData;
e = (0.3:0.0137:1.7)';
q = mag_pair_cdf_rev1(d.Tn, d.Ta, e, 1.0 + 0*e);
for o = {struct("interp", "upper"), struct("el_round_deg", 0.01), struct("el_round_deg", 0.1), ...
        struct("el_round_deg", 0.1, "interp", "upper")}
    r = mag_pair_cdf_rev1(d.Tn, d.Ta, e, 1.0 + 0*e, o{1});
    verifyGreaterThanOrEqual(tc, r.pct{:, 3:end - 1}, q.pct{:, 3:end - 1} - 1e-9, jsonencode(o{1}));
    verifyEqual(tc, r.pct.el_tx_deg, e);
end
% Same cell, same result.
r = mag_pair_cdf_rev1(d.Tn, d.Ta, [0.712; 0.788; 0.75], [1; 1; 1], struct("el_round_deg", 0.1));
verifyEqual(tc, r.pct{1, 3:end}, r.pct{2, 3:end});
verifyError(tc, @() mag_pair_cdf_rev1(d.Tn, d.Ta, 0.5, 1, struct("interp", "nearest")), "mag_pair:nearestRemoved");
end

function testInterpolationError(tc)
% Off-node result against a table built with that elevation as a node.
d = tc.TestData;
e_off = [0.025 0.525 1.225 1.725];
g_x = unique([e_off - 0.025, e_off, e_off + 0.025]);
Ta_x = gain_dist_table_rev1(d.sys_a.tx, d.sys_a.scan, struct("el_grid_deg", g_x));
Tn_x = gain_dist_table_rev1(d.sys_n.tx, d.sys_n.scan, struct("el_grid_deg", g_x));
for e = e_off
    % Ta_x and Tn_x hold e as a node; d.Ta and d.Tn interpolate it.
    exact = mag_pair_cdf_rev1(Ta_x, Tn_x, e, e, struct("pct", d.pct));
    interp = mag_pair_cdf_rev1(d.Ta, d.Tn, e, e, struct("pct", d.pct));
    diff_dB = interp.pct{1, 3:2 + numel(d.pct)} - exact.pct{1, 3:2 + numel(d.pct)};
    verifyLessThanOrEqual(tc, abs(diff_dB(1:4)), 0.15, sprintf("Interpolation error at %.3f deg", e));
    verifyLessThanOrEqual(tc, abs(diff_dB(5)), 0.5, sprintf("99.999%% interpolation error at %.3f deg", e));
end
end

function testInterpolationWeightsPointTheRightWay(tc)
% Red team F2: swapped interpolation weights (a for 1 - a) passed every
% test but the pinned reference, because testInterpolationError checks
% midpoints only. Near a node the result must sit near that node: in the
% grid intervals where the node percentiles differ by 0.5 dB or more, a
% point one tenth of the way in is closer to the near node.
d = tc.TestData;
g = d.Tn.el_grid_deg;
q = mag_pair_cdf_rev1(d.Tn, d.Ta, g(:), 0*g(:));
v = q.pct{:, 3:2 + numel(q.pct_list)};
dv = abs(diff(v, 1, 1));
[~, o] = sort(max(dv, [], 2), "descend");
iv = o(1:min(5, numel(o)));
iv = iv(max(dv(iv, :), [], 2) >= 0.5);
verifyNotEmpty(tc, iv);
for i = iv'
    for a = [0.1 0.9]
        e = g(i) + a*(g(i + 1) - g(i));
        r = mag_pair_cdf_rev1(d.Tn, d.Ta, e, 0).pct{1, 3:2 + numel(q.pct_list)};
        big = dv(i, :) >= 0.5;
        near = (a < 0.5)*v(i, :) + (a > 0.5)*v(i + 1, :);
        far = (a < 0.5)*v(i + 1, :) + (a > 0.5)*v(i, :);
        verifyLessThan(tc, abs(r(big) - near(big)), abs(r(big) - far(big)), ...
            sprintf("interval %.2f to %.2f deg, a = %.1f", g(i), g(i + 1), a));
    end
end
end

function testEqualRatesWarn(tc)
d = tc.TestData;
verifyWarning(tc, @() mag_pair_cdf_rev1(d.Ta, d.Ta, 0, 0), "mag_pair:equalRates");
end

function testInputErrors(tc)
d = tc.TestData;
verifyError(tc, @() mag_pair_cdf_rev1(d.Ta, d.Tn, 3, 0), "mag_pair:outOfGrid");
verifyError(tc, @() mag_pair_cdf_rev1(d.Ta, d.Tn, [0 1], 0), "mag_pair:sizeMismatch");
verifyError(tc, @() mag_pair_cdf_rev1(d.Ta, d.Tn, 0, 0, struct("tx_modes", "none")), "mag_pair:unknownMode");
verifyError(tc, @() mag_pair_cdf_rev1(d.Ta, d.Tn, 0, 0, struct("bogus", 1)), "mag_pair:unknownOption");
Tb = d.Tn;
Tb.bin_dB = 0.1;
verifyError(tc, @() mag_pair_cdf_rev1(d.Ta, Tb, 0, 0), "mag_pair:binMismatch");
verifyError(tc, @() mag_pair_cdf_rev1(d.Ta, d.Tn, 0, 0, struct("conv_method", "bogus")), "mag_pair:badConv");
verifyError(tc, @() mag_pair_cdf_rev1(d.Ta, d.Tn, 0, 0, struct("el_round_deg", -1)), "mag_pair:badRound");
% Review R6: clear errors for repeated percentiles and an empty pair list;
% an angle within the 1e-12 tolerance past the grid edge runs.
verifyError(tc, @() mag_pair_cdf_rev1(d.Ta, d.Tn, 0, 0, struct("pct", [99 99])), "mag_pair:duplicatePct");
verifyError(tc, @() mag_pair_cdf_rev1(d.Ta, d.Tn, zeros(0, 1), zeros(0, 1)), "mag_pair:noPairs");
g = d.Tn.el_grid_deg(end);
R = mag_pair_cdf_rev1(d.Ta, d.Tn, 0, g + 1e-13);
verifyEqual(tc, R.pct{1, 3:end}, mag_pair_cdf_rev1(d.Ta, d.Tn, 0, g).pct{1, 3:end});
% Sector scans are allowed (plan step 5), but every mode of one table must
% share the sector geometry.
s = d.sys_a.scan;
s.sector_deg = 90;
verifyError(tc, @() gain_dist_table_rev1(d.sys_a.tx, {d.sys_a.scan, s}), "gain_dist:mixedSectors");
verifyError(tc, @() gain_dist_table_rev1(d.sys_a.tx, d.sys_a.scan, struct("bogus", 1)), "gain_dist:unknownOption");
end

function testWeightsSumToOne(tc)
d = tc.TestData;
verifyEqual(tc, sum(d.Tn_all.W, 2), ones(74, 1), "AbsTol", 1e-12);
verifyEqual(tc, double(reshape(sum(d.Tn.P_tx, 1), numel(d.grid), [])), ...
    ones(numel(d.grid), numel(d.Tn.tilt_deg)), "AbsTol", 1e-5);
end

function testParallelEqualsSerial(tc)
assumeTrue(tc, ~isempty(ver("parallel")) && canUseParallelPool, "No parallel pool available.");
if isempty(gcp("nocreate"))
    parpool("Processes", 2);
end
parfevalOnAll(@addpath, 0, fileparts(which("gain_dist_table_rev1")));
d = tc.TestData;
o = struct("name", "WSR-88D", "el_grid_deg", -0.1:0.05:0.1);
T0 = gain_dist_table_rev1(d.sys_n.tx, d.sys_n.scan, o);
o.parallel_flag = 1;
T1 = gain_dist_table_rev1(d.sys_n.tx, d.sys_n.scan, o);
verifyEqual(tc, T1.P_tx, T0.P_tx);
e = [-0.1 0 0.07 0.1 0.07 0.02 0.03 0.02];      % repeats: groups and intervals with several pairs
R0 = mag_pair_cdf_rev1(d.Ta_small, T0, e, fliplr(e));
R1 = mag_pair_cdf_rev1(d.Ta_small, T0, e, fliplr(e), struct("parallel_flag", 1));
verifyEqual(tc, R1.pct, R0.pct);
end

function testFftEqualsDirect(tc)
% Plan T4: the batched FFT convolution gives the percentiles of the direct
% one exactly, and the envelope curve to round-off, in both directions.
d = tc.TestData;
[et, er] = batch_pairs();
o = struct("per_mode", true, "return_curve", true);
od = o;
od.conv_method = "direct";
for k = 1:2
    if k == 1
        [A, B] = deal(d.Tn, d.Ta);
    else
        [A, B] = deal(d.Ta, d.Tn);
    end
    R1 = mag_pair_cdf_rev1(A, B, et, er, o);
    R0 = mag_pair_cdf_rev1(A, B, et, er, od);
    verifyEqual(tc, R1.pct, R0.pct);
    verifyEqual(tc, R1.mode_pct, R0.mode_pct);
    verifyEqual(tc, R1.worst_mode, R0.worst_mode);
    for i = 1:numel(et)
        verifyEqual(tc, R1.curve{i}.gain_dBi, R0.curve{i}.gain_dBi);
        verifyEqual(tc, R1.curve{i}.exceed_pct, R0.curve{i}.exceed_pct, "AbsTol", 1e-9);
    end
end
end

function testGroupedEqualsOneByOne(tc)
% Plan T5: grouping pairs by elevation (and sharing slices, FFTs and the
% interpolation setup) does not change any pair's result.
d = tc.TestData;
[et, er] = batch_pairs();
R = mag_pair_cdf_rev1(d.Tn, d.Ta, et, er);
for i = 1:numel(et)
    Ri = mag_pair_cdf_rev1(d.Tn, d.Ta, et(i), er(i), struct("conv_method", "direct"));
    verifyEqual(tc, R.pct(i, :), Ri.pct, sprintf("pair %d", i));
end
end

function testRoundingCellsMatchCellEnds(tc)
% Angles on the cell grid are used as they are; off-grid angles get the
% larger of their two cell ends (no grid node lies inside a 0.01 deg cell
% of the 0.05 deg table), and the table reports the angles given.
d = tc.TestData;
[et, er] = batch_pairs();
et = round(et/0.01)*0.01;
er = round(er/0.01)*0.01;
R0 = mag_pair_cdf_rev1(d.Tn, d.Ta, et, er);
R1 = mag_pair_cdf_rev1(d.Tn, d.Ta, et, er, struct("el_round_deg", 0.01));
verifyEqual(tc, R1.pct{:, 3:end}, R0.pct{:, 3:end}, "AbsTol", 1e-9);
e2 = et + 0.003;                                 % off the 0.01 grid, inside the cell [et, et + 0.01]
R2 = mag_pair_cdf_rev1(d.Tn, d.Ta, e2, er, struct("el_round_deg", 0.01));
Rhi = mag_pair_cdf_rev1(d.Tn, d.Ta, min(et + 0.01, 2.5), er);
verifyGreaterThanOrEqual(tc, R2.pct{:, 3:end - 1}, max(R0.pct{:, 3:end - 1}, Rhi.pct{:, 3:end - 1}) - 1e-9);
verifyEqual(tc, R2.pct.el_tx_deg, e2);
end

% =========================================================================
function [et, er] = batch_pairs()
% 24 pairs inside the -0.5 to 2.5 deg grid: repeated angles on both sides,
% several in one table interval, and nodes.
rs = RandStream("mrg32k3a", "Seed", 7);
et = [round(-0.4 + 2.8*rand(rs, 12, 1), 3); 0.5; 0.5; 0.5; 1.0];
er = [round(-0.4 + 2.8*rand(rs, 12, 1), 3); 0.2; 0.21; 0.2; 2.5];
et = [et; et(1:4) + 0.01; et(5:8)];
er = [er; er(1:4); er(5:8) + 0.013];
end

% =========================================================================
function check_vs_mc(tc, el_tx, el_rx)
d = tc.TestData;
R = mag_pair_cdf_rev1(d.Ta, d.Tn, el_tx, el_rx, struct("pct", d.pct));
got = R.pct{1, 3:2 + numel(d.pct)};
mc = radar_mutual_gain_cdf_rev1(d.sys_a, d.sys_n, struct("seed", 7, "mode", "random", ...
    "n_samples", d.n_mc, "el12_deg", el_tx, "el21_deg", el_rx));
g = sort(mc.g12_dBi);
n = numel(g);
for j = 1:numel(d.pct)
    p = d.pct(j)/100;
    h = 3.5*sqrt(n*p*(1 - p));                  % order-statistic band, about 99.95% two-sided
    lo = g(max(1, floor(n*p - h)));
    hi = g(min(n, ceil(n*p + h)));
    verifyGreaterThanOrEqual(tc, got(j), lo - 1e-9, ...
        sprintf("%.3f%% at %.2f/%.2f deg: %.3f below MC band [%.3f %.3f]", d.pct(j), el_tx, el_rx, got(j), lo, hi));
    verifyLessThanOrEqual(tc, got(j), hi + R.bias_dB + 1e-9, ...
        sprintf("%.3f%% at %.2f/%.2f deg: %.3f above MC band [%.3f %.3f] + bias", d.pct(j), el_tx, el_rx, got(j), lo, hi));
end
end
