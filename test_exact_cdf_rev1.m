function tests = test_exact_cdf_rev1
%TEST_EXACT_CDF_REV1 Tests for plan step 4: exact CDFs in the original model.
%   Covers mag_exact_cdf_rev1 (O1), the wrapper options cdf_method,
%   el_nexrad_deg and el_atc_deg (O1, O5), and mode summaries from a
%   distribution for the 74-mode envelope (O2). References are independent
%   code paths: the random-mode Monte Carlo of radar_mutual_gain_cdf_rev1
%   (order-statistic bands), the two-part model, and the Rev 5 envelope.
%
%   Run: results = runtests("test_exact_cdf_rev1");
tests = functiontests(localfunctions);
end

% =========================================================================
function setupOnce(tc)
o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, "nexrad_pattern", "nssl_measured"));
scan_n = nexrad_vcp_rev1(212, sails=3, avset_top_deg=6.4);
tc.TestData.sys_a = o.sys_asr11;
tc.TestData.sys_n = struct("name", "WSR-88D", "tx", o.sys_nexrad.tx, "rx", o.sys_nexrad.tx, "scan", scan_n);
tc.TestData.pct = [90 99 99.9 99.99 99.999];
tc.TestData.n_mc = 4e6;
tc.TestData.small = struct("seed", 3, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, ...
    "vcp", 212, "sails", 3, "avset_top_deg", 6.4);
end

% =========================================================================
function testEngineMatchesMonteCarloStatgain(tc)
s1 = struct("tx", struct("g_max_dBi", 40), "scan", struct("type", "mech_rotating", "rate_deg_s", 30));
s2 = struct("tx", struct("g_max_dBi", 34), "scan", struct("type", "mech_rotating", "rate_deg_s", 47.3));
check_vs_mc(tc, s1, s2, 0, 0, "statgain");
check_vs_mc(tc, s1, s2, 3, -1, "statgain, raised path");
end

function testEngineMatchesMonteCarloNexradAsr(tc)
d = tc.TestData;
check_vs_mc(tc, d.sys_n, d.sys_a, 0, 0, "NEXRAD-ASR-11 0/0");
check_vs_mc(tc, d.sys_n, d.sys_a, 1.7, 2.0, "NEXRAD-ASR-11 1.7/2.0");
end

function testEngineMatchesMonteCarloSector(tc)
% Radar 1 sweeps a 60 deg sector: centered on the path, then aimed away.
s1 = struct("tx", struct("g_max_dBi", 40), "scan", struct("type", "mech_rotating", ...
    "rate_deg_s", 30, "sector_deg", 60, "sector_center_deg", 0));
s2 = struct("tx", struct("g_max_dBi", 34), "scan", struct("type", "mech_rotating", "rate_deg_s", 47.3));
check_vs_mc(tc, s1, s2, 0, 0, "sector on path");
s1.scan.sector_center_deg = 100;
check_vs_mc(tc, s1, s2, 0, 0, "sector aimed away");
end

function testMatchesTwoPartModel(tc)
% Same slices, same convolution: on grid nodes the two exact methods agree.
d = tc.TestData;
r = mag_exact_cdf_rev1(d.sys_a, d.sys_n, struct("el12_deg", 2.0, "el21_deg", 1.7, "pct", d.pct));
grid = 1.5:0.05:2.5;
Ta = gain_dist_table_rev1(d.sys_a.tx, d.sys_a.scan, struct("name", "ASR-11", "el_grid_deg", grid));
Tn = gain_dist_table_rev1(d.sys_n.tx, d.sys_n.scan, struct("name", "WSR-88D", "el_grid_deg", grid));
R = mag_pair_cdf_rev1(Ta, Tn, 2.0, 1.7, struct("pct", d.pct));
verifyEqual(tc, r.pct12.gain_dBi', R.pct{1, 3:2 + numel(d.pct)}, "AbsTol", 1e-9);
end

function testWrapperExactMatchesMonteCarlo(tc)
d = tc.TestData;
o = d.small;
o.n_samples = d.n_mc;
m = nexrad_asr11_mag_rev1(o);
o.cdf_method = "exact";
x = nexrad_asr11_mag_rev1(o);
verifyEmpty(tc, x.cdf.g21_dBi);
verifyEqual(tc, x.cdf.method, "exact");
for j = find(m.cdf.pct21.pct >= 90)'
    lo = m.cdf.pct21.ci95_lo_dBi(j) - 0.02;
    hi = m.cdf.pct21.ci95_hi_dBi(j) + x.cdf.bias_dB + 0.02;
    g = x.cdf.pct21.gain_dBi(j);
    verifyTrue(tc, g >= lo && g <= hi, sprintf("%.3f%%: exact %.3f, MC band [%.3f %.3f]", ...
        m.cdf.pct21.pct(j), g, lo, hi));
end
% Reference-level time shares. Exact gains are rounded up by at most one
% bin per antenna, so the exact share at level v lies between the true
% share at v and the true share at v - bias_dB. Near a plateau (NEXRAD
% main beam + ASR back lobe, 45.49 dBi) that is a factor of about four, so
% the check uses the Monte Carlo share at both levels, with binomial noise
% (4 sigma).
v = m.ref_levels.mag_dBi;
p_lo = arrayfun(@(u) mean(m.cdf.g21_dBi >= u), v);
p_hi = arrayfun(@(u) mean(m.cdf.g21_dBi >= u - x.cdf.bias_dB), v);
noise = @(p) 4*sqrt(max(p, 1e-6).*(1 - p)/d.n_mc);
got = x.ref_levels.time_exceeded_pct/100;
verifyGreaterThanOrEqual(tc, got, p_lo - noise(p_lo));
verifyLessThanOrEqual(tc, got, p_hi + noise(p_hi));
% The time series does not depend on the CDF method.
verifyEqual(tc, x.scan_max.nexrad_volume, m.scan_max.nexrad_volume);
end

function testWrapperTakeoffAnglesDefaultToPath(tc)
o = tc.TestData.small;
o.el_path_deg = 0.7;
a = nexrad_asr11_mag_rev1(o);
o = rmfield(o, "el_path_deg");
o.el_nexrad_deg = 0.7;
o.el_atc_deg = 0.7;
b = nexrad_asr11_mag_rev1(o);
verifyEqual(tc, b.cdf.pct21, a.cdf.pct21);
verifyEqual(tc, b.scan_max.asr_rev, a.scan_max.asr_rev);
verifyEqual(tc, b.ref_levels, a.ref_levels);
end

function testWrapperAsymmetricTakeoffMatchesEngine(tc)
% Wrapper radar 1 = NEXRAD, radar 2 = ATC radar: el12 = takeoff at the
% NEXRAD, el21 = takeoff at the ATC radar. Same seed, same samples.
o = tc.TestData.small;
o.n_samples = 2e5;
o.el_nexrad_deg = 1.7;
o.el_atc_deg = 2.0;
w = nexrad_asr11_mag_rev1(o);
e = radar_mutual_gain_cdf_rev1(w.sys_nexrad, w.sys_asr11, struct("seed", o.seed, "mode", "random", ...
    "n_samples", o.n_samples, "el12_deg", 1.7, "el21_deg", 2.0, "pct", w.opts.pct));
verifyEqual(tc, w.cdf.g21_dBi, e.g21_dBi);
% Exact wrapper at the same angles agrees with the exact engine.
o.cdf_method = "exact";
x = nexrad_asr11_mag_rev1(o);
r = mag_exact_cdf_rev1(w.sys_nexrad, w.sys_asr11, struct("el12_deg", 1.7, "el21_deg", 2.0, "pct", w.opts.pct));
verifyEqual(tc, x.cdf.pct21, r.pct21);
end

function testModeSummaryFromDistribution(tc)
d = tc.TestData;
r = mag_exact_cdf_rev1(d.sys_n, d.sys_a, struct("pct", d.pct));
s1 = mag_mode_summary_rev1(r.dist21, d.pct);
s2 = mag_mode_summary_rev1(r, d.pct);
verifyEqual(tc, s2, s1);
verifyEqual(tc, s1.q_dBi, r.pct21.gain_dBi);
verifyEqual(tc, s1.n, Inf);
verifyTrue(tc, all(diff(s1.exceed) <= 1e-12), "exceedance must not increase");
end

function testExactEnvelopeMatchesRev5(tc)
% O2: the 74-mode envelope from exact CDFs against the Rev 5 Monte Carlo
% envelope (docs/figures/mode_envelope_rev5.csv: 20.156, 43.949, 45.492,
% 68.857 dBi), within the Rev 5 order-statistic half-widths (0.15, 0.09,
% 0.002, 2.1 dB) plus the upward grid bias.
c = vcp_mode_cases_rev1("lowest");      % the 74 modes of the Rev 5 study
c = c(2:end);
pct = [99 99.9 99.99 99.999];
summ = cell(1, numel(c));
names = strings(1, numel(c));
for i = 1:numel(c)
    o = struct("seed", 1, "n_samples", 1e3, "duration_s", 30, "n_trials", 1, "cdf_method", "exact", ...
        "nexrad_pattern", "nssl_measured");   % Rev 5 pattern
    f = setdiff(fieldnames(c{i}), "label");
    for k = 1:numel(f)
        o.(f{k}) = c{i}.(f{k});
    end
    r = nexrad_asr11_mag_rev1(o);
    summ{i} = mag_mode_summary_rev1(r.cdf, pct);
    names(i) = string(c{i}.label);
end
env = mag_mode_envelope_rev1(summ, names, pct);
got = env.pct.envelope_dBi';
ref = [20.156 43.949 45.492 68.857];
tol = [0.15 0.09 0.002 2.1];
verifyLessThanOrEqual(tc, got, ref + tol + 0.1 + 1e-9);
verifyGreaterThanOrEqual(tc, got, ref - tol - 1e-9);
end

function testInputErrors(tc)
d = tc.TestData;
arr = struct("tx", struct("type", "array_planar"), "scan", struct("type", "array_multiface"));
verifyError(tc, @() mag_exact_cdf_rev1(arr, d.sys_a), "mag_exact:scanType");
verifyError(tc, @() mag_exact_cdf_rev1(d.sys_n, d.sys_a, struct("bogus", 1)), "mag_exact:unknownOption");
verifyError(tc, @() nexrad_asr11_mag_rev1(struct("seed", 1, "cdf_method", "fast")), "radar_coupling:badCdfMethod");
end

% =========================================================================
function check_vs_mc(tc, s1, s2, el12, el21, label)
d = tc.TestData;
w = warning("off", "radar_coupling:commensurateRates");
cleanup = onCleanup(@() warning(w));
r = mag_exact_cdf_rev1(s1, s2, struct("el12_deg", el12, "el21_deg", el21, "pct", d.pct));
mc = radar_mutual_gain_cdf_rev1(s1, s2, struct("seed", 9, "mode", "random", ...
    "n_samples", d.n_mc, "el12_deg", el12, "el21_deg", el21));
for dirn = ["12", "21"]
    g = sort(mc.("g" + dirn + "_dBi"));
    got = r.("pct" + dirn).gain_dBi;
    n = numel(g);
    for j = 1:numel(d.pct)
        p = d.pct(j)/100;
        h = 3.5*sqrt(n*p*(1 - p));
        lo = g(max(1, floor(n*p - h)));
        hi = g(min(n, ceil(n*p + h)));
        verifyTrue(tc, got(j) >= lo - 1e-9 && got(j) <= hi + r.bias_dB + 1e-9, ...
            sprintf("%s, g%s at %.3f%%: exact %.3f, MC band [%.3f %.3f]", label, dirn, d.pct(j), got(j), lo, hi));
    end
end
end
