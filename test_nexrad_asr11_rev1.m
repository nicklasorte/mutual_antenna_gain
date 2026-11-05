function tests = test_nexrad_asr11_rev1
%TEST_NEXRAD_ASR11_REV1 Tests for nexrad_vcp_rev1, the offaxis_table antenna
%   type and nexrad_asr11_mag_rev1. Run with runtests("test_nexrad_asr11_rev1").
%
%   Reference values:
%     Level II VCP sequences decoded with MetPy 1.7 (see nexrad_vcp_rev1):
%       KLOT20240507_171514 VCP 212 SAILSx1, sum(360/rate) = 325.068 s
%       KOAX20240426_150526 VCP 212 SAILSx1, 294.275 s
%       KHGX20240708_013308 VCP 212 SAILSx1, 309.691 s
%       KTBW20241009_035113 VCP 112 SAILSx1, 382.272 s
%     ROC VCP 112 Theory and ConOps (2019) Table 3: 20 scans without SAILS,
%       23 with SAILS, 15 with AVSET at 6.4 deg.
tests = functiontests(localfunctions);
end

% =========================================================================
% nexrad_vcp_rev1
% =========================================================================
function testVcp212SailsMatchesKLOT(tc)
s = nexrad_vcp_rev1(212, sails=1);
el = [0.5 0.5 0.9 0.9 1.3 1.3 0.5 0.5 1.8 2.4 3.1 4.0 5.1 6.4 8.0 10.0 12.5 15.6 19.5]';
r = [21.149 12.524 21.149 12.524 23.027 12.524 21.149 12.524 26.389 27.334 28.224 ...
    26.400 26.400 26.400 28.422 28.411 28.740 28.740 28.740]';
verifyEqual(tc, s.tilt_deg, el, "AbsTol", 1e-12);
verifyEqual(tc, s.rate_deg_s, r, "AbsTol", 1e-12);
verifyEqual(tc, s.volume_s, 325.068, "AbsTol", 0.01);   % rates stored to 3 decimals
verifyEqual(tc, s.cuts.kind(7:8), ["sails"; "sails"]);
end

function testVcp212SailsMatchesKOAXandKHGX(tc)
% Faster Doppler cuts move the SAILS insert after the 1.8 deg cut, as the
% RDA did at KOAX and KHGX.
el = [0.5 0.5 0.9 0.9 1.3 1.3 1.8 0.5 0.5 2.4 3.1 4.0 5.1 6.4 8.0 10.0 12.5 15.6 19.5]';
for c = [17.106 14.458; 294.275 309.691]
    s = nexrad_vcp_rev1(212, sails=1, cd_rate_deg_s=c(1));
    verifyEqual(tc, s.tilt_deg, el, "AbsTol", 1e-12);
    verifyEqual(tc, s.volume_s, c(2), "AbsTol", 0.01);
end
end

function testVcp112SailsCountsAndTime(tc)
s = nexrad_vcp_rev1(112, sails=1);
verifyEqual(tc, height(s.cuts), 23);                % ConOps Table 3
verifyEqual(tc, s.volume_s, 382.272, "AbsTol", 0.01);   % KTBW sum; placement differs (see help)
verifyEqual(tc, height(nexrad_vcp_rev1(112).cuts), 20);
verifyEqual(tc, height(nexrad_vcp_rev1(112, avset_top_deg=6.4).cuts), 15);
end

function testVcpBaseCutCounts(tc)
% ROC NEXRAD Technical Information elevations; split cuts add one cut each
% at 0.5, 0.9 and 1.3 deg (two each for VCP 112).
v = [35 215 212 12 112 31 32 34];
n_el = [9 15 14 14 14 5 5 7];
n_cut = [12 18 17 17 20 8 7 10];
for i = 1:numel(v)
    s = nexrad_vcp_rev1(v(i));
    verifyEqual(tc, numel(unique(s.tilt_deg)), n_el(i), sprintf("VCP %d elevations", v(i)));
    verifyEqual(tc, height(s.cuts), n_cut(i), sprintf("VCP %d cuts", v(i)));
    verifyEqual(tc, s.cuts.dur_s, 360./s.rate_deg_s, "AbsTol", 1e-12);
    verifyLessThanOrEqual(tc, max(s.rate_deg_s), 36);   % pedestal limit, ROC_TECH
end
verifyEqual(tc, max(nexrad_vcp_rev1(35).tilt_deg), 6.4);
verifyEqual(tc, max(nexrad_vcp_rev1(215).tilt_deg), 19.5);
end

function testVcpAvsetDropsHighCuts(tc)
s = nexrad_vcp_rev1(215, avset_top_deg=6.4);
verifyEqual(tc, max(s.tilt_deg), 6.4);
verifyEqual(tc, height(s.cuts), 12);
s = nexrad_vcp_rev1(212, avset_top_deg=10);
verifyEqual(tc, max(s.tilt_deg), 10);
end

function testVcpMrleInsertsLowestGroups(tc)
b = nexrad_vcp_rev1(212);
s = nexrad_vcp_rev1(212, mrle=4);
m = s.cuts(s.cuts.kind == "mrle", :);
verifyEqual(tc, m.el_deg, [0.5 0.5 0.9 0.9 1.3 1.3 1.8]');
verifyEqual(tc, height(s.cuts), height(b.cuts) + 7);
verifyEqual(tc, s.volume_s, b.volume_s + sum(m.dur_s), "AbsTol", 1e-9);
% Inserted near the middle of the base volume
i0 = find(s.cuts.kind == "mrle", 1);
t_before = sum(s.cuts.dur_s(1:i0-1));
verifyLessThan(tc, abs(t_before - b.volume_s/2), 20);
end

function testVcpSailsx3EvenlySpaced(tc)
s = nexrad_vcp_rev1(212, sails=3);
verifyEqual(tc, sum(s.cuts.kind == "sails"), 6);
t_start = cumsum([0; s.cuts.dur_s(1:end-1)]);
t05 = t_start(s.cuts.el_deg == 0.5 & s.cuts.waveform == "CS");   % 4 per volume
verifyEqual(tc, numel(t05), 4);
gaps = diff([t05; s.volume_s]);
verifyLessThan(tc, (max(gaps) - min(gaps))/mean(gaps), 0.6);
end

function testVcpCdRateOverride(tc)
s = nexrad_vcp_rev1(215, cd_rate_deg_s=17.106);
verifyEqual(tc, s.rate_deg_s(s.cuts.waveform == "CD"), 17.106*ones(3, 1));
verifyEqual(tc, s.rate_deg_s(s.cuts.waveform == "CS"), [11.470 13.370 15.919]');   % unchanged
end

function testVcpClearAirLongShortPulseMatchDecoded(tc)
% VCPs 31, 32 and 34 decoded from Level II Message Type 5 (see help).
el = {[0.5 0.5 1.5 1.5 2.5 2.5 3.5 4.5]', ...
    [0.5 0.5 1.5 1.5 2.5 3.5 4.5]', ...
    [0.5 0.5 0.9 0.9 1.3 1.3 1.8 2.4 3.1 4.5]'};
r = {[5.043 5.065 5.043 5.065 5.043 5.065 5.065 5.065]', ...
    [4.966 4.548 4.966 4.548 4.065 4.065 4.065]', ...
    [5.043 8.492 5.043 8.492 5.449 8.492 5.878 8.492 8.492 8.492]'};
t = [569.538 568.980 524.441];                     % sum of 360/rate (s)
pw = ["long" "short" "long"];
v = [31 32 34];
for i = 1:3
    s = nexrad_vcp_rev1(v(i));
    verifyEqual(tc, s.tilt_deg, el{i}, "AbsTol", 1e-12);
    verifyEqual(tc, s.rate_deg_s, r{i}, "AbsTol", 1e-12);
    verifyEqual(tc, s.volume_s, t(i), "AbsTol", 0.01);
    verifyEqual(tc, s.volume_s, sum(360./r{i}), "AbsTol", 1e-9);
    verifyEqual(tc, s.pulse_width, pw(i));
end
end

function testVcpClearAirOptionLimits(tc)
% ROC: SAILSx1 on VCP 34; no SAILS, MRLE or AVSET on 31 or 32.
s = nexrad_vcp_rev1(34, sails=1);
verifyEqual(tc, height(s.cuts), 12);
verifyEqual(tc, sum(s.cuts.kind == "sails"), 2);
verifyEqual(tc, s.volume_s, 524.441 + 360/5.043 + 360/8.492, "AbsTol", 0.01);
verifyError(tc, @() nexrad_vcp_rev1(34, sails=2), "radar_coupling:vcpSails");
for v = [31 32]
    verifyError(tc, @() nexrad_vcp_rev1(v, sails=1), "radar_coupling:vcpSails");
    verifyError(tc, @() nexrad_vcp_rev1(v, mrle=2), "radar_coupling:vcpMrle");
    verifyError(tc, @() nexrad_vcp_rev1(v, avset_top_deg=6.4), "radar_coupling:vcpAvset");
end
for v = [31 32 34 112]
    verifyError(tc, @() nexrad_vcp_rev1(v, cd_rate_deg_s=15), "radar_coupling:vcpCdRate");
end
end

function testVcpLimitsAndSourceReported(tc)
for v = [12 31 32 34 35 112 212 215]
    s = nexrad_vcp_rev1(v);
    verifyTrue(tc, all(isfield(s.limits, ["sails" "mrle" "avset" "cd_rate"])));
    verifyTrue(tc, startsWith(s.source_file, "K") && endsWith(s.source_file, "_V06"));
    verifyTrue(tc, any(s.pulse_width == ["short" "long"]));
end
end

function testVcpInputErrors(tc)
verifyError(tc, @() nexrad_vcp_rev1(11), "radar_coupling:badVcp");
verifyError(tc, @() nexrad_vcp_rev1(215, sails=2), "radar_coupling:vcpSails");
verifyError(tc, @() nexrad_vcp_rev1(35, mrle=2), "radar_coupling:vcpMrle");
verifyError(tc, @() nexrad_vcp_rev1(212, mrle=5), "radar_coupling:vcpMrle");
verifyError(tc, @() nexrad_vcp_rev1(212, mrle=2, sails=1), "radar_coupling:vcpSailsMrle");
verifyError(tc, @() nexrad_vcp_rev1(35, avset_top_deg=6.4), "radar_coupling:vcpAvset");
verifyError(tc, @() nexrad_vcp_rev1(212, avset_top_deg=5), "radar_coupling:vcpAvset");
verifyError(tc, @() nexrad_vcp_rev1(112, cd_rate_deg_s=15), "radar_coupling:vcpCdRate");
verifyError(tc, @() nexrad_vcp_rev1(212, cd_rate_deg_s=40), "radar_coupling:vcpCdRate");
end

function testVcpRunsInModelTimeWeighted(tc)
% Random mode draws each cut with probability dur/volume: the fraction of
% samples at 19.5 deg (one sample per tilt state, read back through a
% pattern that is -60 dB everywhere except within 0.01 deg of the axis
% pointing straight up) is checked indirectly via the gain at the path.
s = nexrad_vcp_rev1(35);
ant = struct("type", "offaxis_table", "g_max_dBi", 40, ...
    "offaxis_deg", [0 0.6 0.61 180], "gain_dB", [0 0 -60 -60]);
sys1 = struct("tx", ant, "scan", s);
sys2 = struct("tx", struct("type", "offaxis_table", "g_max_dBi", 10, ...
    "offaxis_deg", [0 180], "gain_dB", [0 0]), ...
    "scan", struct("type", "mech_rotating", "rate_deg_s", 75));
out = radar_mutual_gain_cdf_rev1(sys1, sys2, struct("seed", 7, "n_samples", 2e6));
% Main beam (flat to 0.6 deg off axis) reaches the horizon path only on
% the 0.5 deg cuts, while az within acosd(cos(0.6)/cos(0.5))... of the
% link: fraction = sum(dur at 0.5)/volume * (2*az_half/360).
az_half = acosd(cosd(0.6)/cosd(0.5));
p_expect = sum(s.cuts.dur_s(s.cuts.el_deg == 0.5))/s.volume_s*2*az_half/360;
p_hat = mean(out.g12_dBi > 40);
verifyEqual(tc, p_hat, p_expect, "RelTol", 0.05);
end

% =========================================================================
% offaxis_table antenna type
% =========================================================================
function testOffaxisTableHandValues(tc)
ant = struct("type", "offaxis_table", "g_max_dBi", 45, ...
    "offaxis_deg", [0 1 10 180], "gain_dB", [0 -12 -40 -60]);
p = antenna_prepare_rev1(ant);
z = zeros(3, 1);
st = struct("mech_az", [0; 0.5; 5], "tilt", z, "steer_az", z, "steer_el", z);
g = antenna_gain_rev1(p, 0, 0, st);
verifyEqual(tc, g, [45; 45 - 6; 45 - 12 - 28*4/9], "AbsTol", 1e-9);
% Tilt only: off-axis angle = tilt for a horizon target straight ahead
st = struct("mech_az", 0, "tilt", 1, "steer_az", 0, "steer_el", 0);
verifyEqual(tc, antenna_gain_rev1(p, 0, 0, st), 33, "AbsTol", 1e-9);
% 3 dB beamwidth from the table: -3 dB at 0.25 deg -> 0.5 deg
verifyEqual(tc, p.az_bw_deg, 0.5, "AbsTol", 1e-9);
end

function testOffaxisTableErrors(tc)
bad = {struct("offaxis_deg", [0 90], "gain_dB", [0 -10]), ...      % stops short of 180
    struct("offaxis_deg", [0 180], "gain_dB", [-1 -10]), ...       % no 0 dB peak
    struct("offaxis_deg", [0 90 90 180], "gain_dB", [0 -1 -2 -3]), ... % not strictly increasing
    struct("offaxis_deg", [0 180], "gain_dB", [0 2])};             % above peak
for i = 1:numel(bad)
    a = bad{i};
    a.type = "offaxis_table";
    a.g_max_dBi = 40;
    verifyError(tc, @() antenna_prepare_rev1(a), "radar_coupling:badOffaxisTable");
end
end

% =========================================================================
% nexrad_asr11_mag_rev1
% =========================================================================
function testAssumptionsAllCited(tc)
out = quick_run(struct());
C = out.citations;
src = unique(out.assumptions.source(out.assumptions.source ~= ""));
verifyTrue(tc, all(ismember(src, C.key)), "every assumption source has a citation");
web = ~ismember(C.key, ["WP_REV3" "LOE_20261104"]);
verifyTrue(tc, all(startsWith(C.url(web), "https://")), "every public source has an https URL");
verifyEqual(tc, numel(unique(C.key)), height(C));
end

function testNexradPatternEnvelopes(tc)
% NSSL_1998 envelope: -35 dB at 2 deg, -48 at 6, -55 at 20 and beyond.
p = antenna_prepare_rev1(quick_run(struct()).sys_nexrad.tx);
off = [0 0.4625 2 4 6 20 90 180];
verifyEqual(tc, p.off_interp(off), [0 -3 -35 -41.5 -48 -55 -55 -55], "AbsTol", 1e-3);  % 0.005 deg table grid
% ROC_TECH spec: -29 dB beyond the main beam, -40 dB beyond 10 deg.
p = antenna_prepare_rev1(quick_run(struct("nexrad_pattern", "roc_spec")).sys_nexrad.tx);
verifyEqual(tc, p.off_interp([5 9.9 10.5 180]), [-29 -29 -40 -40], "AbsTol", 1e-9);
verifyEqual(tc, p.g_max_dBi, 45.5);
end

function testAsr11PatternValues(tc)
p = antenna_prepare_rev1(quick_run(struct()).sys_asr11.tx);
z = zeros(4, 1);
st = struct("mech_az", z, "tilt", z, "steer_az", z, "steer_el", z);
g = antenna_gain_rev1(p, 0, [0; 10; 20; 0], st);
verifyEqual(tc, g(1), 34, "AbsTol", 1e-9);                       % M.1464-2 Radar C peak
verifyEqual(tc, g(2) - g(3), 20*log10(sind(20)/sind(10)), "AbsTol", 1e-9);  % CSC^2 slope
st.mech_az = [0.725; 180; 90; 0];
g = antenna_gain_rev1(p, 0, 0, st);
verifyEqual(tc, g(1), 31, "AbsTol", 0.1);   % 1.45 deg 3 dB beamwidth (M.1851-2 Table 4 K gives -3.07 dB)
verifyEqual(tc, g(2:3), [3.5; 3.5], "AbsTol", 1e-9);            % remote sidelobes +3.5 dBi
end

function testAsrRateVariancePerTrial(tc)
% One rate per trial, inside 75 +/- 0.5 deg/s (M.2069-0 Annex 1).
out = quick_run(struct("duration_s", 120, "n_trials", 3));
R = out.scan_max.asr_rev;
for j = 1:3
    r = unique(R.rate_deg_s(R.trial == j));
    verifyNumElements(tc, r, 1);
    verifyLessThanOrEqual(tc, abs(r - 75), 0.5);
end
verifyNumElements(tc, unique(R.rate_deg_s), 3);
end

function testAsrRateVariancePerRev(tc)
% New rate every revolution; the revolution period follows the rate.
out = quick_run(struct("duration_s", 300, "asr_jitter_mode", "per_rev", "asr_rate_jitter_deg_s", 0.5));
R = out.scan_max.asr_rev(out.scan_max.asr_rev.complete, :);
verifyLessThanOrEqual(tc, max(abs(R.rate_deg_s - 75)), 0.5);
verifyGreaterThan(tc, numel(unique(R.rate_deg_s)), 0.9*height(R));
period = diff(R.t_start_s);
verifyEqual(tc, period, 360./R.rate_deg_s(1:end-1), "AbsTol", out.opts.dt_s);
verifyGreaterThanOrEqual(tc, min(period), 360/75.5 - out.opts.dt_s);
verifyLessThanOrEqual(tc, max(period), 360/74.5 + out.opts.dt_s);
end

function testScanMaximaConsistent(tc)
out = quick_run(struct("duration_s", 300, "asr_rate_jitter_deg_s", 0));
S = out.scan_max;
% The three groupings partition the same samples: equal overall maxima.
verifyEqual(tc, max(S.asr_rev.max_g21_dBi), max(S.nexrad_cut.max_g21_dBi), "AbsTol", 1e-12);
verifyEqual(tc, max(S.asr_rev.max_g21_dBi), max(S.nexrad_volume.max_g21_dBi), "AbsTol", 1e-12);
% Each volume maximum is the maximum of its cut maxima.
for v = unique(S.nexrad_cut.volume)'
    c = S.nexrad_cut(S.nexrad_cut.volume == v, :);
    verifyEqual(tc, S.nexrad_volume.max_g21_dBi(S.nexrad_volume.volume == v), max(c.max_g21_dBi), "AbsTol", 1e-12);
end
% Tx = rx patterns: g12 = g21.
verifyEqual(tc, S.asr_rev.max_g12_dBi, S.asr_rev.max_g21_dBi, "AbsTol", 1e-9);
% Revolutions last 4.8 s (12.5 rpm).
verifyEqual(tc, S.t_rev_s, 4.8, "AbsTol", 1e-12);
d = diff(S.asr_rev.t_start_s(S.asr_rev.complete));
verifyEqual(tc, d, 4.8 + zeros(size(d)), "AbsTol", out.opts.dt_s);
verifyLessThanOrEqual(tc, max(S.asr_rev.max_g21_dBi), 45.5 + 34);
end

function testScanMaximaPhysicalBounds(tc)
% Every ASR-11 revolution sweeps its main beam across the NEXRAD, so its
% maximum is at least 34 dBi + the NEXRAD floor (45.5 - 55 dBi). Every
% NEXRAD cut sweeps its beam across the ASR-11, so its maximum is at least
% the NEXRAD gain at that cut's tilt + the ASR-11 floor (+3.5 dBi). The
% 0.002 s step costs < 0.2 dB at either peak.
out = quick_run(struct("duration_s", 600));
R = out.scan_max.asr_rev(out.scan_max.asr_rev.complete, :);
verifyGreaterThanOrEqual(tc, min(R.max_g21_dBi), 34 + 45.5 - 55 - 0.2);
C = out.scan_max.nexrad_cut(out.scan_max.nexrad_cut.complete, :);
p = antenna_prepare_rev1(out.sys_nexrad.tx);
lower = 45.5 + p.off_interp(C.el_deg) + 3.5 - 0.2;
verifyGreaterThanOrEqual(tc, C.max_g21_dBi, lower);
% Cuts at 0.5 deg reach main beam to main beam when the ASR-11 beam is on
% the NEXRAD as the NEXRAD beam passes (not every cut); the overall
% maximum stays below the two peak gains.
verifyLessThan(tc, max(C.max_g21_dBi), 45.5 + 34);
end

function testTimeSeriesMatchesRandomCdf(tc)
out = quick_run(struct("duration_s", 3600, "n_trials", 2, "n_samples", 2e6));
q_r = out.cdf.pct21.gain_dBi(out.cdf.pct21.pct == 99);
q_t = out.ts_pct.gain_dBi(out.ts_pct.pct == 99);
verifyEqual(tc, q_t, q_r, "AbsTol", 0.5);
q_r = out.cdf.pct21.gain_dBi(out.cdf.pct21.pct == 50);
q_t = out.ts_pct.gain_dBi(out.ts_pct.pct == 50);
verifyEqual(tc, q_t, q_r, "AbsTol", 0.2);
end

function testPaperBaselineUsesF1245ConstantRate(tc)
out = quick_run(struct("nexrad_scan", "constant", "nexrad_pattern", "f1245", "asr_pattern", "f1245"));
verifyEqual(tc, out.sys_nexrad.tx.g_max_dBi, 46);
verifyEqual(tc, out.sys_asr11.tx.g_max_dBi, 34);
C = out.scan_max.nexrad_cut;
verifyEqual(tc, unique(C.el_deg), 0);
d = diff(C.t_start_s(C.complete));
verifyEqual(tc, d, 20 + zeros(size(d)), "AbsTol", out.opts.dt_s);   % 3 rpm
end

function testRefLevelsHandValues(tc)
% Defaults: NSSL NEXRAD envelope, M.1851-2 CSC^2 ASR-11, horizon path,
% NEXRAD lowest cut 0.5 deg, ASR-11 tilt 0.
out = quick_run(struct());
R = out.ref_levels;
g_n_mb = 45.5 - 12*(0.5/0.925)^2;            % NEXRAD main beam on path at 0.5 deg tilt
verifyEqual(tc, R.nexrad_gain_dBi(1), g_n_mb, "AbsTol", 1e-3);
verifyEqual(tc, R.asr11_gain_dBi(1), 34, "AbsTol", 1e-9);
verifyEqual(tc, R.mag_dBi(1), g_n_mb + 34, "AbsTol", 1e-3);          % main beam to main beam
verifyEqual(tc, R.mag_dBi(3), g_n_mb + 3.5, "AbsTol", 1e-3);         % + ASR-11 remote sidelobe
verifyEqual(tc, R.mag_dBi(4), 34 + 45.5 - 35, "AbsTol", 1e-6);       % + NEXRAD first sidelobe (NSSL)
verifyEqual(tc, R.mag_dBi(5), 34 + 45.5 - 55, "AbsTol", 1e-6);       % + NEXRAD back (NSSL)
% ASR-11 near sidelobe sits between its floor and its -20 dB taper bound
verifyGreaterThanOrEqual(tc, R.asr11_gain_dBi(2), 3.5);
verifyLessThanOrEqual(tc, R.asr11_gain_dBi(2), 34 - 20);
% Ordering and consistency with the CDF
verifyEqual(tc, R.mag_dBi(1), max(R.mag_dBi));
verifyGreaterThanOrEqual(tc, R.time_exceeded_pct, 0);
verifyLessThanOrEqual(tc, R.time_exceeded_pct(1), min(R.time_exceeded_pct(2:end)));
% Every complete ASR-11 revolution reaches its main beam + NEXRAD back
% lobe, less the < 0.2 dB sampling loss of the 0.002 s step
Rv = out.scan_max.asr_rev(out.scan_max.asr_rev.complete, :);
verifyGreaterThanOrEqual(tc, Rv.max_g21_dBi, R.mag_dBi(5) - 0.2);
verifyEqual(tc, R.asr_rev_reached_pct(5), 100);
end

function testSensitivityTable(tc)
% Peak-gain MB-MB is 45.5 + 34 = 79.5 dBi; every other choice is lower,
% and the delta column is the difference from it.
out = quick_run(struct());
S = out.sensitivity;
verifyEqual(tc, S.mag_dBi(1), 79.5, "AbsTol", 1e-9);
verifyEqual(tc, S.delta_vs_peak_mbmb_dB, S.mag_dBi - 79.5, "AbsTol", 1e-12);
verifyEqual(tc, S.mag_dBi(2), out.ref_levels.mag_dBi(1), "AbsTol", 1e-12);
verifyLessThanOrEqual(tc, S.mag_dBi(2:end), S.mag_dBi(1));
q = S.mag_dBi(4:end);
verifyTrue(tc, issorted(q), "percentile rows increase");
verifyEqual(tc, S.time_exceeded_pct(1), 0);
end

function testModeEnvelopeBoundsAnyMixture(tc)
% Two synthetic modes; the envelope percentile bounds the percentile of
% every mixture of them.
rs = RandStream("mrg32k3a", "Seed", 3);
a = randn(rs, 2e5, 1);
b = 1 + 3*randn(rs, 2e5, 1);
env = mag_mode_envelope_rev1({a, b}, ["A", "B"], [90 99 99.9]);
verifyEqual(tc, env.pct.envelope_dBi, max(prctile(a, [90 99 99.9]'), prctile(b, [90 99 99.9]')), ...
    "AbsTol", 1e-12);
verifyEqual(tc, env.pct.worst_mode, ["B"; "B"; "B"]);
for w = [0 0.1 0.5 0.9 1]
    k = round(w*2e5);
    mix = [a(1:k); b(1:2e5-k)];
    verifyLessThanOrEqual(tc, prctile(mix, [90 99 99.9]'), env.pct.envelope_dBi + 1e-9);
end
verifyGreaterThanOrEqual(tc, env.exceed.envelope_exceed_pct, 0);
verifyTrue(tc, issorted(flipud(env.exceed.envelope_exceed_pct)));
verifyError(tc, @() mag_mode_envelope_rev1({a}, ["A", "B"]), "radar_coupling:envelopeInput");
end

function testCdfArray(tc)
% CDF as [percent, MAG dBi], 0 to 100, consistent with the samples.
out = quick_run(struct("cdf_pct", 0:5:100));
A = out.cdf.cdf_array21;
verifySize(tc, A, [21 2]);
verifyEqual(tc, A(:, 1), (0:5:100)');
verifyEqual(tc, A(1, 2), min(out.cdf.g21_dBi));
verifyEqual(tc, A(end, 2), max(out.cdf.g21_dBi));
verifyTrue(tc, issorted(A(:, 2)));
verifyEqual(tc, A(A(:, 1) == 50, 2), prctile(out.cdf.g21_dBi, 50));
d = quick_run(struct()).cdf.cdf_array21;          % default grid includes the tail
verifyTrue(tc, all(ismember([0 50 99 99.9 99.99 99.999 100], d(:, 1))));
verifyError(tc, @() quick_run(struct("cdf_pct", [0 50 40])), "radar_coupling:badCdfPct");
verifyError(tc, @() quick_run(struct("cdf_pct", [0 101])), "radar_coupling:badCdfPct");
end

function testPeakDecomposition(tc)
% Each revolution maximum is the NEXRAD gain plus the ASR-11 gain at the
% sample where it occurs, and the offsets are on -180..180 deg.
out = quick_run(struct("duration_s", 300, "keep_ts1", true));
R = out.scan_max.asr_rev;
verifyEqual(tc, R.max_g21_dBi, R.nexrad_gain_at_max_dBi + R.asr11_gain_at_max_dBi, "AbsTol", 1e-9);
verifyLessThanOrEqual(tc, abs([R.nexrad_az_off_at_max_deg; R.asr11_az_off_at_max_deg]), 180);
verifyLessThanOrEqual(tc, R.nexrad_gain_at_max_dBi, 45.5 + 1e-9);
verifyLessThanOrEqual(tc, R.asr11_gain_at_max_dBi, 34 + 1e-9);
T1 = out.scan_max.ts1;
verifyEqual(tc, T1.g21_dBi, T1.asr11_gain_dBi + T1.nexrad_gain_dBi, "AbsTol", 1e-9);
k = find(T1.t_s == R.t_at_max_s(1), 1);
verifyEqual(tc, T1.g21_dBi(k), R.max_g21_dBi(1), "AbsTol", 1e-12);
verifyEmpty(tc, quick_run(struct()).scan_max.ts1);
end

function testPlateauCheckMatchesSimulation(tc)
% VCP 212 SAILSx3: 8 cuts at 0.5 deg per 416.6 s volume, so about
% 8 x 4.8 / 416.6 = 9.2% of ASR-11 revolutions reach 45.49 dBi.
out = quick_run(struct("vcp", 212, "sails", 3, "duration_s", 3600, "asr_rate_jitter_deg_s", 0));
P = out.scan_max.plateaus;
i05 = P.el_deg == 0.5;
verifyEqual(tc, P.cuts_per_volume(i05), 8);
verifyEqual(tc, P.level_dBi(i05), 45.5 - 12*(0.5/0.925)^2 + 3.5, "AbsTol", 1e-3);
verifyEqual(tc, P.predicted_rev_pct(i05), 100*8*4.8/out.scan_max.t_volume_s, "AbsTol", 1e-9);
verifyEqual(tc, P.simulated_rev_pct(i05), P.predicted_rev_pct(i05), "AbsTol", 3.5);
i09 = P.el_deg == 0.9;
verifyEqual(tc, P.level_dBi(i09), 45.5 - 12*(0.9/0.925)^2 + 3.5, "AbsTol", 1e-3);
verifyEqual(tc, P.simulated_rev_pct(i09), P.predicted_rev_pct(i09), "AbsTol", 4);
end

function testPercentileConfidenceInterval(tc)
% Order-statistic 95% interval brackets the estimate; for a normal sample
% its width at the median is close to 2 x 1.96 x sqrt(pi/2)/sqrt(n).
out = quick_run(struct("n_samples", 2e5));
P = out.cdf.pct21;
verifyTrue(tc, all(P.ci95_lo_dBi <= P.gain_dBi & P.gain_dBi <= P.ci95_hi_dBi));
rs = RandStream("mrg32k3a", "Seed", 9);
x = randn(rs, 1e6, 1);
gs = sort(x);
n = numel(x); half = 1.96*sqrt(n*0.25);
w = gs(ceil(n/2 + half)) - gs(floor(n/2 - half));
verifyEqual(tc, w, 2*1.96*sqrt(pi/2)/sqrt(n), "RelTol", 0.15);
end

function testShortRunHasNoCompleteCuts(tc)
% 10 s is shorter than any VCP 212 cut pair: no complete cut, empty
% per-elevation table, no error.
out = quick_run(struct("duration_s", 10));
verifyEqual(tc, height(out.scan_max.by_elevation), 0);
verifyEqual(tc, out.scan_max.by_elevation.Properties.VariableNames, ...
    {'el_deg', 'n_cuts', 'median_max_g21_dBi', 'max_max_g21_dBi'});
end

function testPlotFlagDraws(tc)
f0 = findall(groot, "Type", "figure");
quick_run(struct("plot_flag", true));
f = setdiff(findall(groot, "Type", "figure"), f0);
tc.addTeardown(@() close(f));
verifyNumElements(tc, f, 1);
verifyNumElements(tc, findall(f, "Type", "ConstantLine"), 5);
end

function testReproducibleWithSeed(tc)
a = quick_run(struct("seed", 11));
b = quick_run(struct("seed", 11));
c = quick_run(struct("seed", 12));
verifyEqual(tc, a.scan_max.asr_rev, b.scan_max.asr_rev);
verifyEqual(tc, a.cdf.g21_dBi, b.cdf.g21_dBi);
verifyNotEqual(tc, a.cdf.g21_dBi, c.cdf.g21_dBi);
end

function testOptsErrors(tc)
verifyError(tc, @() nexrad_asr11_mag_rev1(struct()), "radar_coupling:missingSeed");
verifyError(tc, @() quick_run(struct("nexrad_pattern", "gaussian")), "radar_coupling:badNexradPattern");
verifyError(tc, @() quick_run(struct("asr_pattern", "csc")), "radar_coupling:badAsrPattern");
verifyError(tc, @() quick_run(struct("nexrad_scan", "sector")), "radar_coupling:badNexradScan");
verifyError(tc, @() quick_run(struct("dt_s", 0.05)), "radar_coupling:badDt");
verifyError(tc, @() quick_run(struct("asr_tilt_deg", 3)), "radar_coupling:asrTilt");
verifyError(tc, @() quick_run(struct("asr_jitter_mode", "per_scan")), "radar_coupling:badJitterMode");
verifyError(tc, @() quick_run(struct("asr_rate_jitter_deg_s", -1)), "radar_coupling:badJitter");
end

% =========================================================================
% Antenna concerns, 2.7 GHz Retuning LOE meeting of November 4, 2026
% =========================================================================
function testCombinedAvsetSailsAllowed(tc)
% ROC VCP Info Sheet: AVSET is compatible with SAILS and MRLE.
s = nexrad_vcp_rev1(212, sails=3, avset_top_deg=6.4);
b = nexrad_vcp_rev1(212, sails=3);
hi = b.cuts.el_deg > 6.4;
verifyEqual(tc, s.volume_s, b.volume_s - sum(b.cuts.dur_s(hi)), "AbsTol", 1e-9);
share = 100*sum(s.cuts.dur_s(s.cuts.el_deg == 0.5))/s.volume_s;
verifyEqual(tc, share, 100*183.07/s.volume_s, "AbsTol", 0.01);   % Step 5: 183.07 s at 0.5 deg
verifyGreaterThan(tc, share, 100*sum(b.cuts.dur_s(b.cuts.el_deg == 0.5))/b.volume_s);
m = nexrad_vcp_rev1(215, mrle=3, avset_top_deg=6.4);
verifyLessThan(tc, max(m.tilt_deg), 6.4 + 1e-9);
end

function testLobePairsSumAndHandValue(tc)
out = quick_run(struct("vcp", 212, "sails", 3, "n_samples", 2e4, "duration_s", 300));
L = out.lobe_pairs;
E = L.Properties.UserData;
verifyEqual(tc, height(L), 16);
verifyEqual(tc, sum(L.time_pct), 100, "AbsTol", 1e-6);
verifyEqual(tc, sum(L.time_pct_ts), 100, "AbsTol", 1e-6);
% Main-beam edges: NEXRAD main beam meets the -35 dB first sidelobe at
% 0.925*sqrt(35/12) deg; the ASR-11 main beam ends at the -24.5 dB level.
verifyEqual(tc, E.nexrad(1), 0.925*sqrt(35/12), "AbsTol", 2e-3);
verifyGreaterThan(tc, E.asr11(1), 1.45/2);
verifyLessThan(tc, E.asr11(1), E.asr11(2));
% Hand count of main beam to main beam time: NEXRAD cut at tilt e is in
% its main beam for |daz| <= acos(cos(edge)/cos(e)); ASR-11 for |daz| <= edge.
s = out.sys_nexrad.scan;
w = s.cuts.dur_s/s.volume_s;
c = cosd(E.nexrad(1))./cosd(s.cuts.el_deg);
fN = sum(w(c <= 1).*2.*acosd(c(c <= 1))/360);
fA = 2*E.asr11(1)/360;
mm = L.nexrad_region == "main beam" & L.asr11_region == "main beam";
verifyEqual(tc, L.time_pct(mm), 100*fN*fA, "RelTol", 2e-3);
% Largest MB-MB MAG is main beam to main beam on the path (Step 4)
verifyEqual(tc, L.max_mag_dBi(mm), out.ref_levels.mag_dBi(1), "AbsTol", 0.01);
verifyLessThanOrEqual(tc, L.max_mag_ts_dBi(isfinite(L.max_mag_ts_dBi)), ...
    L.max_mag_dBi(isfinite(L.max_mag_ts_dBi)) + 0.01);   % 0.002 deg grid of the exact maxima
end

function testPersistenceConsistent(tc)
out = quick_run(struct("vcp", 212, "sails", 3, "duration_s", 600));
P = out.scan_max.persistence;
R = out.scan_max.asr_rev(out.scan_max.asr_rev.complete, :);
for i = 1:height(P)
    n_hit = sum(R.max_g21_dBi >= P.threshold_dBi(i));
    verifyEqual(tc, P.rev_pct(i), 100*n_hit/height(R), "AbsTol", 1e-9);
    verifyEqual(tc, P.n_runs(i)*P.mean_run_revs(i), n_hit, "AbsTol", 1e-9);
    verifyLessThanOrEqual(tc, P.max_run_revs(i), n_hit);
end
verifyTrue(tc, issorted(flip(P.rev_pct)));
% NEXRAD/ASR-11: each cut crosses the path once, so a 70 dBi revolution
% (NEXRAD main beam needed) is never followed by another.
verifyLessThanOrEqual(tc, P.max_run_revs(P.threshold_dBi == 70), 1);
end

function testSectorScan(tc)
a = quick_run(struct("vcp", 35));
b = quick_run(struct("vcp", 35, "nexrad_sector_deg", 360));
verifyEqual(tc, b.cdf.g21_dBi, a.cdf.g21_dBi);      % full sector is the default
n = quick_run(struct("vcp", 35, "nexrad_sector_deg", 20, "n_samples", 1e5));
verifyGreaterThan(tc, n.cdf.pct21.gain_dBi(n.cdf.pct21.pct == 99), a.cdf.pct21.gain_dBi(a.cdf.pct21.pct == 99));
% NEXRAD main-beam time share scales with 360/sector
fm = @(o) sum(o.lobe_pairs.time_pct(o.lobe_pairs.nexrad_region == "main beam"));
verifyEqual(tc, fm(n)/fm(a), 360/20, "RelTol", 1e-3);
verifyEqual(tc, n.sys_nexrad.scan.sector_deg, 20);
verifyError(tc, @() quick_run(struct("nexrad_sector_deg", 0)), "radar_coupling:badSector");
verifyError(tc, @() quick_run(struct("nexrad_sector_deg", 400)), "radar_coupling:badSector");
end

function testBackLobeAndRemoteSidelobeKnobs(tc)
b = quick_run(struct("nexrad_back_dB", -45));
r = b.ref_levels;
verifyEqual(tc, r.mag_dBi(r.combination == "ASR-11 main beam + NEXRAD back lobe"), 34 + 45.5 - 45, "AbsTol", 1e-6);
a = quick_run(struct("asr_remote_sl_dBi", 6.5));
r = a.ref_levels;
verifyEqual(tc, r.asr11_gain_dBi(r.combination == "NEXRAD main beam + ASR-11 back lobe"), 6.5, "AbsTol", 1e-6);
d = quick_run(struct());
verifyEqual(tc, d.ref_levels.mag_dBi(5), 24.5, "AbsTol", 1e-6);        % defaults unchanged
verifyError(tc, @() quick_run(struct("nexrad_back_dB", -45, "nexrad_pattern", "roc_spec")), "radar_coupling:backNeedsNssl");
verifyError(tc, @() quick_run(struct("nexrad_back_dB", -30)), "radar_coupling:badBack");
verifyError(tc, @() quick_run(struct("asr_remote_sl_dBi", 6.5, "asr_pattern", "f1245")), "radar_coupling:remoteNeedsM1851");
verifyError(tc, @() quick_run(struct("asr_remote_sl_dBi", 12)), "radar_coupling:badRemote");
end

function testPriorPracticeAssumptionsRecorded(tc)
A = quick_run(struct()).assumptions;
k = startsWith(A.parameter, "Prior practice") | startsWith(A.parameter, "Current FAA practice");
verifyGreaterThanOrEqual(tc, sum(k), 5);
verifyTrue(tc, all(A.source(k) == "LOE_20261104"));
verifyEqual(tc, A.value(A.parameter == "NEXRAD far-field distance"), 2*8.53^2*2850e6/physconst("LightSpeed"), "RelTol", 1e-12);
end

function testAlignmentDwellHandValues(tc)
% Equal rates, offsets U(+/-0.5) each: delta is triangular on +/-1 deg/s,
% so P(|delta| < a) = 1 - (1 - a)^2.
D = alignment_dwell_rev1(75, 75, 1.45, 1.45);
verifyEqual(tc, D.ratio, [1 1]);
verifyEqual(tc, D.window_deg, 2.9);
verifyEqual(tc, D.episode0_s, Inf);
a = 2.9./D.p_exceed.t_s;
verifyEqual(tc, D.p_exceed.p_episode_longer, 1 - (1 - a).^2, "AbsTol", 1e-12);
% 2:1: delta = d1 - 2*d2, density 1/2 at 0, so P ~ a for small a.
D2 = alignment_dwell_rev1(75, 37.5, 1.45, 1.45);
verifyEqual(tc, D2.ratio, [2 1]);
verifyEqual(tc, D2.joint_period_s, 9.6, "AbsTol", 1e-12);
verifyEqual(tc, D2.window_deg, 4.35, "AbsTol", 1e-12);
verifyEqual(tc, D2.p_exceed.p_episode_longer(3), 4.35/3600, "RelTol", 1e-9);
% No jitter, nominal drift 0.01 deg/s: episode = W / delta
D3 = alignment_dwell_rev1(75, 75.01, 1.45, 1.45, jitter_deg_s=[0 0]);
verifyEqual(tc, D3.episode0_s, 290, "RelTol", 1e-9);
verifyEqual(tc, D3.p_exceed.p_episode_longer', double([60 600 3600 86400 2592000] < 290));
% A cap at the cut length removes longer episodes
D4 = alignment_dwell_rev1(75, 12.524, 1.45, 0.925, jitter_deg_s=[0.5 0], max_episode_s=28.745);
verifyTrue(tc, all(D4.p_exceed.p_episode_longer(D4.p_exceed.t_s >= 60) == 0));
verifyError(tc, @() alignment_dwell_rev1(75, 75, 1.45, 1.45, jitter_deg_s=[80 0]), "radar_coupling:badJitter");
end

function testAlignmentDwellBruteForce(tc)
% Count joint periods with both beams inside their 3 dB beamwidths at the
% same instant; compare with T_ep / T_j.
bw = 1.45;
cases = {75, 75.01, -2, 400; 75, 37.495, 1.5, 600};
for k = 1:2
    [r1, r2, ph0, T] = cases{k, :};
    D = alignment_dwell_rev1(r1, r2, bw, bw, jitter_deg_s=[0 0]);
    t = (0:0.0005:T)';
    a1 = wrapTo180(r1*t);
    a2 = wrapTo180(ph0 + r2*t);                 % radar 2 points at radar 1 when a2 = 0 here
    hit = abs(a1) <= bw/2 & abs(a2) <= bw/2;
    n_ev = numel(unique(floor(t(hit)/D.joint_period_s)));
    verifyEqual(tc, n_ev, D.events0, "AbsTol", 2, sprintf("case %d", k));
end
end

function testFirstTrialRepeatsTrial(tc)
% Trial j uses substream 1000 + j, so first_trial = 3, n_trials = 1
% repeats trial 3 of a three-trial run exactly.
a = quick_run(struct("n_trials", 3, "duration_s", 30));
b = quick_run(struct("first_trial", 3, "n_trials", 1, "duration_s", 30, "keep_ts1", true));
A = a.scan_max.asr_rev(a.scan_max.asr_rev.trial == 3, :);
verifyEqual(tc, b.scan_max.asr_rev, A);
verifyEqual(tc, max(b.scan_max.ts1.g21_dBi), max(A.max_g21_dBi), "AbsTol", 1e-12);
verifyError(tc, @() quick_run(struct("first_trial", 0)), "radar_coupling:badFirstTrial");
end

function testSensitivityFunction(tc)
b = struct("seed", 3, "vcp", 212, "sails", 3, "n_samples", 2e4, "duration_s", 120, "n_trials", 1);
c = struct("group", {"ASR-11 tilt", "ASR-11 remote sidelobe"}, "label", {"+2.5 deg", "+6.5 dBi"}, ...
    "opts", {struct("asr_tilt_deg", 2.5), struct("asr_remote_sl_dBi", 6.5)});
T = mag_sensitivity_rev1(b, c);
verifyEqual(tc, height(T), 3);
verifyEqual(tc, T.label(1), "Baseline");
verifyEqual(tc, T{1, ["d99_dB", "d99_9_dB", "d99_99_dB", "d99_999_dB"]}, [0 0 0 0]);
verifyEqual(tc, T.mbmb_peak_dBi, 79.5 + [0; 0; 0], "AbsTol", 1e-9);
% Same seed, same pointing samples: a higher ASR-11 floor can only raise each MAG.
verifyGreaterThan(tc, T.p99_9_dBi(3), T.p99_9_dBi(1));
verifyError(tc, @() mag_sensitivity_rev1(b, struct("label", "x")), "radar_coupling:badCase");
end

% =========================================================================
function out = quick_run(o)
% Small run for tests; fields of o override.
d = struct("seed", 5, "n_samples", 2e4, "duration_s", 60, "n_trials", 1);
f = fieldnames(o);
for i = 1:numel(f)
    d.(f{i}) = o.(f{i});
end
out = nexrad_asr11_mag_rev1(d);
end
