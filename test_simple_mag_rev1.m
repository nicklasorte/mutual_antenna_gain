function tests = test_simple_mag_rev1
%TEST_SIMPLE_MAG_REV1 Tests for simple_mag_nexrad_asr11_rev1.
%   Hand-worked values for the two patterns and the MAG distribution, and a
%   comparison with the full model (nexrad_asr11_mag_rev1) at the same
%   settings.
%
%   DRAFT//PRE-DECISIONAL
tests = functiontests(localfunctions);
end

function testPatternPeaksAndLevels(tc)
o = simple_mag_nexrad_asr11_rev1();
p = o.pattern;
g = @(col, a) p.(col)(abs(p.az_off_deg - a) < 1e-9);
verifyEqual(tc, g("g_nexrad_dBi", 0), 45.5 - 12*(0.5/0.925)^2, "AbsTol", 1e-9, "main beam at the 0.5 deg tilt");
verifyEqual(tc, g("g_asr11_dBi", 0), 34, "AbsTol", 1e-12);
verifyEqual(tc, g("g_nexrad_dBi", 5), 16.5, "AbsTol", 1e-12, "29 dB below peak out to 10 deg");
verifyEqual(tc, g("g_nexrad_dBi", 90), 5.5, "AbsTol", 1e-12, "40 dB below peak beyond 10 deg");
verifyEqual(tc, g("g_asr11_dBi", 3), 9.5, "AbsTol", 1e-12, "first sidelobe");
verifyEqual(tc, g("g_asr11_dBi", 4), 3.5, "AbsTol", 1e-12, "remote sidelobe");
verifyEqual(tc, g("g_asr11_dBi", 1), 34 - 12*(1/1.45)^2, "AbsTol", 1e-12);
verifyEqual(tc, o.peak_dBi, 45.5 - 12*(0.5/0.925)^2 + 34, "AbsTol", 1e-9);
end

function testDistributionHandValues(tc)
o = simple_mag_nexrad_asr11_rev1();
verifyEqual(tc, sum(o.dist.prob), 1, "AbsTol", 1e-12);
verifyTrue(tc, all(diff(o.dist.exceed) <= 1e-15), "exceedance falls");
verifyEqual(tc, o.pct.mag_dBi(o.pct.pct == 50), 9.0, "AbsTol", 1e-3, "far sidelobes: 5.5 + 3.5 dBi");
verifyLessThanOrEqual(tc, max(o.dist.mag_dBi), o.peak_dBi + 1e-3, "nothing above both main beams");
% Share of time at or above 70 dBi, by hand: only both main beams reach
% it. For a WSR-88D pointing a, the ASR-11 must be within
% b(a) = 1.45 sqrt((gn(a) + 34 - 70)/12) of the path, so the share is the
% integral of 2 b(a) over a, divided by 360^2.
gn = @(a) 45.5 - 12*(acosd(cosd(0.5)*cosd(a))/0.925).^2;
b = @(a) 1.45*sqrt(max(gn(a) + 34 - 70, 0)/12);
hand = 100*integral(@(a) 2*b(a), -1, 1, "AbsTol", 1e-12)/360^2;
verifyEqual(tc, o.exceed.time_pct(o.exceed.thr_dBi == 70), hand, "RelTol", 0.01);
end

function testMatchesFullModelAtZeroTilt(tc)
% The full model's constant-rate WSR-88D (tilt 0) with the ROC
% specification pattern and the M.1851-2 ASR-11 pattern.
o = simple_mag_nexrad_asr11_rev1(nexrad_tilt_deg=0);
verifyEqual(tc, o.peak_dBi, 79.5, "AbsTol", 1e-9);
w = nexrad_asr11_mag_rev1(struct("seed", 1, "nexrad_scan", "constant", "n_samples", 1e3, "duration_s", 60, "n_trials", 1));
F = w.cdf.pct21;
% The full model's CDF is high by at most 0.3 dB from its grid; the simple
% ASR-11 pattern differs from M.1851-2 by up to 3 dB from 1.4 to 3.55 deg,
% which moves the 99.99% value (ASR-11 main-beam skirt on the WSR-88D
% sidelobe) by about 0.5 dB.
for p = [50 90 99 99.9 99.999]
    verifyEqual(tc, o.pct.mag_dBi(o.pct.pct == p), F.gain_dBi(F.pct == p), "AbsTol", 0.3, sprintf("%g%%", p));
end
verifyEqual(tc, o.pct.mag_dBi(o.pct.pct == 99.99), F.gain_dBi(F.pct == 99.99), "AbsTol", 0.6, "99.99%");
end

function testRatesAndTimeSeries(tc)
v = nexrad_vcp_rev1(212, sails=3, avset_top_deg=6.4);
o = simple_mag_nexrad_asr11_rev1(nexrad_rate="vcp212");
verifyEqual(tc, o.params.value(o.params.parameter == "WSR-88D rotation rate"), v.cuts.rate_deg_s(1), "AbsTol", 1e-3);
o = simple_mag_nexrad_asr11_rev1(align_s=30);
[mx, im] = max(o.ts.mag_dBi);
verifyEqual(tc, mx, o.peak_dBi, "AbsTol", 1e-9);
verifyEqual(tc, o.ts.t_s(im), 30, "AbsTol", 1e-9);
verifyEqual(tc, o.repeat_s, 120, "AbsTol", 1e-9, "18 and 75 deg/s repeat every 120 s");
r1 = simple_mag_nexrad_asr11_rev1(align_s=NaN, seed=3);
r2 = simple_mag_nexrad_asr11_rev1(align_s=NaN, seed=3);
verifyEqual(tc, r1.ts.mag_dBi, r2.ts.mag_dBi, "same seed, same start phases");
end

function testFigure(tc)
f = fullfile(tempname, "simple_mag.png");
mkdir(fileparts(f));
simple_mag_nexrad_asr11_rev1(fig_file=f);
verifyTrue(tc, isfile(f));
rmdir(fileparts(f), "s");
end
