function tests = test_macd_monte_carlo_rev1
%TEST_MACD_MONTE_CARLO_REV1 Tests for macd_monte_carlo_rev1.
%   Hand-worked shares for a box pattern, agreement with the exact
%   distribution of simple_mag_nexrad_asr11_rev1, checkpoints, and parallel
%   against serial.
%
%   DRAFT//PRE-DECISIONAL
tests = functiontests(localfunctions);
end

function testBoxSharesByHand(tc)
% Box pattern: 30 dBi within 1 deg of the path. With continuous pointing,
% each antenna is in its main beam 2/360 of the time, so both are
% (1/180)^2 and exactly one is 2 (1/180)(179/180).
a = box_ant();
n = 2e7;
o = macd_monte_carlo_rev1(a, a, n_draws=n, chunk=1e6, bin_dB=0.5);
e = o.exceed;
p60 = e.exceed(find(e.mag_dBi >= 59.9, 1));
p30 = e.exceed(find(e.mag_dBi >= 29.9, 1)) - p60;
verifyEqual(tc, p60, (1/180)^2, "AbsTol", 5*sqrt((1/180)^2/n), "both main beams");
q = 2*(1/180)*(179/180);
verifyEqual(tc, p30, q, "AbsTol", 5*sqrt(q*(1 - q)/n), "one main beam");
end

function testMatchesExactDistribution(tc)
% 1e8 draws: the 99.9% and 99.99% intervals hold the exact values of the
% simple version (45.301 and 50.472 dBi, Section 1.1).
[nex, asr] = section11_ants();
o = macd_monte_carlo_rev1(nex, asr, n_draws=1e8, chunk=1e7, pcts=[99.9 99.99]);
s = simple_mag_nexrad_asr11_rev1();
for p = [99.9 99.99]
    r = o.pct(o.pct.pct == p, :);
    ex = s.pct.mag_dBi(s.pct.pct == p);
    verifyGreaterThanOrEqual(tc, ex, r.ci_lo_dBi - 0.002, sprintf("%g%% lower", p));
    verifyLessThanOrEqual(tc, ex, r.ci_hi_dBi + 0.002, sprintf("%g%% upper", p));
end
end

function testCheckpointsAndModes(tc)
% Checkpoints give one row per draw count and percentage; a two-mode
% antenna (tilt 0 and 10, probabilities 1 and 3) is in its main beam a
% quarter as often.
g = @(az, tilt) 30*(az <= 1 & tilt == 0);
a = struct("g_max_dBi", 30, "gain_fn", g, "backlobe_dBi", 0, ...
    "modes", table([0; 10], [1; 3], 'VariableNames', ["tilt_deg", "prob"]));
b = struct("g_max_dBi", 0, "gain_fn", @(az, t) 0*az, "backlobe_dBi", 0);
n = 4e6;
o = macd_monte_carlo_rev1(a, b, n_draws=n, chunk=1e6, checkpoints=[1e6 2e6], pcts=[50 99], bin_dB=0.5);
verifyEqual(tc, unique(o.pct.n_draws)', [1e6 2e6 4e6]);
verifyEqual(tc, height(o.pct), 6);
e = o.exceed; p = e.exceed(find(e.mag_dBi >= 29.9, 1));
q = 0.25/180;
verifyEqual(tc, p, q, "AbsTol", 5*sqrt(q/n));
end

function testParallelEqualsSerial(tc)
assumeNotEmpty(tc, gcp("nocreate"), "Needs an open parallel pool.");
[nex, asr] = section11_ants();
s = macd_monte_carlo_rev1(nex, asr, n_draws=4e6, chunk=1e6, parallel_flag=0);
p = macd_monte_carlo_rev1(nex, asr, n_draws=4e6, chunk=1e6, parallel_flag=1);
verifyEqual(tc, p.exceed, s.exceed);
verifyEqual(tc, p.pct, s.pct);
end

function a = box_ant()
a = struct("g_max_dBi", 30, "gain_fn", @(az, ~) 30*(az <= 1), "backlobe_dBi", 0);
end

function [nex, asr] = section11_ants()
nex_gain = @(az, tilt) 45.5 + max(-12*(acosd(min(cosd(tilt).*cosd(az), 1))/0.925).^2, ...
    -29*(acosd(min(cosd(tilt).*cosd(az), 1)) <= 10) - 40*(acosd(min(cosd(tilt).*cosd(az), 1)) > 10));
asr_gain = @(az, ~) max(34 - 12*(abs(az)/1.45).^2, 9.5*(abs(az) <= 3.55) + 3.5*(abs(az) > 3.55));
nex = struct("g_max_dBi", 45.5, "gain_fn", nex_gain, "backlobe_dBi", 5.5, ...
    "modes", table(0.5, 1, 'VariableNames', ["tilt_deg", "prob"]));
asr = struct("g_max_dBi", 34, "gain_fn", asr_gain, "backlobe_dBi", 3.5);
end
