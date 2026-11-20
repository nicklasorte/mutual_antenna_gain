function tests = test_macd_rev1
%TEST_MACD_REV1 Tests for macd_rev1, the MATLAB port of the C# MACD.
%   Hand-worked histograms, convolution, and percentage points, and a check
%   that fine settings reproduce the exact distribution of
%   simple_mag_nexrad_asr11_rev1.
%
%   DRAFT//PRE-DECISIONAL
tests = functiontests(localfunctions);
end

function testSamplingMatchesCsharp(tc)
% numBearings = (int)(scanwidth / (BW/4)) + 1, start at center - width/2.
a = box_ant(30, 2, 0);
o = macd_rev1(a, a);
verifyEqual(tc, o.n_samples, [721 721], "360 / 0.5 + 1");
verifyEqual(tc, o.samples{1}.rel_az_deg(1), 180, "first pointing at -180, normalized to 180");
verifyEqual(tc, diff(o.samples{1}.rel_az_deg(2:3)), 0.5, "AbsTol", 1e-12);
end

function testHistogramAndConvolutionByHand(tc)
% Box pattern: 30 dBi within 1 deg of the path, 0 dBi elsewhere. At a
% 0.5 deg step the pointings within 1 deg are -1, -0.5, 0, 0.5, 1: 5 of 721.
a = box_ant(30, 2, 0);
o = macd_rev1(a, a);
h = o.hist{1};
verifyEqual(tc, h.prob(h.bin_lo_dBi == 30), 5/721, "AbsTol", 1e-15);
verifyEqual(tc, h.prob(h.bin_lo_dBi == 0), 716/721, "AbsTol", 1e-15);
% Antenna 2's bearing to antenna 1 is 180 deg (GetParameters), and its
% scan starts and ends at 180 deg, so the repeated end pointing
% (numBearings = 360/step + 1) lands on the path: 6 of 721 in the main beam.
h2 = o.hist{2};
verifyEqual(tc, h2.prob(h2.bin_lo_dBi == 30), 6/721, "AbsTol", 1e-15, "repeated end pointing on the path");
m = o.mag;
verifyEqual(tc, m.prob(m.bin_lo_dBi == 60), 5*6/721^2, "AbsTol", 1e-15, "both main beams");
verifyEqual(tc, m.prob(m.bin_lo_dBi == 30), (5*715 + 716*6)/721^2, "AbsTol", 1e-15, "one main beam");
verifyEqual(tc, sum(m.prob), 1, "AbsTol", 1e-12);
end

function testPercentagePoints(tc)
a = box_ant(30, 2, 0);
o = macd_rev1(a, a, pcts=[50 99.9 99.999]);
% P(MAG < 30) = 716*715/721^2 = 0.98477, so 50% falls in the 0 dBi bin
% and 99.9% in the 30 dBi bin (P(MAG < 60) = 1 - 30/721^2 = 0.999942).
verifyEqual(tc, o.pct.upper_edge_dBi', [1 31 61]);
q = 0.999; below = 716*715/721^2; pk = (5*715 + 716*6)/721^2;
verifyEqual(tc, o.pct.interp_dBi(2), 30 + (q - below)/pk, "AbsTol", 1e-12);
end

function testModeAveraging(tc)
% Two modes, tilt 0 (box beam on the path) and tilt 10 (beam never on the
% path), probabilities 0.25 and 0.75: the 30 dBi share is a quarter of 5/721.
g = @(az, tilt) 30*(az <= 1 & tilt == 0);
a = struct("g_max_dBi", 30, "bw_deg", 2, "gain_fn", g, "backlobe_dBi", 0, ...
    "modes", table([0; 10], [1; 3], 'VariableNames', ["tilt_deg", "prob"]));
o = macd_rev1(a, box_ant(30, 2, 0));
h = o.hist{1};
verifyEqual(tc, h.prob(h.bin_lo_dBi == 30), 0.25*5/721, "AbsTol", 1e-15);
end

function testBinFloorFromTruncatedBacklobe(tc)
% (int)(-5.5) = -5 in C#: the lowest bin starts at -5, and a -5.5 dBi
% backlobe gain is clamped into it.
a = box_ant(20, 2, -5.5);
o = macd_rev1(a, a);
verifyEqual(tc, o.hist{1}.bin_lo_dBi(1), -5);
verifyEqual(tc, o.hist{1}.prob(1), 716/721, "AbsTol", 1e-15);
end

function testFineSettingsMatchExact(tc)
% With a fine step and 0.001 dB bins, the port reproduces the exact
% percentiles of the simple version (which convolves on a 0.0001 deg grid).
nex_gain = @(az, tilt) 45.5 + max(-12*(acosd(min(cosd(tilt).*cosd(az), 1))/0.925).^2, ...
    -29*(acosd(min(cosd(tilt).*cosd(az), 1)) <= 10) - 40*(acosd(min(cosd(tilt).*cosd(az), 1)) > 10));
asr_gain = @(az, ~) max(34 - 12*(abs(az)/1.45).^2, 9.5*(abs(az) <= 3.55) + 3.5*(abs(az) > 3.55));
nex = struct("g_max_dBi", 45.5, "bw_deg", 0.925, "gain_fn", nex_gain, "backlobe_dBi", 5.5, ...
    "modes", table(0.5, 1, 'VariableNames', ["tilt_deg", "prob"]));
asr = struct("g_max_dBi", 34, "bw_deg", 1.45, "gain_fn", asr_gain, "backlobe_dBi", 3.5);
o = macd_rev1(nex, asr, step_div=400, bin_dB=0.001, pcts=[50 99 99.9 99.99 99.999]);
s = simple_mag_nexrad_asr11_rev1();
ref = s.pct.mag_dBi(ismember(s.pct.pct, [50 99 99.9 99.99 99.999]));
verifyEqual(tc, o.pct.interp_dBi, ref, "AbsTol", 0.05);
end

function a = box_ant(g, bw, bl)
a = struct("name", "box", "g_max_dBi", g, "bw_deg", bw, ...
    "gain_fn", @(az, ~) g*(az <= bw/2 + 1e-9) + bl*(az > bw/2 + 1e-9), "backlobe_dBi", bl);
end
