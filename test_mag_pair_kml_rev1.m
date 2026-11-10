function tests = test_mag_pair_kml_rev1
%TEST_MAG_PAIR_KML_REV1 Tests for mag_pair_kml_rev1 (KMZ and KML output).
%   Run: results = runtests("test_mag_pair_kml_rev1");
tests = functiontests(localfunctions);
end

function setupOnce(tc)
o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1, "cdf_method", "exact"));
g = -0.5:0.05:1;
Tn = gain_dist_table_rev1(o.sys_nexrad.tx, o.sys_nexrad.scan, struct("name", "WSR-88D", "el_grid_deg", g));
Ta = gain_dist_table_rev1(o.sys_asr11.tx, o.sys_asr11.scan, struct("name", "ASR-11", "el_grid_deg", g));
% Test geometry only (not real sites).
pairs = table([35.0; 35.0; 36.0], [-97.0; -97.0; -98.0], [35.2; 34.9; 36.1], [-97.3; -96.8; -98.2], ...
    ["Test N1"; "Test N1"; "Test N2"], ["Test A1"; "Test A2"; "Test A3"], ...
    'VariableNames', ["lat_tx", "lon_tx", "lat_rx", "lon_rx", "name_tx", "name_rx"]);
tc.TestData.R = mag_pair_cdf_rev1(Tn, Ta, [0; 0.3; 0.8], [0.1; -0.2; 0.5]);
tc.TestData.pairs = pairs;
tc.TestData.dir = tempname;
mkdir(tc.TestData.dir);
end

function teardownOnce(tc)
rmdir(tc.TestData.dir, "s");
end

function testKmzHoldsDocFirstAndAllPairs(tc)
d = tc.TestData;
f = mag_pair_kml_rev1(fullfile(d.dir, "pairs.kmz"), d.pairs, d.R);
verifyTrue(tc, isfile(f));
u = fullfile(d.dir, "unz");
names = unzip(f, u);
[~, b, e] = fileparts(string(names));
verifyEqual(tc, b(1) + e(1), "doc.kml");
verifyTrue(tc, all(ismember(["pairs.kml", "sites.kml"], b + e)));
G = readgeotable(fullfile(u, "pairs.kml"));
verifyEqual(tc, height(G), 3);
S = readgeotable(fullfile(u, "sites.kml"));
verifyEqual(tc, height(S), 5);                   % Test N1 appears in two pairs
k = fileread(fullfile(u, "pairs.kml"));
verifyEqual(tc, numel(strfind(k, "DRAFT//PRE-DECISIONAL")), 3);
verifySubstring(tc, k, compose("%.2f", d.R.pct.p99_99(2)));
verifySubstring(tc, fileread(fullfile(u, "doc.kml")), "<href>pairs.kml</href>");
end

function testKmlWritesTwoFiles(tc)
d = tc.TestData;
f = mag_pair_kml_rev1(fullfile(d.dir, "out.kml"), d.pairs, d.R, struct("pct", 99.9));
verifyEqual(tc, numel(f), 2);
verifyTrue(tc, all(isfile(f)));
verifySubstring(tc, fileread(f(1)), "99.9%:");
end

function testInputErrors(tc)
d = tc.TestData;
f = fullfile(d.dir, "x.kmz");
verifyError(tc, @() mag_pair_kml_rev1(f, d.pairs(:, 1:3), d.R), "mag_kml:missingColumns");
verifyError(tc, @() mag_pair_kml_rev1(f, d.pairs(1:2, :), d.R), "mag_kml:sizeMismatch");
verifyError(tc, @() mag_pair_kml_rev1(f, d.pairs, d.R, struct("pct", 12)), "mag_kml:badPct");
verifyError(tc, @() mag_pair_kml_rev1(fullfile(d.dir, "x.txt"), d.pairs, d.R), "mag_kml:badExt");
verifyError(tc, @() mag_pair_kml_rev1(f, d.pairs, d.R, struct("bogus", 1)), "mag_kml:unknownOption");
end
