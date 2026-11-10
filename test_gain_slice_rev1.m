function tests = test_gain_slice_rev1
%TEST_GAIN_SLICE_REV1 Tests for gain_slice_prep_rev1 and gain_slice_rev1.
%   Reference: the azimuth-grid histogram (antenna_gain_rev1 on 360,000
%   azimuths, gains rounded up to the bin edge), the method of
%   gain_dist_table_rev1 at commit 8c980b7. Every antenna type in
%   antenna_prepare_rev1 is covered, with full rotation and sectors.
%
%   Run: results = runtests("test_gain_slice_rev1");
tests = functiontests(localfunctions);
end

% =========================================================================
function setupOnce(tc)
o = nexrad_asr11_mag_rev1(struct("seed", 1, "n_samples", 1e3, "duration_s", 60, "n_trials", 1));
tc.TestData.ants = struct( ...
    "name", {"WSR-88D NSSL table", "ASR-11 rect csc2", "statgain", "f1245", "parabolic", "rect", "measured"}, ...
    "ant", {o.sys_nexrad.tx, o.sys_asr11.tx, struct("type", "statgain", "g_max_dBi", 40), ...
    struct("type", "f1245", "g_max_dBi", 46), ...
    struct("type", "parabolic", "g_max_dBi", 38, "az_bw_deg", 1.5, "az_sll_dB", -26, "envelope", "peak"), ...
    struct("type", "rect", "g_max_dBi", 33, "az_bw_deg", 1.3, "az_sll_dB", -26, "el_bw_deg", 5, ...
    "el_pattern", "rect", "el_sll_dB", -20, "envelope", "peak"), ...
    struct("type", "measured", "g_max_dBi", 30, "az_deg", [-180 -10 -1 0 1 10 180], ...
    "az_gain_dB", [-40 -30 -3 0 -3 -30 -40], "el_deg", [-90 -5 0 5 90], ...
    "el_gain_dB", [-30 -6 0 -6 -30], "floor_dB", -45)});
tc.TestData.method = ["closed_form", "separable", "closed_form", "closed_form", "closed_form", "separable", "separable"];
tc.TestData.geo = [0.5 0.3; 1.3 1.7; 0 -1; 6.4 2; 19.5 0.5; 0.5 0.5];
tc.TestData.secs = {[], [-30 30], [100 160], [-200 -150], [-5 355]};
tc.TestData.levels = [0.5 1e-1 1e-2 1e-3 1e-4 1e-5];
tc.TestData.dg = 0.05;
end

% =========================================================================
function testMethodChoice(tc)
d = tc.TestData;
for i = 1:numel(d.ants)
    S = gain_slice_prep_rev1(d.ants(i).ant);
    verifyEqual(tc, S.method, d.method(i), d.ants(i).name);
end
end

function testMatchesAzimuthGrid(tc)
% Closed form bounds the pattern from above in each cell and has exact
% time shares, so it may differ from the grid by one bin, never more.
% Separable uses the same azimuth samples, so it must equal the grid.
d = tc.TestData;
for i = 1:numel(d.ants)
    S = gain_slice_prep_rev1(d.ants(i).ant);
    for gi = 1:size(d.geo, 1)
        for si = 1:numel(d.secs)
            lim = d.secs{si};
            h = gain_slice_rev1(S, d.geo(gi, 1), d.geo(gi, 2), lim);
            [kr, pr, gmax_r] = grid_ref(S, d.geo(gi, 1), d.geo(gi, 2), lim, d.dg);
            dq = pct_diff(h, kr, pr, d.levels);
            msg = sprintf("%s, tilt %.1f, el %.1f, sector %s", d.ants(i).name, d.geo(gi, 1), ...
                d.geo(gi, 2), mat2str(lim));
            verifyEqual(tc, sum(h.p), 1, "AbsTol", 1e-12, "mass: " + msg);
            if S.method == "separable"
                verifyEqual(tc, dq, zeros(size(dq)), "AbsTol", 1e-9, "separable vs grid: " + msg);
            else
                verifyLessThanOrEqual(tc, abs(dq), d.dg + 1e-9, "closed form vs grid: " + msg);
            end
            verifyGreaterThanOrEqual(tc, h.gmax, gmax_r - 1e-9, "peak: " + msg);
        end
    end
end
end

function testFullSectorEqualsRotation(tc)
d = tc.TestData;
for i = 1:numel(d.ants)
    S = gain_slice_prep_rev1(d.ants(i).ant);
    a = gain_slice_rev1(S, 0.5, 1.0);
    b = gain_slice_rev1(S, 0.5, 1.0, [-180 180]);
    verifyEqual(tc, b.k_lo, a.k_lo, d.ants(i).name);
    verifyEqual(tc, b.p, a.p, "AbsTol", 1e-12, d.ants(i).name);
end
end

function testSectorAimedAwayLowersPeak(tc)
% A sector that points away from the station cannot reach the main beam.
S = gain_slice_prep_rev1(tc.TestData.ants(1).ant);
on = gain_slice_rev1(S, 0.5, 0.5, [-10 10]);
off = gain_slice_rev1(S, 0.5, 0.5, [90 110]);
verifyGreaterThan(tc, on.gmax - off.gmax, 40);
end

function testSectorInsideOnePatternSegment(tc)
% Red team F1 (step 8): a sector whose off-axis range falls inside one
% segment of a coarse table (no breakpoint inside) used to crash in
% refine_cells. It must match the grid like any other slice.
ant = struct("type", "offaxis_table", "g_max_dBi", 41.1464, ...
    "offaxis_deg", [0 0.446 0.8958 1.9947 2.4744 2.7664 90 180], ...
    "gain_dB", [0 -3 -10.9467 -29.6758 -42.6005 -57.3186 -59.6935 -61.4425]);
S = gain_slice_prep_rev1(ant);
lim = [-181.98, -125.45];
h = gain_slice_rev1(S, 4.6886, 1.2, lim);
[kr, pr] = grid_ref(S, 4.6886, 1.2, lim, tc.TestData.dg);
verifyLessThanOrEqual(tc, abs(pct_diff(h, kr, pr, tc.TestData.levels)), tc.TestData.dg + 1e-9);
verifyEqual(tc, sum(h.p), 1, "AbsTol", 1e-12);
end

function testSeparableSectorPeakIsExact(tc)
% Review R5: for a separable sector slice, gmax was the largest azimuth
% sample, up to 0.007 dB below the true peak. It now includes az_off = 0
% when inside the sector and the sector ends.
S = gain_slice_prep_rev1(tc.TestData.ants(2).ant);
verifyEqual(tc, S.method, "separable");
for c = [0.31 -2.17 2.4]
    lim = c + [-3.0003, 2.9997];               % holds az_off = 0, off the sample grid
    h = gain_slice_rev1(S, 1.0, 0.7, lim);
    st = struct("mech_az", 0, "tilt", 1.0, "steer_az", 0, "steer_el", 0);
    verifyEqual(tc, h.gmax, antenna_gain_rev1(S.pa, 0, 0.7, st), "AbsTol", 1e-12);
end
lim = [10.00037, 30.00037];                      % main beam outside: peak at the near end
h = gain_slice_rev1(S, 1.0, 0.7, lim);
st = struct("mech_az", -10.00037, "tilt", 1.0, "steer_az", 0, "steer_el", 0);
verifyGreaterThanOrEqual(tc, h.gmax, antenna_gain_rev1(S.pa, 0, 0.7, st) - 1e-12);
end

function testZenithIsOneGain(tc)
S = gain_slice_prep_rev1(tc.TestData.ants(1).ant);
h = gain_slice_rev1(S, 0.5, 90);
verifyEqual(tc, h.p, 1);
end

function testGridFallbackForArrays(tc)
% Planar arrays use the grid method: a coarse step keeps the test fast.
S = gain_slice_prep_rev1(struct("type", "array_planar", "g_max_dBi", 38), struct("az_step_deg", 0.05));
verifyEqual(tc, S.method, "grid");
h = gain_slice_rev1(S, 0, 1);
[kr, pr] = grid_ref(S, 0, 1, [], 0.05);
verifyEqual(tc, pct_diff(h, kr, pr, tc.TestData.levels), zeros(1, 6), "AbsTol", 1e-9);
end

function testInputErrors(tc)
S = gain_slice_prep_rev1(tc.TestData.ants(1).ant);
verifyError(tc, @() gain_slice_rev1(S, 0.5, 0, [10 10]), "gain_slice:badLim");
verifyError(tc, @() gain_slice_rev1(S, 0.5, 0, [0 400]), "gain_slice:badLim");
verifyError(tc, @() gain_slice_prep_rev1(tc.TestData.ants(1).ant, struct("bogus", 1)), "gain_slice:unknownOption");
verifyError(tc, @() gain_slice_prep_rev1(tc.TestData.ants(2).ant, struct("method", "closed_form")), "gain_slice:badMethod");
end

% =========================================================================
function [k_lo, p, gmax] = grid_ref(S, tilt, el, lim, dg)
az = S.az;
if ~isempty(lim)
    lo = mod(lim(1) + 180, 360) - 180;
    az = az(mod(az - lo, 360) <= min(diff(lim), 360));
end
z = zeros(size(az));
g = antenna_gain_rev1(S.pa, 0, el, struct("mech_az", -az, "tilt", z + tilt, "steer_az", z, "steer_el", z));
k = ceil(g/dg - 1e-9);
k_lo = min(k);
p = accumarray(k - k_lo + 1, 1)/numel(k);
gmax = max(g);
end

function dq = pct_diff(h, kr, pr, levels)
% Difference (bins x bin width is applied by the caller's tolerance; here
% in dB) of the gain exceeded with probability <= each level.
dg = 0.05;
k0 = min(h.k_lo, kr);
n = max(h.k_lo + numel(h.p), kr + numel(pr)) - k0;
pn = zeros(n, 1);
pg = zeros(n, 1);
pn(h.k_lo - k0 + (1:numel(h.p))) = h.p;
pg(kr - k0 + (1:numel(pr))) = pr;
En = cumsum(pn, "reverse") - pn;
Eg = cumsum(pg, "reverse") - pg;
dq = zeros(1, numel(levels));
for j = 1:numel(levels)
    dq(j) = (find(En <= levels(j) + 1e-13, 1) - find(Eg <= levels(j) + 1e-13, 1))*dg;
end
end
