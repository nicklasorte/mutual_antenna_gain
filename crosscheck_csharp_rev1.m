%CROSSCHECK_CSHARP_REV1 Compare the MATLAB radar coupling models with the C# sources.
%   Ports the C# reference code from the 2.7 GHz project (AntennaModels
%   namespace) into local functions named cs_*, line for line, and compares
%   them with antenna_prepare_rev1, antenna_gain_rev1 and
%   radar_mutual_gain_cdf_rev1. Prints one row per check.
%
%   C# files ported:
%     AntennaModel_StatGain.cs      Statgain()
%     AntennaModel_ITU_1245_1.cs    CalculateAntGain_dBi() (no MATLAB model yet)
%     AperatureAntenna.cs           DetermineGain(), DetermineMaskGain()
%     CosecantSquaredAntenna.cs     DetermineGain4GroundRadar()
%     PlanarArrayAntenna.cs         DetermineArrayGain() (psi_y bug flagged)
%     AntennaGainCalculations.cs    CalculateAnt3dBBeamwidth() search
%     AntennaGainDistributions.cs   DetermineAntennaGainDistribution(), MACD
%     ScanRateAnalysis.cs           worst-case SRA (per paper Appendix C)

clear;
seed = 20261101;
rows = {};

%% 1. Statgain pattern, 10 <= G <= 60 dBi (TM-13-489 range used by MATLAB)
th = (0:0.01:180)';
st = states(numel(th));
G = 10:0.25:60;
d = zeros(size(G));
for k = 1:numel(G)
    p = antenna_prepare_rev1(struct("g_max_dBi", G(k)));
    d(k) = max(abs(antenna_gain_rev1(p, th, 0, st) - cs_statgain(G(k), th)));
end
rows(end+1, :) = {"Statgain pattern vs C# Statgain()", max(d), 1e-9};

%% 2. M.1851-2 rectangular aperture, theoretical and masks (no pedestal)
%     C# AperatureAntenna types Uniform, Cos, Cos2, Cos3, Cos4; MATLAB picks
%     the taper from az_sll_dB with M.1851-2 Table 9 (test values).
cs_types = ["Uniform", "Cos", "Cos2", "Cos3", "Cos4"];
sll = [-15 -25 -35 -42 -47];
bw = 2;
a = (0:0.005:90)';
st = states(numel(a));
for env = ["theoretical", "peak", "average"]
    worst = 0;
    for n = 1:5
        ant = struct("type", "rect", "g_max_dBi", 0, "az_bw_deg", bw, "az_sll_dB", sll(n), ...
            "el_bw_deg", 10, "el_sll_dB", -15, "envelope", env, "floor_dB", -200);
        p = antenna_prepare_rev1(ant);
        g_m = antenna_gain_rev1(p, a, 0, st);
        switch env
            case "theoretical"
                g_c = cs_aperture(cs_types(n), bw, a, false, false);
                ok = g_c > p.az_cut.floor_dB + 1;      % MATLAB floors nulls
            case "peak"
                g_c = cs_aperture(cs_types(n), bw, a, true, true);
                ok = true(size(a));
            case "average"
                g_c = cs_aperture(cs_types(n), bw, a, true, false);
                ok = true(size(a));
        end
        worst = max(worst, max(abs(g_m(ok) - g_c(ok))));
    end
    rows(end+1, :) = {"M.1851 rect aperture (" + env + ") vs C# AperatureAntenna", worst, 0.02}; %#ok<SAGROW>
end

%% 3. M.1851-2 CSC^2 ground radar elevation pattern (Eqs 25-27)
el = (-10:0.005:40)';
st = states(numel(el));
st.tilt = 2 + zeros(size(el));
ant = struct("type", "rect", "g_max_dBi", 0, "az_bw_deg", 1.35, "az_sll_dB", -26.2, ...
    "el_bw_deg", 4.8, "el_pattern", "csc2", "csc2_end_deg", 30, "csc2_floor_dB", -55, ...
    "envelope", "peak", "floor_dB", -200);
p = antenna_prepare_rev1(ant);
g_m = antenna_gain_rev1(p, 0, el, st);
g_c = cs_csc2_ground(0, 4.8, 2, el, 30, -55);
rows(end+1, :) = {"M.1851 CSC^2 elevation vs C# CosecantSquaredAntenna", max(abs(g_m - g_c)), 1e-6};

%% 4. Planar array, uniform, isotropic element: azimuth cut, steered 0 and 30 deg
lam = physconst("LightSpeed")/3e9;
arr = struct("type", "array_planar", "nx", 32, "ny", 16, "dx_m", lam/2, "dy_m", lam/2, ...
    "freq_MHz", 3000, "elem_gain_dBi", 0, "elem_cos_exp", 0, "taper", "uniform", ...
    "floor_dB", -200, "front_back_dB", -200);
p = antenna_prepare_rev1(arr);
az = (-89:0.01:89)';
worst = 0;
for s = [0 30]
    st = states(numel(az));
    st.steer_az = s + zeros(size(az));
    g_m = antenna_gain_rev1(p, az, 0, st) - p.g_max_dBi;
    g_c = cs_planar(32, 0.5, 16, 0.5, s, 0, abs(az), 180*(az < 0), true);
    nrm_err = max(g_c);              % C# normalizes at broadside, not the beam
    g_c = g_c - nrm_err;
    ok = g_c > -40;                  % MATLAB interpolates |AF|^2 near nulls
    worst = max(worst, max(abs(g_m(ok) - g_c(ok))));
end
rows(end+1, :) = {"Planar array azimuth cut shape vs C# PlanarArrayAntenna", worst, 0.05};
rows(end+1, :) = {"INFO: C# steered-beam normalization error at 30 deg steer (dB)", nrm_err, NaN};
% Elevation cut: C# never sets psi_y_coeff (Initialize assigns psi_x_coeff
% twice), so its vertical array factor is flat. Report the error size.
elv = (0:0.05:60)';
st = states(numel(elv));
g_m = antenna_gain_rev1(p, 0, elv, st) - p.g_max_dBi;
g_bug = cs_planar(32, 0.5, 16, 0.5, 0, 0, elv, 90 + zeros(size(elv)), false);
g_fix = cs_planar(32, 0.5, 16, 0.5, 0, 0, elv, 90 + zeros(size(elv)), true);
ok = g_fix > -40;
rows(end+1, :) = {"Planar array elevation cut vs C# (psi_y fixed)", max(abs(g_m(ok) - g_fix(ok))), 0.05};
ok = g_fix > -60;
rows(end+1, :) = {"INFO: C# psi_y bug, max elevation-cut error above -60 dB (dB)", max(g_bug(ok) - g_fix(ok)), NaN};

%% 5. MACD (unsynchronized) vs C# DetermineMutualCouplingDistribution
pct = [50 90 99 99.9];
for pair = [33.5 34; 34 34; 34 46]'
    ref = cs_macd(pair(1), pair(2), pct);
    s1 = struct("tx", struct("g_max_dBi", pair(1)), "scan", struct("type", "mech_rotating", "rate_deg_s", 75));
    s2 = struct("tx", struct("g_max_dBi", pair(2)), "scan", struct("type", "mech_rotating", "rate_deg_s", 90.37));
    o = struct("seed", seed, "mode", "random", "n_samples", 5e6);
    out = radar_mutual_gain_cdf_rev1(s1, s2, o);
    q = prctile(out.g21_dBi, pct);
    cont = continuous_macd(pair(1), pair(2), pct);
    rows(end+1, :) = {sprintf("MACD %g/%g dBi pct %s: max |MATLAB - continuous|", pair, mat2str(pct)), ...
        max(abs(q(:) - cont(:))), 0.1}; %#ok<SAGROW>
    rows(end+1, :) = {sprintf("INFO: MACD %g/%g dBi, C# bw/4 sampling error, worst pct", pair), ...
        max(abs(ref.exact(:) - cont(:))), NaN}; %#ok<SAGROW>
    rows(end+1, :) = {sprintf("INFO: MACD %g/%g dBi, C# 1 dB binning error at 99%%", pair), ...
        ref.binned(3) - ref.exact(3), NaN}; %#ok<SAGROW>
end

%% 6. Synchronized (same-rate) method vs C# swept grid
for pair = [34 34; 34 46]'
    s1 = struct("tx", struct("g_max_dBi", pair(1)), "scan", struct("type", "mech_rotating", "rate_deg_s", 75));
    s2 = struct("tx", struct("g_max_dBi", pair(2)), "scan", struct("type", "mech_rotating", "rate_deg_s", 75));
    o = struct("seed", seed, "mode", "random", "n_samples", 1e5, "scan_sync", "synchronized");
    out = radar_mutual_gain_cdf_rev1(s1, s2, o);
    [grid_on, grid_off] = cs_sync_grid(pair(1), pair(2));
    x = (0:0.001:360)';
    cont = max(cs_statgain(pair(1), abs(wrapTo180(x))) + cs_statgain(pair(2), abs(wrapTo180(x - 180))));
    rows(end+1, :) = {sprintf("Sync %g/%g dBi: MATLAB sweep max - continuous max", pair), ...
        out.sync.max_g21_dBi - cont, 1e-6}; %#ok<SAGROW>
    rows(end+1, :) = {sprintf("INFO: Sync %g/%g dBi, C# bw/4 grid shortfall (worst of 2 phases)", pair), ...
        cont - min(grid_on, grid_off), NaN}; %#ok<SAGROW>
end

%% 7. Scan rate analysis (worst case) vs a brute-force scan-by-scan check
%     Paper Table 4.4-1 rates: ASR victim 12.5 rpm, WSR-88D interferer 3 rpm.
[sra, brute] = cs_sra_check(34, 12.5, 46, 3, 60);
rows(end+1, :) = {"SRA consecutive misses: C# formula - brute force", sra.n_consec - brute.n_consec, 0};
rows(end+1, :) = {"SRA scans between misses: C# - brute force", sra.n_between - brute.n_between, 0};

%% 8. ITU-R F.1245 (paper Section 4.4 model) and the paper MAG values
th = (0:0.01:180)';
G = 10:0.25:60;                      % covers both D/lambda branches
d = zeros(size(G));
step48 = zeros(size(G));
off48 = abs(th - 48) > 1e-6;         % F.1245 is discontinuous at 48 deg
for k = 1:numel(G)
    p = antenna_prepare_rev1(struct("type", "f1245", "g_max_dBi", G(k)));
    d(k) = max(abs(antenna_gain_rev1(p, th(off48), 0, states(sum(off48))) - cs_f1245(G(k), th(off48))));
    step48(k) = abs(diff(cs_f1245(G(k), [48 - 1e-9; 48])));
end
rows(end+1, :) = {"F.1245 pattern vs C# AntennaModel_ITU_1245_1 (48 deg excluded)", max(d), 1e-9};
rows(end+1, :) = {"INFO: F.1245 step at 48 deg, largest over 10-60 dBi (dB)", max(step48), NaN};
for pair = [34 34 24.425; 34 46 36.425; 46 46 33.425]'
    s1 = struct("tx", struct("type", "f1245", "g_max_dBi", pair(1)), "scan", struct("type", "mech_rotating", "rate_deg_s", 75));
    s2 = struct("tx", struct("type", "f1245", "g_max_dBi", pair(2)), "scan", struct("type", "mech_rotating", "rate_deg_s", 75));
    o = struct("seed", seed, "mode", "random", "n_samples", 1e5, "scan_sync", "synchronized");
    out = radar_mutual_gain_cdf_rev1(s1, s2, o);
    rows(end+1, :) = {sprintf("F.1245 sync %g/%g dBi - paper Table 4.6-1 MAG", pair(1:2)), ...
        out.sync.max_g21_dBi - pair(3), 1e-3}; %#ok<SAGROW>
end
rows(end+1, :) = {"INFO: F.1245 34/34 99% MACD (C# bw/4 sampling)", cs_macd(34, 34, 99, @cs_f1245).exact, NaN};

%% Report
T = cell2table(rows, "VariableNames", ["check", "value", "tol"]);
T.pass = strings(height(T), 1);
is_chk = ~isnan(T.tol);
T.pass(is_chk) = string(abs(T.value(is_chk)) <= T.tol(is_chk));
T.pass(~is_chk) = "info";
disp(T)

% =========================================================================
function st = states(n)
z = zeros(n, 1);
st = struct("mech_az", z, "tilt", z, "steer_az", z, "steer_el", z);
end

% ---- AntennaModel_StatGain.cs, Statgain() -------------------------------
function g = cs_statgain(G, th)
gm = 10^(G/10);
g = zeros(size(th));
if G >= 48
    tm = 50*sqrt(0.25*G + 7)/sqrt(gm); tr = 27.466*gm^-0.3; tb = 48;
    side = @(t) 29 - 25*log10(t); back = -13;
elseif G >= 22
    tm = 50*sqrt(0.25*G + 7)/sqrt(gm); tr = 250/sqrt(gm); tb = 48;
    side = @(t) 53 - G/2 - 25*log10(t); back = 11 - G/2;
else
    tm = 50*sqrt(0.25*G + 7)/sqrt(gm); tr = 250/sqrt(gm); tb = 131.8257*gm^-0.2;
    side = @(t) 53 - G/2 - 25*log10(t); back = 0;
end
i1 = th <= tm; i2 = th > tm & th <= tr; i3 = th > tr & th <= tb; i4 = th > tb;
g(i1) = G - 4e-4*gm*th(i1).^2;
g(i2) = 0.75*G - 7;
g(i3) = side(th(i3));
g(i4) = back;
end

% ---- AntennaModel_ITU_1245_1.cs -----------------------------------------
function g = cs_f1245(G, th)
dl = 10^((G - 7.7)/20);
G1 = 2 + 15*log10(dl);
pm = 20*sqrt(G - G1)/dl;
pr = 12.02*dl^-0.6;
pmax = max(pm, pr);
th = abs(th);
g = zeros(size(th));
main = th < pm;
g(main) = G - 0.0025*(dl*th(main)).^2;
if dl > 100
    i = ~main & th < pmax; g(i) = G1;
    i = th >= pmax & th < 48; g(i) = 29 - 25*log10(th(i));
    g(th >= 48) = -13;
else
    i = ~main & th < 48; g(i) = 39 - 5*log10(dl) - 25*log10(th(i));
    g(th >= 48) = -3 - 5*log10(dl);
end
end

% ---- AperatureAntenna.cs --------------------------------------------------
function g = cs_aperture(type, theta3, th, for_mask, for_peak)
switch type
    case "Uniform", K = 50.8; A = 8.584; B = 2.876; fl = -30; pk = -5.75; av = -12.16; p2a = -3.72; nrm = 0;
    case "Cos",     K = 68.8; A = 17.51; B = 2.33;  fl = -50; pk = -14.4; av = -20.6;  p2a = -4.32; nrm = 3.92;
    case "Cos2",    K = 83.2; A = 26.882; B = 1.962; fl = -60; pk = -22.3; av = -29.0; p2a = -4.6;  nrm = 6.02;
    case "Cos3",    K = 95.0; A = 35.84; B = 1.756; fl = -70; pk = -31.5; av = -37.6;  p2a = -4.2;  nrm = 7.44;
    case "Cos4",    K = 106.0; A = 45.88; B = 1.56; fl = -80; pk = -39.4; av = -42.5;  p2a = -2.6;  nrm = 8.52;
end
gain = @(t) 20*log10(abs(cs_fmu(type, pi*K*sind(t)/theta3))) + nrm;
% First null: 0.1 deg search, as in SetMaskParameters().
t = (0:0.1:89.9)';
gt = gain(t);
k = find(gt(2:end-1) < gt(1:end-2) & gt(2:end-1) < gt(3:end), 1) + 1;
first_null = t(k + 1);           % C# stores offset_deg of the next sample
g = gain(th);
if ~for_mask
    return
end
m = -A*log(B*abs(th)/theta3);
m(th == 0) = g(th == 0);
if for_peak
    keep = g >= pk & th < first_null;
else
    keep = g >= av & th < first_null;
    m(~keep) = m(~keep) + p2a;
end
m(~keep) = max(m(~keep), fl);
m(keep) = g(keep);
g = m;
end

function f = cs_fmu(type, mu)
switch type
    case "Uniform"
        f = sin(mu)./mu;
    case "Cos"
        f = (pi/2)*cos(mu)./((pi/2)^2 - mu.^2);
    case "Cos2"
        f = (pi^2/2)*(sin(mu)./mu)./(pi^2 - mu.^2);
    case "Cos3"
        f = (3*pi*cos(mu)/8).*(1./((pi/2)^2 - mu.^2) - 1./(9*(pi/2)^2 - mu.^2));
    case "Cos4"
        f = 3*pi^4*sin(mu)./(2*mu.*(mu.^2 - pi^2).*(mu.^2 - 4*pi^2));
end
% C# MathFunctions.Sinc(0) and the closed forms at mu = 0: take the limit.
lim = struct("Uniform", 1, "Cos", 2/pi, "Cos2", 0.5, "Cos3", 4/(3*pi), "Cos4", 3/8);
f(mu == 0) = lim.(type);
end

% ---- CosecantSquaredAntenna.cs, ground radar ------------------------------
function g = cs_csc2_ground(gmax, theta3, tilt, th, th_end, floor_dB)
th_null = tilt - theta3/0.88;
th_start = theta3/2 + tilt;
sincx = @(x) sin(x)./x;
g = floor_dB + zeros(size(th));
i = th > th_null & th < th_start;
mu = 50.8*pi*sind(th(i) - tilt)/theta3;
s = sincx(mu); s(mu == 0) = 1;
g(i) = 20*log10(abs(s));
i = th >= th_start & th <= th_end;
arg = 50.8*pi*sind(th_start - tilt)/theta3;
g(i) = 20*log10(abs((1./sind(th(i)))/(1/sind(th_start)))) + 20*log10(abs(sincx(arg)));
g = g + gmax;
% MATLAB applies the floor to the whole cut (M.1851-2 Eq. 25 floor G0).
g = max(g, floor_dB + gmax);
end

% ---- PlanarArrayAntenna.cs -------------------------------------------------
function g = cs_planar(nx, dxr, ny, dyr, oth, ophi, th, phi, fix_psi_y)
cx = 2*pi*dxr;
cy = 2*pi*dyr;
if ~fix_psi_y
    cy = 0;                      % C# Initialize() never sets psi_y_coeff
end
af = @(t, p) afn(nx, cx*(sind(t).*cosd(p) - sind(oth)*cosd(ophi))) ...
    .* afn(ny, cy*(sind(t).*sind(p) - sind(oth)*sind(ophi)));
nrm = 10*log10(af(0, 0).^2/(nx*ny));
g = 10*log10(af(th, phi).^2/(nx*ny)) - nrm;
end

function a = afn(n, psi)
a = sin(n*psi/2)./sin(psi/2);
a(abs(sin(psi/2)) < 1e-12) = n;
end

% ---- AntennaGainCalculations.cs, 3 dB beamwidth search ---------------------
function bw = cs_beamwidth(G, model)
t = (0.001:0.001:180)';
[~, i] = min(abs(model(G, t) - (G - 3)));
bw = 2*t(i);
end

% ---- AntennaGainDistributions.cs, MACD --------------------------------------
function r = cs_macd(g1, g2, pct, model)
if nargin < 4
    model = @cs_statgain;
end
[x1, lo1] = gain_sweep(g1, 0, model);       % antenna 1, bearing to 2 = 0
[x2, lo2] = gain_sweep(g2, 180, model);     % antenna 2, bearing to 1 = 180
s = sort(reshape(x1 + x2', [], 1));         % convolution of uniform weights
r.exact = s(max(1, ceil(pct/100*numel(s))));
% 1 dB histograms (bin lower edge = (int)backlobe + k), then convolution.
b1 = lo1 + floor(x1 - lo1);
b2 = lo2 + floor(x2 - lo2);
sb = sort(reshape(b1 + b2', [], 1));
r.binned = sb(max(1, ceil(pct/100*numel(sb))));
end

function q = continuous_macd(g1, g2, pct)
% Independent uniform pointing on a 0.01 deg grid (near-continuous).
% Each gain distribution is histogrammed on a 0.01 dB grid and the two are
% convolved, so percentiles carry at most ~0.01 dB of grid error.
t = (0:0.001:180)';
[ha, a0] = fine_hist(cs_statgain(g1, t));
[hb, b0] = fine_hist(cs_statgain(g2, t));
F = cumsum(conv(ha, hb));
x = a0 + b0 + (0:numel(F)-1)'*0.01;
q = zeros(size(pct));
for k = 1:numel(pct)
    q(k) = x(find(F >= pct(k)/100, 1));
end
end

function [h, c0] = fine_hist(g)
c0 = floor(min(g)*100)/100;
idx = round((g - c0)/0.01) + 1;
h = accumarray(idx, 1)/numel(g);
end

function [g, lo] = gain_sweep(G, bear_other, model)
bw = cs_beamwidth(G, model);
inc = bw/4;
n = floor(360/inc) + 1;
brg = -180 + (0:n-1)'*inc;
rel = abs(wrapTo180(bear_other - brg));
g = model(G, rel);
lo = fix(model(G, 180));
end

% ---- Synchronized method, paper Section 3.1 --------------------------------
function [gmax_on, gmax_off] = cs_sync_grid(g1, g2)
% Both antennas start at one angle and step together by bw/4 of the
% narrower beam. With the link bearing on the grid the grid max is exact;
% half a step off the grid is the worst case.
inc = min(cs_beamwidth(g1, @cs_statgain), cs_beamwidth(g2, @cs_statgain))/4;
x = (0:inc:360)';
f = @(off) max(cs_statgain(g1, abs(wrapTo180(x + off))) + ...
    cs_statgain(g2, abs(wrapTo180(x + off - 180))));
gmax_on = f(0);
gmax_off = f(inc/2);
end

% ---- ScanRateAnalysis.cs, worst case ----------------------------------------
function [sra, brute] = cs_sra_check(g_v, rpm_v, g_i, rpm_i, mag_thresh_dBi)
% Victim boresight on the interferer every victim scan (worst case).
tv = 60/rpm_v;
dth = wrapTo180(6*rpm_i*tv);                       % Eq. C-4, normalized
g_thr = mag_thresh_dBi - g_v;                      % Eq. C-2
t = (0.001:0.001:180)';
gi = cs_statgain(g_i, t);
[~, k] = min(abs(gi - g_thr));
rqd = t(k);                                        % CalculateOffset4GivenGain
sra.n_consec = fix(2*rqd/abs(dth)) + 1;            % Eq. C-9
n = 1;
while abs(wrapTo180(n*dth)) >= rqd && n < 1e5      % search loop, Section C.4
    n = n + 1;
end
sra.n_between = n;
% Brute force: interferer pointing offset at each victim look, every start
% phase; a look is missed when interferer gain >= g_thr.
K = 2000;
phase = (-rqd:rqd/200:rqd)';
runs = zeros(size(phase));
for j = 1:numel(phase)
    miss = cs_statgain(g_i, abs(wrapTo180(phase(j) + (0:K)'*dth))) >= g_thr;
    d = diff([false; miss; false]);
    runs(j) = max(find(d == -1) - find(d == 1));
end
brute.n_consec = max(runs);
miss0 = cs_statgain(g_i, abs(wrapTo180((0:K)'*dth))) >= g_thr;
brute.n_between = find(miss0(2:end), 1);
end
