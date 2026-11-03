%EXAMPLE_RADAR_MUTUAL_GAIN_REV1 Mutual antenna gain CDFs, 2.7-3.1 GHz radars.
%   Three cases against an ATC reflector radar (victim, Radar A):
%     1. Radar B, reflector, rotating at a different rate
%     2. Rotating active phased array, random beam steering in elevation
%     3. Four-face staring phased array, beam steered over 360 deg
%   Antenna and scan values come from Rec. ITU-R M.1464-2 Tables 1 and 2.
%   Reflector radars use the default Statgain pattern (NTIA TM-13-489);
%   arrays use Rec. ITU-R M.1851-2 Section 6.2. Values the recommendations
%   don't give are marked TODO and must be replaced with system data before use.
%   Geometry is location-agnostic: the radars see each other along a
%   horizontal line-of-sight path (opts.el12_deg = opts.el21_deg = 0).
%   Polarization mismatch is not modeled.

clear; close all;

% ---- Reflector radars use the default Statgain pattern (NTIA TM-13-489
%      Table 4-1), which needs only the peak gain. For an M.1851-2 fan-beam
%      or CSC^2 pattern instead, set .type = "rect" and its fields; see
%      antenna_prepare_rev1.

% ---- Radar A, M.1464-2 Table 1: 33.5 dBi, 75 deg/s ----
radar_a = struct("g_max_dBi", 33.5);            % type omitted -> "statgain"

sys1.name = "Radar A (victim)";
sys1.tx = radar_a;
sys1.scan = struct("type", "mech_rotating", "rate_deg_s", 75, ...
    "tilt_deg", 0, ...           % TODO: beam-axis tilt
    "rate_jitter_deg_s", 0.5);   % M.2069-0 Annex 1: +/-0.5 deg/s drive accuracy

% ---- Case 1: Radar B, M.1464-2 Table 1: 34 dBi, 90 deg/s ----
radar_b = struct("g_max_dBi", 34);
sys_b.name = "Radar B";
sys_b.tx = radar_b;
sys_b.scan = struct("type", "mech_rotating", "rate_deg_s", 90, "tilt_deg", 0, ...
    "rate_jitter_deg_s", 0.5);

% ---- Case 2: rotating AESA (M.1464-2 Table 2, Radar L/M class: active
%      array, up to 40-43 dBi, 30-360 deg/s, elevation steering 0-90) ----
%      Array defaults (antenna_prepare_rev1): lambda/2 at 2900 MHz, unit-cell
%      element, Taylor -25 dB, TM-13-489 floor, sized to g_max_dBi.
aesa = struct("type", "array_planar", "g_max_dBi", 40);   % M.1464-2 Radar L/M: up to 40-43 dBi
aesa.face_tilt_deg = 10;                         % TODO: face back-tilt
aesa_tx = aesa;
aesa_tx.taylor_sll_dB = -17;                     % M.1464-2 Radar F1: 17 dB on transmit
sys_aesa.name = "Rotating AESA";
sys_aesa.tx = aesa_tx;
sys_aesa.rx = aesa;
sys_aesa.scan = struct("type", "mech_rotating", "rate_deg_s", 60, ...  % M.1464-2 Radar M: 30-360
    "beam_mode", "random", "beam_az_lim_deg", [0 0], ...
    "beam_el_lim_deg", [0 30], ...               % M.1464-2 Radar K: 0 to 30 deg
    "dwell_s", 0.05);                            % M.1464-2 2.3: tens to hundreds of ms

% ---- Case 3: four-face staring array (M.1464-2 Radar F1: 41 dBi,
%      4 faces, irregular 360 deg coverage, 17/25 dB tx/rx sidelobes) ----
mf = struct("type", "array_planar");             % defaults: 41 dBi, horizon
mf_tx = mf;  mf_tx.taylor_sll_dB = -17;
sys_mf.name = "4-face array";
sys_mf.tx = mf_tx;
sys_mf.rx = mf;
sys_mf.scan = struct("type", "array_multiface", ...  % defaults: 4 faces, random beam,
    "face_az0_deg", 45, ...                      % horizon only, 50 ms dwell
    "beam_el_lim_deg", [0 30]);                  % TODO: beam schedule

% ---- Run ----
opts = struct("seed", 20261101, "mode", "random", "n_samples", 5e6);
cases = {sys_b, sys_aesa, sys_mf};
res = cell(size(cases));
for i = 1:numel(cases)
    res{i} = radar_mutual_gain_cdf_rev1(sys1, cases{i}, opts);
    fprintf("\n%s -> %s\n", cases{i}.name, sys1.name);
    disp(res{i}.pct21)       % path 2 -> 1: interferer tx into victim rx
end

% Time series for Case 1: event durations and the M.2069-0 Fig. 1 view.
opts_ts = struct("seed", 20261101, "mode", "time_series", ...
    "duration_s", 3600, "n_trials", 8, "event_thresh_dBi", [20 30 40 50]);
ts1 = radar_mutual_gain_cdf_rev1(sys1, sys_b, opts_ts);
fprintf("\nRadar B -> Radar A, time series (beat %.0f s; nominal rates repeat every %.0f s)\n", ...
    ts1.geom.beat_period_s, ts1.geom.joint_period_s);
disp(ts1.pct21)
disp(ts1.events12)

% ---- Equal rates: synchronized versus unsynchronized ASR-9 pair ----
%      2.7 GHz link budget working paper Table 4.4-1: ASR-9, 34 dBi,
%      12.5 rpm. Unsynchronized = MACD (unknown relative phase);
%      synchronized = both antennas start at the same angle (offset 0).
asr9 = struct("name", "ASR-9", "tx", struct("g_max_dBi", 34), ...
    "scan", struct("type", "mech_rotating", "rate_rpm", 12.5));
o_un = struct("seed", 20261101, "mode", "random", "n_samples", 5e6);
w = warning("off", "radar_coupling:equalRates");
un = radar_mutual_gain_cdf_rev1(asr9, asr9, o_un);
warning(w);
o_sy = o_un;
o_sy.scan_sync = "synchronized";
o_sy.sync_offset_deg = 0;
sy = radar_mutual_gain_cdf_rev1(asr9, asr9, o_sy);
fprintf("\nASR-9 pair, 12.5 rpm, Statgain: unsynchronized 99%% = %.2f dBi, 99.9%% = %.2f dBi; synchronized max = %.2f dBi\n", ...
    un.pct12.gain_dBi(un.pct12.pct == 99), un.pct12.gain_dBi(un.pct12.pct == 99.9), sy.sync.max_g12_dBi);

% Same pair with the working paper's antenna model (Table 4.4-1: ITU-R
% F.1245). The synchronized max reproduces the paper's 24.425 dBi MAG.
asr9.tx.type = "f1245";
w = warning("off", "radar_coupling:equalRates");
un = radar_mutual_gain_cdf_rev1(asr9, asr9, o_un);
warning(w);
sy = radar_mutual_gain_cdf_rev1(asr9, asr9, o_sy);
fprintf("ASR-9 pair, 12.5 rpm, F.1245: unsynchronized 99%% = %.2f dBi, 99.9%% = %.2f dBi; synchronized max = %.3f dBi\n", ...
    un.pct12.gain_dBi(un.pct12.pct == 99), un.pct12.gain_dBi(un.pct12.pct == 99.9), sy.sync.max_g12_dBi);

% ---- Plots: % of time coupling is exceeded (M.2069-0 Fig. 10 layout) ----
figure("Name", "Mutual gain exceedance");
hold on; grid on;
names = cellfun(@(s) s.name, cases, "UniformOutput", false);
for i = 1:numel(res)
    semilogy(res{i}.cdf21.gain_dBi, res{i}.cdf21.exceed_pct, "LineWidth", 1.5);
end
semilogy(ts1.cdf21.gain_dBi, ts1.cdf21.exceed_pct, "--", "LineWidth", 1.2);
set(gca, "YScale", "log");
ylim([1e-4 100]);
xlabel("Mutual antenna gain (dBi)");
ylabel("Time exceeded (%)");
legend([names, {"Radar B, time series"}], "Location", "southwest");
title("Interferer tx + victim rx gain, 2.7-3.1 GHz");

figure("Name", "Coupling time series");
idx = ts1.ts.t_s <= 600;
plot(ts1.ts.t_s(idx), ts1.ts.g21_dBi(idx));
grid on;
xlabel("Time (s)");
ylabel("Mutual antenna gain (dBi)");
title(sprintf("Radar B to Radar A, rates 90 and 75 deg/s (beat %.0f s)", ts1.geom.beat_period_s));
